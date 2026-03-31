##########################################################################################
## SETUP ##
##########################################################################################

library(data.table)
library(qs2)
library(here)
library(ggplot2)
library(ggthemes)
library(usethis, include.only = "use_data")
library(collapse, include.only = c("descr", "GRP"))

source(here::here("inst", "00_constants.r"))

theme_set(
  theme_pander(
    base_size = 16, boxes = TRUE, nomargin = FALSE
  ) +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.9, 0.9),
      axis.title.x = element_text(margin = margin(t = 0.5, unit = "cm")),
      axis.title.y = element_text(margin = margin(r = 0.5, unit = "cm"))
    )
)


##########################################################################################
## LOAD ADNI DATA ##
##########################################################################################

## Access to ADNI data requires a Memorandum of Understanding.
## See https://adni.loni.usc.edu for details.

## ADNI data are stored in a private directory outside of the package. The
## path to the directory is stored in the environment variable ADNI_PATH.

## NOTE: The call to unique() when reading in the Berkeley data drops 6
##       rows that appear to be duplicates.
ADNI_PATH <- Sys.getenv("ADNI_PATH")
berk <- unique(fread(file.path(ADNI_PATH, qp("UCBERKELEY_AMY_6MM"))))[SITEID != 381]
names(berk) <- tolower(names(berk))
setkeyv(berk, c("rid", "scandate"))

ptdemog <- fread(file.path(ADNI_PATH, qp("PTDEMOG")))[!is.na(VISDATE) & SITEID != 381]
names(ptdemog) <- tolower(names(ptdemog))
setkeyv(ptdemog, c("rid", "visdate"))

dxsum <- fread(file.path(ADNI_PATH, qp("DXSUM")))[!is.na(EXAMDATE)]
names(dxsum) <- tolower(names(dxsum))
setkeyv(dxsum, c("rid", "examdate"))

## subset to columns needed for estimation
berk <- berk[, .(rid, scandate, centiloids, tracer)]
ptdemog <- ptdemog[, .(rid, visdate, ptgender,
                       ptdob = as.IDate(paste0(ptdob, "/01"), format = "%m/%Y/%d"))]
# with this column subset, need to drop duplicate rows
dxsum <- unique(dxsum[, .(rid, examdate, diagnosis)]) 


##########################################################################################
## SIMULATION PARAMETER: PROPORTION with N SCANS ##
##########################################################################################

## NOTE: Calculate among everyone with a scandate, regardless of missing age or
## centiloid information.

num_scan_props <-
  berk[, .(num_scans = .N), rid] |>
  _[, .N, keyby = num_scans] |>
  _[, P := N / sum(N)]

num_scan_props[]


##########################################################################################
## SIMULATION PARAMETER: TIME BETWEEN SCANS ##
##########################################################################################

## NOTE: Calculate among everyone with a scandate, regardless of missing age or
## centiloid information.

setkeyv(berk, c("rid", "scandate"))
berk[, scan_num := rowid(rid)]
berk[, days_since_last_scan := scandate - shift(scandate, n = 1, type = "lag"),
     keyby = rid]

## Quick summary
descr(berk[scan_num > 1, .(lag_days = days_since_last_scan,
                           lag_years = days_since_last_scan / 365.25)],
      by = GRP(berk[scan_num > 1], ~ scan_num),
      Qprobs = NULL, Ndistinct = FALSE)

## Subset to scans after the first for each patient
berk_multi_scan <- berk[scan_num > 1]

## Inspect time since last scan (in days) by scan number.
## Decision to use a single vector of scan lags is based on visual inspection,
## as these data are simply used to generate plausible amyloid curves.
## NOTE: With the 27Oct2025 download of ADNI, this plot throws a (benign) warning
## because only one participant had scan_num=8.
plot_multiscan <- berk_multi_scan |>
  ggplot(aes(x = days_since_last_scan / 365.25)) +
  geom_density(aes(color = scan_num, group = scan_num)) +
  scale_color_continuous(name = "Scan number") +
  theme_pander()

scan_lag_days <- berk_multi_scan[, density(days_since_last_scan, bw = "SJ")]


##########################################################################################
## CREATE EMPIRICAL ANALYSIS SAMPLE ##
##########################################################################################

# Get first rows from diagnosis and demographic information and merge into the scan data
berkadni <- dxsum[, .SD[1], keyby = rid
                 ][ptdemog[, .SD[1], keyby = rid], on = .(rid)
                 ][berk, on = .(rid)]
setnames(berkadni,
         c("visdate", "examdate", "diagnosis"),
         c("visdate_bl", "examdate_bl", "dx_bl"))

berkadni[, `:=`(age_visdate_bl = (visdate_bl - ptdob) / 365.25,
                age_examdate_bl = (examdate_bl - ptdob) / 365.25,
                age_at_scan = (scandate - ptdob) / 365.25)]

NROW_BERKADNI_INIT       <- berkadni[, .N]
NRID_BERKADNI_INIT       <- berkadni[, uniqueN(rid)]
NSCAN_MISS_CL            <- berkadni[is.na(centiloids), .N]
NRID_MISS_CL             <- berkadni[is.na(centiloids), uniqueN(rid)]

RID_DROP <- berkadni[!is.na(centiloids), .(
  first_nomiss_cent = min(scan_num),
  last_nomiss_cent = max(scan_num)
), rid][first_nomiss_cent > 1 & first_nomiss_cent == last_nomiss_cent][, rid]
NRID_MISS_ONESCAN_NOMISS_CL <- length(RID_DROP)

droptable <- data.table(
  reason  = c("missing centiloids",
              "missing centiloids",
              "missing centiloids + single subsequent non-missing centiloids"),
  rowtype = c("subid", "scan", "subid"),
  value   = c(NRID_MISS_CL, NSCAN_MISS_CL, NRID_MISS_ONESCAN_NOMISS_CL)
)


##########################################################################################
## SIMULATION PARAMETER: AGE AT FIRST SCAN ##
##########################################################################################

berkadni_first_scan <- berkadni[, .SD[1], keyby = rid]

## inspect symmetry of age distribution
berkadni_first_scan[, hist(age_at_scan)]

agedist_first_scan <- berkadni_first_scan[, .(
  mean = mean(age_at_scan),
  sd = sd(age_at_scan)
)]

agedist_first_scan[]


##########################################################################################
## MINIMUM CENTILOIDS ##
##########################################################################################

min_centiloids <- min(berkadni$centiloids, na.rm = TRUE)


##########################################################################################
## COPY OBJECTS TO DATA DIRECTORY ##
##########################################################################################

usethis::use_data(
  num_scan_props,
  scan_lag_days,
  agedist_first_scan,
  min_centiloids,
  droptable,
  overwrite = TRUE,
  internal = FALSE
)


##########################################################################################
## ADD VARIABLES to EMPIRICAL DATASET ##
##########################################################################################

# calculate years since first visit (vidate_bl from ptdemog)
berkadni[, yrs_since_bl := (scandate - visdate_bl) / 365.25]

# rename age_at_scan
setnames(berkadni, "age_at_scan", "age")

# specify the dates to use to match berkadni with the diagnosis data
berkadni[, joindate := scandate]
dxsum[, joindate := examdate]

# retrieve most recent cognitive diagnosis relative to scandate
# dx_scan = diagnosis as of scandate
berkadni <- dxsum[, .(rid, joindate, dx_scan = diagnosis, dx_date = examdate)
                  ][berkadni, on = .(rid, joindate), roll = TRUE]

# clean variables
berkadni[, `:=`(dx_bl_clean = fcase(dx_bl == 1, "CN",
                                    dx_bl == 2, "MCI",
                                    dx_bl == 3, "Dementia",
                                    default = NA_character_),
                dx_scan_clean = fcase(dx_scan == 1, "CN",
                                      dx_scan == 2, "MCI",
                                      dx_scan == 3, "Dementia",
                                      default = NA_character_),
                ptgender = fcase(ptgender == 1, "Male",
                                 ptgender == 2, "Female"))]

berkadni[, .N, keyby = .(dx_bl, dx_bl_clean)]
berkadni[, .N, keyby = .(dx_scan, dx_scan_clean)]
berkadni[, .N, keyby = .(dx_bl_clean, dx_scan_clean)]


##########################################################################################
## WRITE PRIVATE DATA OBJECTS ##
##########################################################################################

# create the directory if it does not exist
if (!dir.exists(PRIVATE_OUTPUT_DIR)) dir.create(PRIVATE_OUTPUT_DIR)
qs_save(berkadni, file.path(PRIVATE_OUTPUT_DIR, "berkadni.qs2"))
