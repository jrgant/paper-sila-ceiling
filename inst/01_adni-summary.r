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
berk <- unique(fread(file.path(ADNI_PATH, qp("UCBERKELEY_AMY_6MM"))))
names(berk) <- tolower(names(berk))

adnimerge <- fread(file.path(ADNI_PATH, qp("ADNIMERGE")))
names(adnimerge) <- tolower(names(adnimerge))

## subset to columns needed for estimation
berk <- berk[, .(rid, scandate, centiloids)]
adnimerge <- adnimerge[, .(rid, examdate, age, ptgender, dx, apoe4)]


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
berk_multi_scan |>
  ggplot(aes(x = days_since_last_scan / 365.25)) +
  geom_density(aes(color = scan_num, group = scan_num)) +
  scale_color_continuous(name = "Scan number") +
  theme_pander()

berk_multi_scan |>
  ggplot(aes(x = days_since_last_scan)) +
  geom_histogram(color = "white") +
  facet_wrap(vars(scan_num)) +
  theme_pander()

scan_lag_days <- berk_multi_scan[, density(days_since_last_scan, bw = "SJ")]


##########################################################################################
## CREATE EMPIRICAL ANALYSIS SAMPLE ##
##########################################################################################

berkadni <- merge(berk,
                  adnimerge[, .(
                    age_bl = first(age),
                    examdate_bl = first(examdate),
                    ptgender = first(ptgender),
                    apoe4 = first(apoe4),
                    dx_bl =  first(dx)
                  ), keyby = rid],
                  by = "rid",
                  all.x = TRUE)
NRID_MISS_AGE <- berkadni[is.na(age_bl), uniqueN(rid)]
berkadni <- berkadni[!is.na(age_bl)]
NRID_MISS_CENTILOIDS <- berkadni[is.na(centiloids), uniqueN(rid)]
NSCAN_MISS_CENTILOIDS <- berkadni[is.na(centiloids), .N]
berkadni <- berkadni[!is.na(centiloids)]

droptable <- data.table(
  reason  = c("miss_baseline_age", "miss_centiloids", "miss_centiloids"),
  rowtype = c("subid", "subid", "scan"),
  value   = c(NRID_MISS_AGE, NRID_MISS_CENTILOIDS, NSCAN_MISS_CENTILOIDS)
)


################################################################################
## SIMULATION PARAMETER: AGE AT FIRST SCAN ##
################################################################################

setkeyv(adnimerge, c("rid", "examdate"))

berkadni[, age_at_scan := round(age_bl + (scandate - examdate_bl) / 365.25, 1)]
setkeyv(berkadni, c("rid", "scandate"))

berkadni_first_scan <- berkadni[, .SD[1], keyby = rid]

## inspect symmetry of age distribution
berkadni_first_scan[, hist(age_at_scan)]

agedist_first_scan <- berkadni_first_scan[, .(
  mean = mean(age_at_scan),
  sd = sd(age_at_scan)
)]

agedist_first_scan[]


################################################################################
## COPY OBJECTS TO DATA DIRECTORY ##
################################################################################

usethis::use_data(
  num_scan_props,
  scan_lag_days,
  agedist_first_scan,
  droptable,
  overwrite = TRUE,
  internal = FALSE
)


##########################################################################################
## ADD VARIABLES to EMPIRICAL DATASET ##
##########################################################################################

berkadni[, yrs_since_bl := (scandate - examdate_bl) / 365.25]
berkadni[, age := round(age_bl + yrs_since_bl, digits = 1)]

berkadni[, joindate := scandate]
adnimerge[, joindate := examdate]

# retrieve most recent cognitive diagnosis relative to scandate
# dx_scan = diagnosis as of scandate
berkadni <- adnimerge[, .(rid, joindate, dx_scan = dx, dx_date = examdate)
                      ][berkadni, on = .(rid, joindate), roll = TRUE]

# clean dx variables
berkadni[, `:=`(dx_bl_clean = fcase(dx_bl == "Dementia", "AD",
                                    dx_bl == "CN", "CN",
                                    dx_bl == "MCI", "MCI",
                                    default = NA_character_),
                dx_scan_clean = fcase(dx_scan == "Dementia", "AD",
                                      dx_scan == "CN", "CN",
                                      dx_scan == "MCI", "MCI",
                                      default = NA_character_))]

berkadni[, .N, .(dx_bl, dx_bl_clean)]
berkadni[, .N, .(dx_scan, dx_scan_clean)]


##########################################################################################
## WRITE PRIVATE DATA OBJECTS ##
##########################################################################################

# create the directory if it does not exist
if (!dir.exists(PRIVATE_OUTPUT_DIR)) dir.create(PRIVATE_OUTPUT_DIR)
qs_save(berkadni, file.path(PRIVATE_OUTPUT_DIR, "berkadni.qs2"))
