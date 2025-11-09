################################################################################
## SETUP ##
################################################################################

library(data.table)
library(here)
library(ggplot2)
library(ggthemes)
library(usethis, include.only = "use_data")
library(collapse, include.only =  "descr")

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


################################################################################
## LOAD ADNI DATA ##
################################################################################

## Access to ADNI data requires a Memorandum of Understanding.
## See https://adni.loni.usc.edu for details.

## ADNI data are stored in a private directory outside of the package. The
## path to the directory is stored in the environment variable ADNI_PATH.

## NOTE: The call to unique() when reading in the Berkeley data drops 6
##       rows that appear to be duplicates.
ADNI_PATH <- Sys.getenv("ADNI_PATH")
berk <- unique(fread(file.path(ADNI_PATH, qp("UCBERKELEY_AMY_6MM"))))
adnimerge <- fread(file.path(ADNI_PATH, qp("ADNIMERGE")))

## names to lowercase
setnames(berk, names(berk), tolower(names(berk)))
setnames(adnimerge, names(adnimerge), tolower(names(adnimerge)))

## subset to columns needed for estimation
berk <- berk[, .(rid, scandate, centiloids)]
adnimerge <- adnimerge[, .(rid, examdate, ptgender, dx, apoe4)]


################################################################################
## PROPORTION with N SCANS ##
################################################################################

num_scan_props <-
  berk[, .(num_scans = .N), rid] |>
  _[, .N, keyby = num_scans] |>
  _[, P := N / sum(N)]

num_scan_props[]


################################################################################
## TIME BETWEEN SCANS ##
################################################################################

setkeyv(berk, c("rid", "scandate"))
berk[, scan_num := rowid(rid)]
berk[, days_since_last_scan := scandate - shift(scandate, n = 1, type = "lag"),
     keyby = rid]

## Quick summary
berk[scan_num > 1, .(mean = mean(days_since_last_scan)), keyby = scan_num]

## Subset to scans after the first for each patient
berk_multi_scan <- berk[scan_num > 1]

## Inspect time since last scan (in days) by scan number.
## Decision to use a single vector of scan lags is based on visual inspection,
## as these data are simply used to generate plausible amyloid curves.
berk_multi_scan |>
  ggplot(aes(x = days_since_last_scan)) +
  geom_density(aes(color = factor(scan_num))) +
  scale_color_few(name = "Scan number") +
  theme_pander()

berk_multi_scan |>
  ggplot(aes(x = days_since_last_scan)) +
  geom_histogram(color = "white") +
  facet_wrap(vars(scan_num)) +
  theme_pander()

scan_lag_days <- berk_multi_scan[, density(days_since_last_scan, bw = "SJ")]


################################################################################
## AGE AT FIRST SCAN ##
################################################################################

setkeyv(adnimerge, c("rid", "examdate"))

berkadni <- merge(berk,
                  adnimerge[, .(
                    age_bl = first(age),
                    examdate_bl = first(examdate),
                    ptgender = first(ptgender),
                    apoe4 = first(apoe4),
                    dx =  first(dx)
                  ), keyby = rid],
                  by = "rid",
                  all.x = TRUE)

berkadni[, age_at_scan := round(age_bl + (scandate - examdate_bl) / 365.25, 1)]
setkeyv(berkadni, c("rid", "scandate"))
berkadni[, rownum := rowid(rid)][]

berkadni_first_scan <- berkadni[rownum == 1 & !is.na(age_at_scan)]

berkadni_first_scan[, hist(age_at_scan)]

agedist_first_scan <- berkadni_first_scan[, .(mean = mean(age_at_scan),
                                              sd = sd(age_at_scan))]

agedist_first_scan[]


################################################################################
## COPY OBJECTS TO DATA DIRECTORY ##
################################################################################

usethis::use_data(
  num_scan_props,
  scan_lag_days,
  agedist_first_scan,
  overwrite = TRUE,
  internal = FALSE
)
