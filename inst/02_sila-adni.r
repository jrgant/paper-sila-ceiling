##########################################################################################
## SETUP and DATA ##
##########################################################################################

library(data.table)
library(silaR)
library(qs2)

source(here::here("inst", "00_constants.r"))

# ADNIMERGE
adnimerge <- fread(file.path(ADNI_PATH, qp("ADNIMERGE")))
adnimerge[, joindate := EXAMDATE]
names(adnimerge) <- tolower(names(adnimerge))
setkeyv(adnimerge, c("rid", "examdate"))

## Extract baseline dataset
adni_age_bl <- adnimerge[, .(
  bl_exam_date = first(examdate),
  age_bl = first(age),
  dx_bl = first(dx_bl)
), keyby = rid]

# AMYLOID PET SCANs (UC BERKELEY)
berk <- fread(file.path(ADNI_PATH, qp("UCBERKELEY_AMY_6MM")))
names(berk) <- tolower(names(berk))

# ANALYSIS DATASET

## Inner join, dropping missing centiloids and pulling in baseline information
AOUT_NMISS_CENTILOIDS <- berk[, sum(is.na(centiloids))]
berkadni <- berk[!is.na(centiloids), .(rid, scandate, centiloids)
                 ][adni_age_bl, on = .(rid), nomatch = NULL]

## Include most recent DX status at first scan date
berkadni_first_scan <- berkadni[, .SD[1], keyby = .(rid, scandate)][, .(rid, scandate)]
berkadni_first_scan[, joindate := scandate]
berkadni_first_scan <- adnimerge[berkadni_first_scan,
                                 on = .(rid, joindate), roll = TRUE
                                 ][, .(rid, scandate, examdate, dx)]

## Create and add relevant variables to analytic dataset
berkadni[, yrs_since_bl := (scandate - bl_exam_date) / 365.25]
berkadni[, age := round(age_bl + yrs_since_bl, digits = 1)]
berkadni[berkadni_first_scan, on = .(rid), j = `:=`(dx_fs = i.dx,
                                                    dx_fs_date = i.examdate)]

berkadni[, `:=`(dx_bl_grp = fcase(dx_bl == "AD", "AD",
                                  dx_bl == "CN", "CN",
                                  dx_bl %in% c("EMCI", "LMCI", "SMC"), "MCI/SMC",
                                  default = NA_character_),
                dx_fs_grp = fcase(dx_fs == "AD", "AD",
                                  dx_fs == "CN", "CN",
                                  dx_fs %in% c("EMCI", "LMCI", "SMC"), "MCI/SMC",
                                  default = NA_character_))]

berkadni[is.na(dx_fs_grp), dx_fs_grp := dx_bl_grp]

berkadni[, .SD[1], keyby = rid][, .N, .(dx_bl_grp, dx_fs_grp)]


##########################################################################################
## RUN SILA ##
##########################################################################################

# IDs for subjects with multiple PET scans
lu_subid_all_multiscan <- berkadni[, .N, keyby = rid][N > 1, rid]
berksub <- berkadni[rid %in% lu_subid_all_multiscan]
setnames(berksub, c("rid", "centiloids"), c("subid", "val"))

# NOTE: Fitting SILA will throw warnings due to perfect fits in some of the drawn
#       samples. Traceable to lm() and LOESS fits.
set.seed(987312344)
fit <- list()
fit[["res"]] <- sila(berksub, dt = SILA_DT, val0 = SILA_VAL0, maxi = SILA_MAX_ITER)
fit[["resfit"]] <- sila_estimate(fit$res$tsila, df = berksub)


##########################################################################################
## WRITE ##
##########################################################################################

if (!dir.exists(PRIVATE_OUTPUT_DIR)) dir.create(PRIVATE_OUTPUT_DIR)

qs_save(fit, file.path(PRIVATE_OUTPUT_DIR, "sila_empirical_sample_berkeley.qs2"))
fwrite(berksub, file.path(PRIVATE_OUTPUT_DIR, "berkeley_scans_formatted.csv"))
