##########################################################################################
## SETUP and DATA ##
##########################################################################################
source(here::here("inst", "00_constants.r"))

library(data.table)
library(silaR)
library(qs2)

ADNI_DL_DATE <- "27Oct2025"

# Function: Helper to paste ADNI source dataset filename into a path
qp <- function(string) {
  paste0(string, "_", ADNI_DL_DATE, ".csv")
}

# ADNIMERGE
adnimerge <- fread(file.path(ADNI_PATH, qp("ADNIMERGE")))
adnimerge[, joindate := EXAMDATE]
setkeyv(adnimerge, c("RID", "EXAMDATE"))

## Extract baseline dataset
adni_age_bl <- adnimerge[, .(
  BL_EXAM_DATE = first(EXAMDATE),
  AGE_bl = first(AGE),
  DX_bl = first(DX_bl)
), keyby = RID]

# AMYLOID PET SCANs (UC BERKELEY)
berk <- fread(file.path(ADNI_PATH, qp("UCBERKELEY_AMY_6MM")))

# ANALYSIS DATASET

## Inner join, dropping missing centiloids and pulling in baseline information
AOUT_NMISS_CENTILOIDS <- berk[, sum(is.na(CENTILOIDS))]
berkadni <- berk[!is.na(CENTILOIDS), .(RID, SCANDATE, CENTILOIDS)
                 ][adni_age_bl, on = .(RID), nomatch = NULL]

## Include most recent DX status at first scan date
berkadni_first_scan <- berkadni[, .SD[1], keyby = .(RID, SCANDATE)][, .(RID, SCANDATE)]
berkadni_first_scan[, joindate := SCANDATE]
berkadni_first_scan <- adnimerge[berkadni_first_scan,
                                 on = .(RID, joindate), roll = TRUE
                                 ][, .(RID, SCANDATE, EXAMDATE, DX)]

## Create and add relevant variables to analytic dataset
berkadni[, yrs_since_bl := (SCANDATE - BL_EXAM_DATE) / 365.25]
berkadni[, age := round(AGE_bl + yrs_since_bl, digits = 1)]
berkadni[berkadni_first_scan, on = .(RID), j = `:=`(DX_fs = i.DX,
                                                    DX_fs_date = i.EXAMDATE)]

berkadni[, `:=`(DX_bl_grp = fcase(DX_bl == "AD", "AD",
                                  DX_bl == "CN", "CN",
                                  DX_bl %in% c("EMCI", "LMCI", "SMC"), "MCI/SMC",
                                  default = NA_character_),
                DX_fs_grp = fcase(DX_fs == "AD", "AD",
                                  DX_fs == "CN", "CN",
                                  DX_fs %in% c("EMCI", "LMCI", "SMC"), "MCI/SMC",
                                  default = NA_character_))]

berkadni[is.na(DX_fs_grp), DX_fs_grp := DX_bl_grp]

berkadni[, .SD[1], keyby = RID][, .N, .(DX_bl_grp, DX_fs_grp)]


##########################################################################################
## RUN SILA ##
##########################################################################################

# IDs for subjects with multiple PET scans
lu_subid_all_multiscan <- berkadni[, .N, keyby = RID][N > 1, RID]
berksub <- berkadni[RID %in% lu_subid_all_multiscan]
setnames(berksub, c("RID", "CENTILOIDS"), c("subid", "val"))

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
