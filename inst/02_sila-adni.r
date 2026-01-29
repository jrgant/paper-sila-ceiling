##########################################################################################
## SETUP and DATA ##
##########################################################################################

library(data.table)
library(silaR)
library(qs2)

source(here::here("inst", "00_constants.r"))
berkadni <- qs_read(file.path(PRIVATE_OUTPUT_DIR, "berkadni.qs2"))


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
fit[["resfit"]] <- sila_estimate(fit$res$tsila, df = berksub, align_event = "all")

##########################################################################################
## WRITE ##
##########################################################################################

qs_save(fit, file.path(PRIVATE_OUTPUT_DIR, "sila_empirical_sample_berkeley.qs2"))
