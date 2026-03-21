##########################################################################################
## SETUP and DATA ##
##########################################################################################

library(data.table)
library(silaR)
library(qs2)

source(here::here("inst", "00_constants.r"))
berkadni <- qs_read(file.path(PRIVATE_OUTPUT_DIR, "berkadni.qs2"))
setkeyv(berkadni, c("rid", "scandate"))


##########################################################################################
## RUN SILA ##
##########################################################################################

# IDs for subjects with multiple PET scans, non-missing age, and at least two non-missing
# centiloid measurements
lu_subid_all_multiscan <- berkadni[, .(
  num_scans = .N,
  two_cl = sum(!is.na(centiloids)) > 1
), by = rid][num_scans > 1 & two_cl == TRUE][, rid]

berksub <- berkadni[rid %in% lu_subid_all_multiscan]
setnames(berksub, c("rid", "centiloids"), c("subid", "val"))

# NOTE: Fitting SILA will throw warnings due to perfect fits in some of the drawn
#       samples. Traceable to lm() and LOESS fits.
set.seed(987312344)
fit <- list()
fit[["res"]] <- sila(berksub, dt = SILA_DT, val0 = SILA_VAL0, maxi = SILA_MAX_ITER)
fit[["resfit"]] <- sila_estimate(fit$res$tsila,
                                 df = berkadni[!is.na(centiloids),
                                               .(subid = rid, val = centiloids, age)],
                                 align_event = "last")


##########################################################################################
## WRITE ##
##########################################################################################

qs_save(fit, file.path(PRIVATE_OUTPUT_DIR, "sila_empirical_sample_berkeley.qs2"))
