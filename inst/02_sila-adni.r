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

berksub <- berkadni[rid %in% lu_subid_all_multiscan
                    ][, .(subid = rid, age, val = centiloids, years_death)]

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
## RUN SILA (IMPUTING SCANS AMONG DECEDENTS) ##
##########################################################################################

last_scan_died <- berksub[!is.na(years_death) & years_death < -1, .(
  age = last(age),
  val = last(val),
  years_death = last(years_death)
), keyby = subid]

last_scan_died[, num_scans_add := floor(abs(years_death) / 2)]

last_scan_died <- last_scan_died[, .SD[rep(seq_len(.N), num_scans_add)],
                                 keyby = subid
                                 ][, new_scan := rowid(subid)][]
setnames(last_scan_died, "age", "age_last_scan")

last_scan_died[, age := age_last_scan + new_scan * 2]

# NOTE: The final imputation carries forward the final centiloid measurement
#   for each decedent to impose the extreme assumption that anyone who died and for whom
#   we don't have a scan every two years until their death, would have flattened out on
#   the amyloid curve.
berksub_impute <- rbind(berksub[, -c("years_death")],
                        last_scan_died[, .(subid, age, val)])

# NOTE: Fitting SILA will throw warnings due to perfect fits in some of the drawn
#       samples. Traceable to lm() and LOESS fits.
set.seed(245292734)
fit_death <- list()
fit_death[["res"]] <- sila(berksub_impute,
                           dt = SILA_DT, val0 = SILA_VAL0, maxi = SILA_MAX_ITER)
fit_death[["resfit"]] <- sila_estimate(
  fit_death$res$tsila,
  df = rbind(berkadni[!is.na(centiloids), .(subid = rid, val = centiloids, age)],
             last_scan_died[, .(subid, age, val)]),
  align_event = "last"
)


##########################################################################################
## WRITE ##
##########################################################################################

qs_save(fit,
        file.path(PRIVATE_OUTPUT_DIR, "sila_empirical_sample_berkeley.qs2"))

qs_save(fit_death,
        file.path(PRIVATE_OUTPUT_DIR, "sila_empirical_sample_impute_berkeley.qs2"))
