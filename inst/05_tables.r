##########################################################################################
## SETUP ##
##########################################################################################

library(data.table)
library(qs2)
library(flextable)
library(collapse, include.only = "descr")

source(here::here("inst", "00_constants.r"))

berkadni <- qs_read(file.path(PRIVATE_OUTPUT_DIR, "berkadni.qs2"))
empsila <- qs_read(file.path(PRIVATE_OUTPUT_DIR, "sila_empirical_sample_berkeley.qs2"))

# IDs for subjects with multiple PET scans, non-missing age, and at least two non-missing
# centiloid measurements
lu_subid_all_multiscan <- berkadni[, .(
  num_scans = .N,
  miss_age = is.na(age_at_scan[1]),
  two_cl = sum(!is.na(centiloids)) > 1
), by = rid][num_scans > 1 & miss_age == FALSE & two_cl == TRUE][, rid]

# identify individuals with more than 1 PET scan
lu_subid_multiscan <- berkadni[, .N, rid][N > 1, rid]


##########################################################################################
## SUMMARIES ##
##########################################################################################

catvars <- c("dx_bl_clean", "dx_scan_clean")
contvars <- c("scandate", "yrs_since_bl", "age")

sum_fs_all <- collapse::descr(
  firstscan[, .SD, .SDcols = c(catvars, contvars)]
)

sum_fs_one <- collapse::descr(
  firstscan[rid %notin% lu_subid_multiscan, .SD, .SDcols = c(catvars, contvars)]
)

sum_fs_multi <- collapse::descr(
  firstscan[rid %in% lu_subid_multiscan, .SD, .SDcols = c(catvars, contvars)]
)

