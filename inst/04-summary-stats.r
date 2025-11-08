##########################################################################################
## SETUP and DATA ##
##########################################################################################
source(here::here("inst", "00_constants.r"))

library(data.table)
library(qs2)

berkadni <- fread(file.path(PRIVATE_OUTPUT_DIR, "berkeley_scans_formatted.csv"))


##########################################################################################
## BASELINE DESCRIPTIVES ##
##########################################################################################

berkadni[, num_scans := .N, keyby = subid]
setkeyv(berkadni, c("subid", "SCANDATE"))

# NOTE: A handful of scans (n=37) occurred shortly before the first exam date
desc_first_scan  <- berkadni[, .SD[1],
                             .SDcols = patterns("age|num_scans|DX_bl_grp|yrs|SCANDATE"),
                             keyby = subid]

desc_fs_sum <- summary(desc_first_scan[, .(first_scan_date = SCANDATE,
                                          age_at_first_scan = age,
                                          num_scans)])

desc_fs_sum |> dim()
