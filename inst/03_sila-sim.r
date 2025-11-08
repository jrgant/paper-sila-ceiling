##########################################################################################
## SETUP ##
##########################################################################################

pkgload::load_all()
library(data.table)
library(doFuture)
library(qs2)
library(silaR)

source(here::here("inst", "00_constants.r"))

simlist <- qs_read(file.path(OUTPUT_DIR, "simulated-datasets.qs2"))

# Function: Subset the data to include individuals with multiple scans within a given
#    simulated dataset. Input must be an element of `simulated-datasets.qs2`.
fmt_sim <- function(simdata) {
  # impose a centiloid cap similar to maximum centiloids found in ADNI
  # SILA does not handle the occasional extremes produced by exponential simulations
  tmp <- simdata[centiloids_measured <= 350]
  lu_multiscan <- tmp[, .N, .(sim, subid)][N > 1]
  multiscan_subset <- simdata[lu_multiscan, on = .(sim, subid)]
  stopifnot(all(multiscan_subset[, N] > 1))
  out <- multiscan_subset[, .(sim, subid, age, val = centiloids_measured)]
  setkeyv(out, c("sim", "subid", "age"))
  out
}

#  Function: Run SILA over simulated datasets
fit_sila <- function(fmt_data, parallel = TRUE) {
  if (parallel == FALSE) {
    plan(sequential)
  } else {
    plan(multisession) # start parallel processing
  }
  fits <- foreach(i = unique(fmt_data$sim),
                  .options.future = list(packages = c("data.table", "silaR"),
                                         seed = TRUE)) %dofuture% {
    tmp <- fmt_data[sim == i]
    currfit <- list()
    currfit[["res"]] <- sila(tmp, dt = SILA_DT, val0 = SILA_VAL0, maxi = SILA_MAX_ITER)
    currfit[["resfit"]] <- sila_estimate(currfit$res$tsila, df = tmp)
    currfit
  }
  if (parallel == TRUE) {
    plan(sequential) # end parallel processing
  }
  fits
}


##########################################################################################
## FORMAT DATASETS ##
##########################################################################################

fmt_simexp_homo   <- fmt_sim(simlist$simexp_homo)
fmt_simexp_hetero <- fmt_sim(simlist$simexp_hetero)
fmt_simlog_homo   <- fmt_sim(simlist$simlog_homo)
fmt_simlog_hetero <- fmt_sim(simlist$simlog_hetero)


##########################################################################################
## RUN SILA ##
##########################################################################################

fitlist <- list()

fitlist[["simexp_homo"]]   <- fit_sila(fmt_simexp_homo)
fitlist[["simexp_hetero"]] <- fit_sila(fmt_simexp_homo)
fitlist[["simlog_homo"]]   <- fit_sila(fmt_simexp_homo)
fitlist[["simlog_hetero"]] <- fit_sila(fmt_simexp_homo)


##########################################################################################
## WRITE SILA ##
##########################################################################################

qs_save(fitlist, file.path(OUTPUT_DIR, "simulated-sila-fits.qs2"))
