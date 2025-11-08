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
  lu_multiscan <- simdata[, .N, .(sim, subid)][N > 1]
  multiscan_subset <- simdata[lu_multiscan, on = .(sim, subid)]
  stopifnot(all(multiscan_subset[, N] > 1))
  out <- multiscan_subset[, .(sim, subid, xvalue, val = centiloids_measured)]
  setkeyv(out, c("sim", "subid", "xvalue"))
  out
}

#  Function
fit_sila <- function(fmt_data) {
  ## print(names(fmt_data))
  ## print(class(fmt_data))
  plan(multisession) # start parallel processing
  fits <- foreach(i = unique(fmt_data$sim),
                  .options.future = list(packages = c("data.table"))) %dofuture% {
    tmp <- fmt_data[sim == i]
    currfit <- list()
    currfit[["res"]] <- sila(tmp, dt = SILA_DT, val0 = SILA_VAL0, maxi = SILA_MAX_ITER)
    currfit[["resfit"]] <- sila_esimate(tmp$res$tsila, df = tmp)
  }
  plan(sequential) # end parallel processing
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

fit_sila(fmt_simexp_homo)
