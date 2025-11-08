#' Simulate amyloid curves
#'
#' @param N Number of amyloid curves to simulate
#' @param centiloid_range A numeric vector of length 2 specifiying the minimum
#'   and maximum possible centiloid values.
#' @param scan_dist A data.frame specifying the probability of a patient having
#'   x PET scans.
#' @param lag_dist An estimated probability density (output by `density()`)
#'   for the number of days between two of a given individual's PET scans.
#' @param age_scan1 A list specifying the mean and standard deviation
#'   for the age distribution at first PET scan. E.g. `list(mean = 50, sd = 5)`.
#' @param age_apos A list specifying the mean and standard deviation for the age
#'   at which a patient's amyloid level becomes positive, based on the
#'   centiloid threshold specified in `Apos_threshold`.
#' @param apos_treshhold The centiloid value at which a scan/individual is
#'   considered amyloid-positive.
#' @param genfun A function that generates the amyloid curve.
#' @param static_vars A named expression specifying input parameters to
#'   `genfun()`, that remain fixed for each individual. These inputs may be
#'   scalars or expressions specifying operations to execute in the initial
#'   dataset created within the this function.
#' @param args A named expression specifying the arguments to be passed to
#'   `genfun`.
#' @param epsilon An expression that, when evaluated, will generate noise
#'   at a given point along the amyloid curve. The resulting measurement
#'   typically stands in as the observed amyloid level.
#'
#' @import data.table
#' @export simulate_curves
simulate_curves <- function(N,
                            centiloid_range = NULL,
                            scan_dist = num_scan_props,
                            lag_dist = scan_lag_days,
                            age_scan1 = agedist_first_scan,
                            age_apos,
                            apos_threshold,
                            genfun = NULL,
                            args = NULL,
                            static_vars = NULL,
                            epsilon = expression(rnorm(.N, sd = 4))) {

  if (length(centiloid_range) != 2) {
    stop("`centiloid_range` should be a vector of length 2 specifying ",
         "the minimum and maximum possible values for the simulated ",
         "centiloid measurements.")
  }

  if (!is.function(genfun)) {
    stop("The argument to `genfun` must be a function.")
  }

  if (!is.expression(args)) {
    stop("The argument to `args` must be an expression.")
  }

  if (!is.expression(static_vars) && !is.call(static_vars)) {
    stop("The argument to `static_vars` must be an expression or a call.")
  }

  if (!is.expression(epsilon)) {
    stop("The argument to `epsilon` must be an expression.")
  }

  # initial dataset with parameters external to the curve-generating function set
  init <- data.table(subid = seq_len(N),
                     num_scans = sample(scan_dist$num_scans,
                                        size = N,
                                        replace = TRUE,
                                        prob = scan_dist$P),
                     age_apos = rnorm(N, mean = age_apos$mean, sd = age_apos$sd),
                     age_scan1 = rnorm(N, age_scan1$mean, age_scan1$sd))

  # handle the case of static_vars being provided as a call
  varnames <- names(static_vars)[names(static_vars) != ""]
  lapply(varnames, \(var) {
    init[, (var) := eval(static_vars[[var]])]
  })

  # transform the initial dataset to long form
  long <- init[rep(subid, num_scans)]
  long[, scan_num := rowid(subid)]

  # draw scan lags from input density (for scans after the first)
  long[scan_num == 1, lag_days := 0L]
  long[scan_num > 1, lag_days := as.integer(
    sample(lag_dist$x,
           size = .N,
           replace = TRUE,
           prob = lag_dist$y) + rnorm(.N, 0, lag_dist$bw)
  )]
  long[, cum_lag_days := cumsum(lag_days), by = subid]

  # calculate scan time points
  # (relative to disease stage, i.e., time at amyloid positivity onset)
  long[, t0 := age_scan1 - age_apos]
  long[, xvalue := t0 + (cum_lag_days / 365.25)]

  # calculate age at scan
  long[, age := age_scan1 + (cum_lag_days / 365.25)]

  # generate observations at xvalue
  long[, centiloids := eval(bquote(genfun(..(args)), splice = TRUE))]
  long[, centiloids_measured := centiloids + eval(epsilon)]

  setkeyv(long, c("subid", "scan_num"))
  long[]
}
