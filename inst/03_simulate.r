##########################################################################################
## SETUP ##
##########################################################################################

pkgload::load_all()
library(data.table)
library(qs2)

source(here::here("inst", "00_constants.r"))

# Print parameter settings and output directory
cat("\n", rep("-", 40), "\n", "Parameter settings:\n", rep("-", 40), "\n", sep = "")
cbind(value = sapply(grep("^(?!OUTPUT)", x = ls(), perl = TRUE, value = TRUE), get))
sapply(ls(pattern = "OUTPUT_DIR"), get)

# Function: A simple wrapper around simulate_curves() that sets the static parameters
#    used by all simulation scenarios.
simquick <- function(...) {
  simulate_curves(
    N = NDAT,
    centiloid_range = c(-20, Inf),
    scan_dist = num_scan_props,
    lag_dist = scan_lag_days,
    age_scan1 = agedist_first_scan,
    age_apos = list(mean = AGE_APOS_MU, sd = AGE_APOS_SD),
    apos_threshold = APOS_THRESHOLD,
    epsilon = e(rnorm(.N, sd = 5)),
    ...
  )
}

# Function: Set a unique seed each time the code is run
#   The seed will be recorded as an attribute on each simulated dataset
set_seed <- function() {
  seedint <- round(runif(1, -.Machine$integer.max, .Machine$integer.max))
  set.seed(seedint)
  seedint
}


##########################################################################################
## SIMULATE EXPONENTIAL AMYLOID CURVE SCENARIOS ##
##########################################################################################

## Scenario 1: Homogeneous inter-individual rates

s1_seed <- set_seed()
expRatesHomo <- runif(NSIM, min = EXP_RATE_MU_MIN, max = EXP_RATE_MU_MAX)

simexp_homo <- rbindlist(lapply(expRatesHomo, \(rate) {
  tmp <- simquick(
    genfun = gen_exponential,
    args = e(x = xvalue, k = k, x0 = x0, offset = offset),
    static_vars = substitute(
      expr = e(
        offset = EXP_OFFSET,
        k = curr_k,
        x0 = - (1 / k) * log(apos_threshold + 20)
      ),
      env = list(curr_k = rate)
    )
  )
}), idcol = "sim")
attr(simexp_homo, "rng_info") <- list(rng_kind = RNGkind(), seed = s1_seed)


## Scenario 2: Heterogeneous inter-individual rates

s2_seed <- set_seed()
expRatesHetero <- runif(NSIM, min = EXP_RATE_MU_MIN, max = EXP_RATE_MU_MAX)

simexp_hetero <- rbindlist(lapply(expRatesHetero, \(rate_mean) {
  tmp <- simquick(
    genfun = gen_exponential,
    args = e(x = xvalue, k = k, x0 = x0, offset = offset),
    static_vars = substitute(
      expr = e(
        offset = EXP_OFFSET,
        k = curr_k,
        x0 = - (1 / k) * log(apos_threshold + 20)
      ),
      env = list(curr_k = rnorm(NDAT, mean = rate_mean, sd = rate_mean * 0.1))
    )
  )
}), idcol = "sim")
attr(simexp_hetero, "rng_info") <- list(rng_kind = RNGkind(), seed = s2_seed)


##########################################################################################
## SIMULATE LOGISTIC AMYLOID CURVE SCENARIOS ##
##########################################################################################

## Scenario 3: Homogeneous inter-individual rates and maxima

s3_seed <- set_seed()
logFunMaxesHomo <- runif(NSIM, min = LOG_FMAX_MU_MIN, max = LOG_FMAX_MU_MAX)
logRatesHomo <- runif(NSIM, min = LOG_RATE_MU_MIN, max = LOG_RATE_MU_MAX)

simlog_homo <- rbindlist(lapply(seq_len(NSIM), \(i) {
  curr_fmax <- logFunMaxesHomo[i]
  curr_rate <- logRatesHomo[i]
  simquick(
  genfun = gen_logistic,
  args = e(x = xvalue, L = L, k = k, x0 = x0),
  static_vars = substitute(
    expr = e(L = fmax, k = frate, x0 = log((L / apos_threshold) - 1) / k),
    env = list(fmax = curr_fmax, frate = curr_rate)
  )
)}), idcol = "sim")
attr(simlog_homo, "rng_info") <- list(rng_kind = RNGkind(), seed = s3_seed)


## Scenario 4: Heterogeneous inter-individual rates and maxima

s4_seed <- set_seed()
logFunMaxesHetero <- runif(NSIM, min = LOG_FMAX_MU_MIN, max = LOG_FMAX_MU_MAX)
logRatesHetero <- runif(NSIM, min = LOG_RATE_MU_MIN, max = LOG_RATE_MU_MIN)

simlog_hetero <- rbindlist(lapply(seq_len(NSIM), \(i) {
  curr_fmax_mean <- logFunMaxesHetero[i]
  curr_rate_mean <- logRatesHetero[i]
  simquick(
    genfun = gen_logistic,
    args = e(x = xvalue, L = L, k = k, x0 = x0),
    static_vars = substitute(
      expr = e(
        L = rnorm(NDAT, mean = fmax_mean, sd = fmax_mean * 0.1),
        k = rnorm(NDAT, mean = rate_mean, sd = rate_mean * 0.1),
        x0 = log((L / apos_threshold) - 1) / k
      ),
      env = list(fmax_mean = curr_fmax_mean, rate_mean = curr_rate_mean)
    )
  )
}), idcol = "sim")
attr(simlog_hetero, "rng_info") <- list(rng_kind = RNGkind(), seed = s4_seed)


##########################################################################################
## WRITE SIMULATED DATASETS ##
##########################################################################################

simlist <- lapply(setNames(nm = ls(pattern = "simexp|simlog")), get)
qs_save(simlist, file = file.path(OUTPUT_DIR, "simulated-datasets.qs2"))
