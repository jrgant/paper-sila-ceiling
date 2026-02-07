##########################################################################################
## SETUP ##
##########################################################################################

pkgload::load_all()
library(data.table)
library(qs2)
library(MASS, include.only = "gamma.shape")

source(here::here("inst", "00_constants.r"))
empfit <- qs_read(file.path(PRIVATE_OUTPUT_DIR, "sila_empirical_sample_berkeley.qs2"))

# Fit an intercept-only model for the distribution of age at amyloid positivity, as
# estimated by the SILA fit in the empirical dataset.
# Shape and rate parameter calculation based on gfit$family$simulate
estaget0 <- as.data.table(
  empfit$resfit
)[, .(estaget0 = first(estaget0)), keyby = subid][, estaget0]

gfit <- glm(estaget0 ~ 1, family = Gamma(link = "identity"))
AGE_APOS_SHAPE <- gamma.shape(gfit)$alpha # 25.31755
AGE_APOS_RATE  <- unname(gamma.shape(gfit)$alpha / coef(gfit)[1]) # 0.3239963

# Print parameter settings and output directory
cat("\n", rep("-", 40), "\n", "Parameter settings:\n", rep("-", 40), "\n", sep = "")
cbind(value = sapply(grep("^(?!OUTPUT|qp)", x = ls(), perl = TRUE, value = TRUE), get))
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
    age_apos_fun = "rgamma",
    age_apos_args = list(shape = AGE_APOS_SHAPE, rate = AGE_APOS_RATE),
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

exp_x0_calc <- function(k, apos) {
  - (1 / k) * log(apos + 20)
}

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
        x0 = exp_x0_calc(k = k, apos = apos_threshold)
      ),
      env = list(curr_k = rate)
    )
  )
}), idcol = "sim")
attr(simexp_homo, "rng_info") <- list(rng_kind = RNGkind(), seed = s1_seed)
attr(simexp_homo, "params") <- data.table(
  sim    = seq_len(NSIM),
  k      = expRatesHomo,
  x0     = exp_x0_calc(expRatesHomo, apos = APOS_THRESHOLD),
  offset = EXP_OFFSET
)

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
      env = list(
        curr_k = rnorm(NDAT, mean = rate_mean, sd = rate_mean * SIGMA_HETERO_MULTIPLIER)
      )
    )
  )
}), idcol = "sim")
attr(simexp_hetero, "rng_info") <- list(rng_kind = RNGkind(), seed = s2_seed)
attr(simexp_hetero, "params") <- data.table(
  sim    = seq_len(NSIM),
  k      = expRatesHetero,
  x0     = exp_x0_calc(expRatesHetero, apo = APOS_THRESHOLD),
  offset = EXP_OFFSET
)


##########################################################################################
## SIMULATE LOGISTIC AMYLOID CURVE SCENARIOS ##
##########################################################################################

log_x0_calc <- function(L, k, apos) {
  log((L / apos) - 1) / k
}

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
    expr = e(L = fmax, k = frate, x0 = log_x0_calc(L, k, apos_threshold)),
    env = list(fmax = curr_fmax, frate = curr_rate)
  )
)}), idcol = "sim")
attr(simlog_homo, "rng_info") <- list(rng_kind = RNGkind(), seed = s3_seed)
attr(simlog_homo, "params") <- data.table(
  sim = seq_len(NSIM),
  L   = logFunMaxesHomo,
  k   = logRatesHomo,
  x0  = log_x0_calc(L = logFunMaxesHomo, k = logRatesHomo, apos = APOS_THRESHOLD)
)

## Scenario 4: Heterogeneous inter-individual rates and maxima

s4_seed <- set_seed()
logFunMaxesHetero <- runif(NSIM, min = LOG_FMAX_MU_MIN, max = LOG_FMAX_MU_MAX)
logRatesHetero <- runif(NSIM, min = LOG_RATE_MU_MIN, max = LOG_RATE_MU_MAX)

simlog_hetero <- rbindlist(lapply(seq_len(NSIM), \(i) {
  curr_fmax_mean <- logFunMaxesHetero[i]
  curr_rate_mean <- logRatesHetero[i]
  simquick(
    genfun = gen_logistic,
    args = e(x = xvalue, L = L, k = k, x0 = x0),
    static_vars = substitute(
      expr = e(
        L = rnorm(NDAT, mean = fmax_mean, sd = fmax_mean * SIGMA_HETERO_MULTIPLIER),
        k = rnorm(NDAT, mean = rate_mean, sd = rate_mean * SIGMA_HETERO_MULTIPLIER),
        x0 = log((L / apos_threshold) - 1) / k
      ),
      env = list(fmax_mean = curr_fmax_mean, rate_mean = curr_rate_mean)
    )
  )
}), idcol = "sim")
attr(simlog_hetero, "rng_info") <- list(rng_kind = RNGkind(), seed = s4_seed)
attr(simlog_hetero, "params") <- data.table(
  sim = seq_len(NSIM),
  L   = logFunMaxesHetero,
  k   = logRatesHetero,
  x0  = log_x0_calc(L = logFunMaxesHetero, k = logRatesHetero, apos = APOS_THRESHOLD)
)


##########################################################################################
## WRITE SIMULATED DATASETS ##
##########################################################################################

simlist <- lapply(setNames(nm = ls(pattern = "simexp|simlog")), get)
qs_save(simlist, file = file.path(OUTPUT_DIR, "simulated-datasets.qs2"))
