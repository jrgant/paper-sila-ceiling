# Control parameters
NSIM <- 100
NDAT <- 2000

# Amyloid positivity centiloid threshold
APOS_THRESHOLD <- 20

# Exponential parameters
EXP_RATE_MU_MIN <- 0.01
EXP_RATE_MU_MAX <- 0.07
EXP_OFFSET      <- (-20)

# Logistic parameters
LOG_FMAX_MU_MIN <- 130
LOG_FMAX_MU_MAX <- 165
LOG_RATE_MU_MIN <- 0.2
LOG_RATE_MU_MAX <- 0.8

# Common parameters
SIGMA_HETERO_MULTIPLIER <- 0.1
NOISE_MU <- 0
NOISE_SD <- 5

## Source for mean age at amyloid positivity (and standard deviation):
##   Betthauser, Bilgel & Koscik et al. (2022-11) Multi-method investigation of factors
##   influencing amyloid onset and impairment in three cohorts, Brain.
AGE_APOS_MU <- 68.5
AGE_APOS_SD <- 8.7

# SILA SETTINGS
SILA_DT        <- 0.25
SILA_VAL0      <- APOS_THRESHOLD
SILA_MAX_ITER  <- 200

# PUBLIC RESULTS OUTPUT
OUTPUT_DIR <- here::here("output")

# PRIVATE ADNI DIRECTORIES
ADNI_DL_DATE <- "27Oct2025"
ADNI_PATH <- Sys.getenv("ADNI_PATH")
PRIVATE_OUTPUT_DIR <- file.path(ADNI_PATH, "sila-ceiling-analysis")

# Function: Helper to paste ADNI source dataset filename into a path
qp <- function(string) {
  paste0(string, "_", ADNI_DL_DATE, ".csv")
}
