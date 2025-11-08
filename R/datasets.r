#' Distribution of the number of PET scans
#'
#' @description
#' A data.table containing the distribution of the number of PET scans
#' available for patients in the U.C. Berkeley amyloid PET scan dataset.
#'
#' @format A data.table with 7 rows and 3 columns.
#'
#' \itemize{
#'   \item num_scans. `<integer>` Number of PET scans.
#'   \item N. `<integer>` Number of patients in the U.C. Berkeley scan data with
#'      n PET scans available for analysis.
#'   \item P. `<double>` The proportion of patients in the U.C. Berkeley scan data
#'      with n PET scans available for analysis.
#' }
#'
#' @details
#'
#' Source: UC Berkeley - Amyloid PET 6mm Res analysis \[ADNI1,GO,2,3,4\]
#'
#' \href{https://ida.loni.usc.edu/explore/jsp/search/search.jsp?project=ADNI#studyFiles}{https://ida.loni.usc.edu/explore/jsp/search/search.jsp?project=ADNI#studyFiles}
"num_scan_props"


#' Age distribution at first PET scan
#'
#' @description
#' A data.table containing two columns specifying the mean age at one's first
#' PET scan, based on the ADNIMERGE and U.C. Berkeley scan datasets.
#'
#' @format A data.table with 1 row and 2 columns.
#'
#' \itemize{
#'   \item mean. `<double>` Mean age at first PET scan.
#'   \item sd. `<double>` Standard deviation of the age at first PET scan.
#' }
#'
#' @details
#'
#' Sources:
#'
#' UC Berkeley - Amyloid PET 6mm Res analysis \[ADNI1,GO,2,3,4\]
#' ADNIMERGE - Key ADNI tables merged into one table \[ADNI1,GO,2,3\]
#'
#' \href{https://ida.loni.usc.edu/explore/jsp/search/search.jsp?project=ADNI#studyFiles}{https://ida.loni.usc.edu/explore/jsp/search/search.jsp?project=ADNI#studyFiles}
"agedist_first_scan"


#' Distribution of time between scans
#'
#' @description
#' The density of days since an individual's prior PET scan (for scans after the
#' first), based on the U.C. Berkeley amyloid PET scan dataset.
#'
#' @format An object of class `density` produced by `stats::density()`.
#'
#' @details
#' Source: UC Berkeley - Amyloid PET 6mm Res analysis \[ADNI1,GO,2,3,4\]
#'
#' \href{https://ida.loni.usc.edu/explore/jsp/search/search.jsp?project=ADNI#studyFiles}{https://ida.loni.usc.edu/explore/jsp/search/search.jsp?project=ADNI#studyFiles}
"scan_lag_days"
