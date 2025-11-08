#' @title Curve-generating functions
#'
#' @param x Point at which to evaluate function
#' @param L Function maximum
#' @param k Rate parameter
#' @param x0 Midpoint of the logistic curve
#'
#' @rdname generate-curves
#' @export
gen_logistic <- function(x, L, k = 1, x0 = 0) {
  (L / (1 + exp(-k * (x - x0))))
}

#' @param offset Vertical shift
#' @rdname generate-curves
#' @export
gen_exponential <- function(x, k, x0, offset) {
  offset + exp(k * (x - x0))
}
