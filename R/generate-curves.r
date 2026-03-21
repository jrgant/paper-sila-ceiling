#' @title Curve-generating functions
#'
#' @param x Point at which to evaluate function
#' @param L Function maximum
#' @param k Rate parameter
#' @param x0 Midpoint of the exponentential and logistic curve
#' @param b Y-intercept for the linear model
#' @param offset For the exponential function, the vertical shift. For the logistic
#'   function, must be <= 0 (the default); setting it allows the logistic function
#'   to elicit values < 0 while maintaining the desired maximum.
#'
#' @rdname generate-curves
#' @export
gen_logistic <- function(x, L, k = 1, x0 = 0, offset = 0) {
  stopifnot(offset <= 0)
  ((L + abs(offset)) / (1 + exp(-k * (x - x0)))) + offset
}

#' @rdname generate-curves
#' @export
gen_exponential <- function(x, k, x0, offset = 0) {
  offset + exp(k * (x - x0))
}

#' @rdname generate-curves
#' @export
gen_linear <- function(x, k, b) {
  k * x + b
}
