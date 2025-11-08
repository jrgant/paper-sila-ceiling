#' Alias and helpers for common operations

#' @describeIn helpers Alias for [base::expression()]
#' @export
e <- base::expression


#' @param object Any object with an attribute named "rng_info".
#' @describeIn helpers Retrieve RNG info for a simulated object
#' @export
rng_info <- function(object) {
  attr(object, "rng_info")
}
