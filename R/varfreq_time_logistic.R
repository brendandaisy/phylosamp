##' Calculate time at which a variant reaches a given prevalence given logistic growth
##'
##' This function calculates the (possibly non-integer) number of time steps (e.g., days) until a variant
##' reaches a given observed prevalence, given a logistic growth rate and initial variant prevalence.
##' It is the inverse of [varfreq_freq_logistic()].
##'
##' @param p_v1 variant prevalence (proportion) for which to calculate the time it is reached
##' @param p0_v1 initial variant prevalence (# introductions / infected population size)
##' @param r_v1 logistic growth rate
##' @param c_ratio coefficient of detection ratio, calculated as the ratio of the coefficients of variant 1 to variant 2; default = 1 (no bias)
##' @return scalar giving the time step at which the variant reaches prevalence `p_v1`. Negative values indicate the prevalence was reached before time 0
##'
##' @author B. K. M. Case
##'
##' @examples
##' varfreq_time_logistic(p_v1 = 0.01, p0_v1 = 1/10000, r_v1 = 0.1, c_ratio = 1)
##'
##' @family logistic growth functions
##' @family variant frequency functions
##'
##' @export


varfreq_time_logistic <- function(p_v1, p0_v1, r_v1, c_ratio = 1) {
  if (!all(is.numeric(p_v1), p_v1 > 0 & p_v1 < 1)) stop("Variant prevalence must be numeric and between 0 and 1.")
  if (!all(is.numeric(p0_v1), p0_v1 > 0 & p0_v1 < 1)) stop("Initial variant prevalence must be numeric and between 0 and 1.")
  if (!all(is.numeric(r_v1), r_v1 != 0)) stop("Growth rate must be numeric and non-zero.")
  if (!all(is.numeric(c_ratio), c_ratio > 0)) stop("Coefficient of detection ratio must be numeric and greater than 0")

  a <- (1 / p0_v1) - 1
  b <- a * (1 / c_ratio)
  t <- -(1 / r_v1) * log(((1 / p_v1) - 1) / b)
  return(t)
}
