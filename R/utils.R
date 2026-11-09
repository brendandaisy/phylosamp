check_proportion <- function(p) {
  name_p <- rlang::enexpr(p)
  if (p < 0 | p > 1) {
    cli::cli_abort(c("`{name_p}` must be between 0 and 1.",
                     "x" = "You input `{name_p}` = {p}"),
                   call = NULL)
  }
}

check_m_N <- function(m, N) {
  if (m > N) {
    cli::cli_warn(
      c(
        "Necessary sample size not achievable.",
        "i" = "The necessary sample size given the input parameters is {round(m)}",
        "x" = "You input {.var N} = {N} as the total size of the outbreak"
      )
    )
    return(NA)
  }
}
check_positive <- function(x) {
  name_x <- rlang::enexpr(x)
  if (x <= 0) {
    cli::cli_abort(c("`{name_x}` must be positive.",
                     "x" = "You input `{name_x}` = {x}"),
                   call = NULL)
  }
}

# Variance of the maximum likelihood estimate of the logistic growth rate given
# one sample per time step from time 0 to t (i.e., the inverse Fisher information)
var_growth_rate <- function(t, p0_v1, r_v1) {
  ts <- 0:t
  q <- stats::plogis(r_v1 * ts - log((1 / p0_v1) - 1))
  return(1 / sum(ts^2 * q * (1 - q)))
}

# Derivative of variant prevalence at time t with respect to the logistic growth rate
dprev_dgrowth <- function(t, p0_v1, r_v1) {
  q <- stats::plogis(r_v1 * t - log((1 / p0_v1) - 1))
  return(t * q * (1 - q))
}

# Derivative of the dominance time (time to 50% prevalence) with respect to the logistic growth rate
ddominance_dgrowth <- function(p0_v1, r_v1) {
  a <- (1 / p0_v1) - 1
  return(-log(a) / r_v1^2)
}

# Absolute error in dominance time corresponding to a relative precision in prevalence at dominance,
# i.e., the time between the variant reaching 0.5*(1 - precision) and 0.5 prevalence
dominance_error <- function(precision, r_v1) {
  return(log((1 + precision) / (1 - precision)) / abs(r_v1))
}
