##' Calculate time needed to estimate variant prevalence given a per-timestep sample size
##'
##' This function calculates the minimum number of time steps (e.g., days) of periodic sampling
##' needed to estimate the current prevalence of a variant with a desired confidence and precision,
##' given a per-timestep sample size.
##' It assumes a periodic sampling strategy, where samples are collected at regular intervals (time steps),
##' and that prevalence is estimated from a logistic growth model fit to all samples available.
##'
##' @param r_v1 logistic growth rate
##' @param p0_v1 initial variant prevalence when sampling begins (# introductions / infected population size)
##' @param n per-timestep (e.g., per day) sample size
##' @param prob desired confidence in variant prevalence estimate. Default = 0.95
##' @param precision desired precision in variant prevalence estimate, relative to the true variant prevalence at the time of estimation. Default = 0.15
##' @param omega probability of sequencing (or other characterization) success. Default = 1
##' @param c_ratio coefficient of detection ratio, calculated as the ratio of the coefficients of variant 1 to variant 2. Default = 1 (no bias)
##' @param tau turnaround time, i.e., the number of time steps between sample collection and when results are available. Default = 0
##' @param max_t maximum number of time steps to consider. Default = 10000
##' @return scalar giving the time step at which the desired precision is first reached (NA if not reached by `max_t`)
##'
##' @author B. K. M. Case
##'
##' @examples
##' vartrack_tmin_prev(r_v1 = 0.12, p0_v1 = 1/10000, n = 100, omega = 0.8)
##'
##' @family variant prevalence estimation functions
##' @family variant tracking functions
##'
##' @export


vartrack_tmin_prev <- function(r_v1, p0_v1, n, prob = 0.95, precision = 0.15, omega = 1,
    c_ratio = 1, tau = 0, max_t = 10000) {
    if (!all(is.numeric(n), n > 0))
        stop("Number of samples per timestep must be numeric and greater than zero.")
    if (!all(is.numeric(prob), prob > 0 & prob < 1))
        stop("Desired confidence must be numeric and between 0 and 1.")
    if (!all(is.numeric(precision), precision > 0 & precision < 1))
        stop("Desired precision must be numeric and between 0 and 1.")
    if (!all(is.numeric(p0_v1), p0_v1 > 0 & p0_v1 < 1))
        stop("Initial variant prevalence must be numeric and between 0 and 1.")
    if (!all(is.numeric(omega), omega > 0 & omega <= 1))
        stop("Sequencing success rate must be numeric and between 0 and 1.")
    if (!all(is.numeric(r_v1), r_v1 != 0))
        stop("Growth rate must be numeric and non-zero.")
    if (!all(is.numeric(c_ratio), c_ratio > 0))
        stop("Coefficient of detection ratio must be numeric and greater than 0.")
    if (!all(is.numeric(tau), tau >= 0))
        stop("Turnaround time must be numeric and non-negative.")
    if (!all(is.numeric(max_t), max_t > tau + 1))
        stop("Maximum time step must be numeric and greater than tau + 1.")

    n_star <- n * omega
    p0_star <- varfreq_obs_freq(p0_v1, c_ratio)
    Z <- stats::qnorm(1 - ((1 - prob)/2))

    # standard error of the prevalence estimate for each candidate time step
    ts <- 0:max_t
    q <- stats::plogis(r_v1 * ts - log((1/p0_star) - 1))
    info <- cumsum(ts^2 * q * (1 - q))
    t_cand <- (floor(tau) + 1):max_t
    se <- sqrt(1/(n_star * info[t_cand - floor(tau) + 1])) * dprev_dgrowth(t_cand, p0_v1, r_v1)

    p_true <- varfreq_freq_logistic(t_cand, p0_v1, r_v1, c_ratio = 1)
    reached <- which(Z * se <= precision * p_true)
    if (length(reached) == 0) {
        warning("Desired precision not reached by max_t. This likely means n is too small to reach the desired precision in finite time. Returning NA")
        return(NA)
    }
    return(t_cand[reached[1]])
}
