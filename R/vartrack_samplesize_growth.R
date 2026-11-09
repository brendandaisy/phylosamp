##' Calculate sample size needed for estimating variant growth advantage assuming periodic sampling
##'
##' This function calculates the per-timestep sample size needed for estimating the logistic
##' growth rate (growth advantage) of a variant given a desired confidence and desired precision
##' in the growth rate estimate, and either a desired maximum time by which to estimate the
##' growth rate or a desired prevalence by which to estimate it.
##' It assumes a periodic sampling strategy, where samples are collected at regular intervals (time steps),
##' and that the growth rate is estimated by maximum likelihood from all samples available.
##'
##' @param r_v1 logistic growth rate
##' @param p0_v1 initial variant prevalence when sampling begins (# introductions / infected population size)
##' @param t time step number (e.g., days) at which the growth rate should be estimated by. Default = NA (either `'t'` or `'p_v1'` should be provided, not both)
##' @param p_v1 the true variant prevalence by which the growth rate should be estimated. Sampling is assumed to occur up to the last time step before the variant's true prevalence reaches this value, regardless of `c_ratio`. Default = NA (either `'t'` or `'p_v1'` should be provided, not both)
##' @param prob desired confidence in growth rate estimate. Default = 0.95
##' @param precision desired precision in growth rate estimate, relative to the growth rate (e.g., 0.15 for an estimate within 15% of `r_v1`). Default = 0.15
##' @param omega probability of sequencing (or other characterization) success. Default = 1
##' @param c_ratio coefficient of detection ratio, calculated as the ratio of the coefficients of variant 1 to variant 2. Default = 1 (no bias)
##' @param tau turnaround time, i.e., the number of time steps between sample collection and when results are available. Default = 0
##' @return scalar of expected sample size per timestep
##'
##' @author B. K. M. Case
##'
##' @examples
##' vartrack_samplesize_growth(r_v1 = 0.84, p0_v1 = 1/10000, t = 8, omega = 0.8)
##'
##' vartrack_samplesize_growth(r_v1 = 0.84, p0_v1 = 1/10000, p_v1 = 0.03, prob = 0.95,
##' precision = 0.15, omega = 0.8, c_ratio = 1.37, tau = 1)
##'
##' @family variant growth estimation functions
##' @family variant tracking functions
##'
##' @export


vartrack_samplesize_growth <- function(r_v1, p0_v1, t = NA, p_v1 = NA, prob = 0.95, precision = 0.15,
    omega = 1, c_ratio = 1, tau = 0) {
    if (is.na(t) & is.na(p_v1))
        stop("Either time of estimation or estimation prevalence must be provided.")
    if (!is.na(t) & !is.na(p_v1))
        stop("Please provide either a desired time of estimation OR desired prevalence.")

    if (!is.na(t)) {
        if (!is.numeric(t))
            stop("Time step must be numeric.")
    }

    if (!is.na(p_v1)) {
        if (!all(is.numeric(p_v1), p_v1 > 0 & p_v1 < 1))
            stop("Variant prevalence must be numeric and between 0 and 1.")
    }

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

    if (is.na(t))
        t <- floor(varfreq_time_logistic(p_v1, p0_v1, r_v1, c_ratio = 1))

    if (t - tau < 1)
        stop("At least two time steps of results must be available before the time of estimation (t - tau >= 1).")

    p0_star <- varfreq_obs_freq(p0_v1, c_ratio)
    v <- var_growth_rate(t - tau, p0_star, r_v1)
    Z <- stats::qnorm(1 - ((1 - prob)/2))
    n <- ((Z/(precision * abs(r_v1)))^2) * v
    n_samples <- n/omega
    return(n_samples)
}
