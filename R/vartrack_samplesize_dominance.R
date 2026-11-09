##' Calculate sample size needed for estimating variant dominance time assuming periodic sampling
##'
##' This function calculates the per-timestep sample size needed for estimating the time at which a
##' variant becomes dominant (reaches 50% true prevalence) given a desired confidence and desired precision,
##' and either a desired time or prevalence by which to estimate the dominance time.
##' It assumes a periodic sampling strategy, where samples are collected at regular intervals (time steps),
##' and that the dominance time is estimated from a logistic growth model fit to all samples available.
##' Precision is expressed in terms of prevalence: an estimate with precision 0.15 falls within the time
##' window during which the variant is between 0.5*(1 - 0.15) and 0.5*(1 + 0.15) true prevalence.
##'
##' @param r_v1 logistic growth rate
##' @param p0_v1 initial variant prevalence when sampling begins (# introductions / infected population size)
##' @param t time step number (e.g., days) at which the dominance time should be estimated by. Default = NA (either `'t'` or `'p_v1'` should be provided, not both)
##' @param p_v1 the true variant prevalence by which the dominance time should be estimated. Sampling is assumed to occur up to the last time step before the variant's true prevalence reaches this value, regardless of `c_ratio`. Default = NA (either `'t'` or `'p_v1'` should be provided, not both)
##' @param prob desired confidence in dominance time estimate. Default = 0.95
##' @param precision desired precision in dominance time estimate, relative to the 50% true variant prevalence at dominance. Default = 0.15
##' @param omega probability of sequencing (or other characterization) success. Default = 1
##' @param c_ratio coefficient of detection ratio, calculated as the ratio of the coefficients of variant 1 to variant 2. Default = 1 (no bias)
##' @param tau turnaround time, i.e., the number of time steps between sample collection and when results are available. Default = 0
##' @return scalar of expected sample size per timestep
##'
##' @author B. K. M. Case
##'
##' @examples
##' vartrack_samplesize_dominance(r_v1 = 0.84, p0_v1 = 1/10000, t = 8, omega = 0.8)
##'
##' vartrack_samplesize_dominance(r_v1 = 0.84, p0_v1 = 1/10000, p_v1 = 0.03, prob = 0.95,
##' precision = 0.15, omega = 0.8, c_ratio = 1.37, tau = 1)
##'
##' @family variant dominance time estimation functions
##' @family variant tracking functions
##'
##' @export


vartrack_samplesize_dominance <- function(r_v1, p0_v1, t = NA, p_v1 = NA, prob = 0.95, precision = 0.15,
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
    if (varfreq_time_logistic(0.5, p0_v1, r_v1, c_ratio = 1) <= 0)
        stop("Variant must reach 50% prevalence after sampling begins (p0_v1 < 0.5 with r_v1 > 0, or p0_v1 > 0.5 with r_v1 < 0).")

    if (is.na(t))
        t <- floor(varfreq_time_logistic(p_v1, p0_v1, r_v1, c_ratio = 1))

    if (t - tau < 1)
        stop("At least two time steps of results must be available before the time of estimation (t - tau >= 1).")

    p0_star <- varfreq_obs_freq(p0_v1, c_ratio)
    v <- var_growth_rate(t - tau, p0_star, r_v1) * ddominance_dgrowth(p0_v1, r_v1)^2
    Z <- stats::qnorm(1 - ((1 - prob)/2))
    n <- ((Z/dominance_error(precision, r_v1))^2) * v
    n_samples <- n/omega
    return(n_samples)
}
