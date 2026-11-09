##' Calculate sample size needed to detect a variant growth advantage assuming periodic sampling
##'
##' This function calculates the per-timestep sample size needed to detect a non-zero logistic
##' growth rate (growth advantage) of a variant with a desired power, using a one-sided Wald test
##' of the null hypothesis of no growth advantage, and either a desired time or prevalence by which
##' to perform the test.
##' The test is one-sided in the direction of `r_v1` (i.e., testing for growth when `r_v1 > 0` and
##' for decline when `r_v1 < 0`).
##' It assumes a periodic sampling strategy, where samples are collected at regular intervals (time steps).
##'
##' @param r_v1 logistic growth rate under the alternative hypothesis
##' @param p0_v1 initial variant prevalence when sampling begins (# introductions / infected population size)
##' @param t time step number (e.g., days) at which the test should be performed. Default = NA (either `'t'` or `'p_v1'` should be provided, not both)
##' @param p_v1 the true variant prevalence by which the test should be performed. Sampling is assumed to occur up to the last time step before the variant's true prevalence reaches this value, regardless of `c_ratio`. Default = NA (either `'t'` or `'p_v1'` should be provided, not both)
##' @template power
##' @template alpha
##' @param omega probability of sequencing (or other characterization) success. Default = 1
##' @param c_ratio coefficient of detection ratio, calculated as the ratio of the coefficients of variant 1 to variant 2. Default = 1 (no bias)
##' @param tau turnaround time, i.e., the number of time steps between sample collection and when results are available. Default = 0
##' @return scalar of expected sample size per timestep
##'
##' @author B. K. M. Case
##'
##' @examples
##' vartrack_samplesize_growth_power(r_v1 = 0.05, p0_v1 = 1/10000, t = 30, omega = 0.8)
##'
##' vartrack_samplesize_growth_power(r_v1 = 0.05, p0_v1 = 1/10000, p_v1 = 0.01, power = 0.9,
##' alpha = 0.01, omega = 0.8, c_ratio = 1.37, tau = 3)
##'
##' @family variant growth estimation functions
##' @family variant tracking functions
##'
##' @export


vartrack_samplesize_growth_power <- function(r_v1, p0_v1, t = NA, p_v1 = NA, power = 0.8, alpha = 0.05,
    omega = 1, c_ratio = 1, tau = 0) {
    if (is.na(t) & is.na(p_v1))
        stop("Either time of testing or testing prevalence must be provided.")
    if (!is.na(t) & !is.na(p_v1))
        stop("Please provide either a desired time of testing OR desired prevalence.")

    if (!is.na(t)) {
        if (!is.numeric(t))
            stop("Time step must be numeric.")
    }

    if (!is.na(p_v1)) {
        if (!all(is.numeric(p_v1), p_v1 > 0 & p_v1 < 1))
            stop("Variant prevalence must be numeric and between 0 and 1.")
    }

    if (!all(is.numeric(power), power > 0 & power < 1))
        stop("Desired power must be numeric and between 0 and 1.")
    if (!all(is.numeric(alpha), alpha > 0 & alpha < 1))
        stop("Alpha level must be numeric and between 0 and 1.")
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
        stop("At least two time steps of results must be available before the time of testing (t - tau >= 1).")

    p0_star <- varfreq_obs_freq(p0_v1, c_ratio)
    v <- var_growth_rate(t - tau, p0_star, r_v1)
    Z_a <- stats::qnorm(1 - alpha)
    Z_b <- stats::qnorm(power)
    n <- (v/r_v1^2) * (Z_a + Z_b)^2
    n_samples <- n/omega
    return(n_samples)
}
