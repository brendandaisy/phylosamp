##' Calculate sample size needed for estimating variant prevalence given a desired confidence
##'
##' This function calculates the sample size needed for estimating variant prevalence
##' given a desired confidence and desired precision in the variant prevalence estimate
##' and sampling strategy.
##'
##' @param p_v1 variant prevalence (proportion). For periodic sampling, the true variant prevalence (regardless of `c_ratio`) at which variant prevalence should be estimated. Default = NA (for periodic sampling, either `'t'` or `'p_v1'` should be provided, not both)
##' @param prob desired confidence in variant prevalence estimate
##' @param precision desired precision in variant prevalence estimate
##' @param omega probability of sequencing (or other characterization) success
##' @param c_ratio coefficient of detection ratio, calculated as the ratio of the coefficients of variant 1 to variant 2. Default = 1 (no bias)
##' @param sampling_freq the sampling frequency (must be either 'xsect' or 'cont')
##' @param t time step number (e.g., days) at which prevalence should be estimated (periodic sampling only). Default = NA
##' @param p0_v1 initial variant prevalence when sampling begins (periodic sampling only)
##' @param r_v1 logistic growth rate (periodic sampling only)
##' @param tau turnaround time, i.e., the number of time steps between sample collection and when results are available (periodic sampling only). Default = 0
##' @return scalar of sample size (per timestep for periodic sampling)
##'
##' @author Shirlee Wohl, Elizabeth C. Lee, Bethany L. DiPrete, and Justin Lessler
##'
##' @examples
##' # Cross-sectional sampling
##' vartrack_samplesize_prev(p_v1 = 0.1, prob = 0.95, precision = 0.25,
##' omega = 0.8, c_ratio = 1, sampling_freq = 'xsect')
##'
##' # Periodic sampling
##' vartrack_samplesize_prev(p_v1 = 0.03, prob = 0.95, precision = 0.15, omega = 0.8,
##' p0_v1 = 1/10000, r_v1 = 0.12, c_ratio = 1, sampling_freq = 'cont')
##'
##' @family variant prevalence estimation functions
##' @family variant tracking functions
##'
##' @export


vartrack_samplesize_prev <- function(p_v1 = NA, prob, precision, omega, c_ratio = 1, sampling_freq,
    t = NA, p0_v1 = NA, r_v1 = NA, tau = 0) {

    if (sampling_freq == "xsect") {
        if(any(is.na(c(p_v1, prob, precision, omega, c_ratio)))){
            stop("Please specify the correct arguments for the xsect method: p_v1, prob, precision, omega, c_ratio.")
        }

        message("Calculating sample size for variant prevalence estimation assuming single cross-sectional sample")
        out <- vartrack_samplesize_prev_xsect(p_v1 = p_v1, prob = prob, precision = precision,
            omega = omega, c_ratio = c_ratio)

    } else if (sampling_freq == "cont") {
        if(any(is.na(c(prob, precision, omega, p0_v1, r_v1, c_ratio, tau))) | sum(is.na(c(t, p_v1))) %in% c(0,2)) {
            stop("Please specify the correct arguments for the cont method: prob, precision, t, p_v1, omega, p0_v1, r_v1, c_ratio, tau.")
        }

        message("Calculating sample size for variant prevalence estimation assuming periodic sampling")
        out <- vartrack_samplesize_prev_cont(prob = prob, precision = precision, t = t, p_v1 = p_v1,
            omega = omega, p0_v1 = p0_v1, r_v1 = r_v1, c_ratio = c_ratio, tau = tau)

    } else {
        stop("Incorrect sampling frequency argument (please specify 'xsect' or 'cont')")
    }
    return(out)
}
