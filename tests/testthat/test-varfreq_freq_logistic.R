
test_that("varfreq_freq_logistic input arguments are valid", {

    expect_error(varfreq_freq_logistic(t = "invalid", p0_v1 = 0.1, r_v1 = 0.1, c_ratio = 1))

    expect_error(varfreq_freq_logistic(t = 5, p0_v1 = 1.1, r_v1 = 0.1, c_ratio = 1))

    expect_error(varfreq_freq_logistic(t = 5, p0_v1 = "invalid", r_v1 = 0.1, c_ratio = 1))

    expect_error(varfreq_freq_logistic(t = 5, p0_v1 = 0.1, r_v1 = "invalid", c_ratio = 1))

    expect_error(varfreq_freq_logistic(t = 5, p0_v1 = 0.1, r_v1 = 0, c_ratio = 1))

    expect_error(varfreq_freq_logistic(t = 5, p0_v1 = 0.1, r_v1 = 0.1, c_ratio = "invalid"))

    expect_error(varfreq_freq_logistic(t = 5, p0_v1 = 0.1, r_v1 = 0.1, c_ratio = 0))
})

test_that("varfreq_freq_logistic return object is valid double", {

    expect_type(varfreq_freq_logistic(t = 5, p0_v1 = 0.1, r_v1 = 0.1, c_ratio = 1),
        "double")
})

test_that("varfreq_freq_logistic variant prevalence increases when growth rate is positive", {

    expect_gt(varfreq_freq_logistic(t = 5, p0_v1 = 1/10000, r_v1 = 0.3, c_ratio = 1),
        1/10000)
})

test_that("varfreq_freq_logistic variant prevalence decreases when growth rate is negative", {

    expect_lt(varfreq_freq_logistic(t = 5, p0_v1 = 1/10000, r_v1 = -0.3, c_ratio = 1),
        1/10000)
})

test_that("varfreq_freq_logistic observed variant prevalence is higher when c_ratio is above 1", {

    expect_gt(varfreq_freq_logistic(t = 4, p0_v1 = 1/10000, r_v1 = 0.2, c_ratio = 1.2),
        varfreq_freq_logistic(t = 4, p0_v1 = 1/10000, r_v1 = 0.2, c_ratio = 1))
})

test_that("varfreq_freq_logistic observed variant prevalence is lower when c_ratio is below 1", {

    expect_lt(varfreq_freq_logistic(t = 10, p0_v1 = 1/1000, r_v1 = 0.11, c_ratio = 0.9),
        varfreq_freq_logistic(t = 10, p0_v1 = 1/1000, r_v1 = 0.11, c_ratio = 1))
})

test_that("varfreq_freq_logistic manuscript results remain valid", {

    expect_equal(round(varfreq_freq_logistic(t = 14, p0_v1 = 1/10000, r_v1 = 0.1,
        c_ratio = 1), 4), 4e-04)
})

test_that("find_tmax enforces valid input arguments", {
    expect_error(find_tmax(p_v1=1.1, p0_v1=0.1, r_v1=0.1))
    expect_error(find_tmax(p_v1="invalid", p0_v1=0.1, r_v1=0.1))
    
    expect_error(find_tmax(p_v1=0.1, p0_v1=1.1, r_v1=-0.1))
    expect_error(find_tmax(p_v1=0.1, p0_v1=0, r_v1=0.1))
    expect_error(find_tmax(p_v1=0.1, p0_v1="invalid", r_v1=-0.1))
    
    expect_error(find_tmax(p_v1=0.5, p0_v1=0.1, r_v1="invalid"))
})

test_that("find_tmax returns Inf when desired prevalence cannot be reached", {
    expect_equal(find_tmax(p_v1=1, p0_v1=0.1, r_v1=0.1), Inf)
    expect_equal(find_tmax(p_v1=0.1, p0_v1=0.2, r_v1=0.1), Inf)
    expect_equal(find_tmax(p_v1=0.2, p0_v1=0.1, r_v1=-0.1), Inf)
})

test_that("find_tmax returns 0 when desired prevalence and initial prevalence are equal", {
    expect_equal(find_tmax(p_v1=0.1, p0_v1=0.1, r_v1=0.1), 0)
    expect_equal(find_tmax(p_v1=0.1, p0_v1=0.1, r_v1=-2), 0)
})

test_that("find_tmax matches with varfreq_freq_logistic", {
    tmax_curve <- min(which(varfreq_freq_logistic(0:100, 1e-4, 0.11) >= 0.2)) - 1
    expect_equal(find_tmax(0.2, 1e-4, 0.11), tmax_curve)
})

test_that("dom_time enforces valid input arguments", {
    expect_error(dom_time(p0_v1=1.1, r_v1=-0.1))
    expect_error(dom_time(p0_v1="invalid", r_v1=1.1))
    
    expect_error(dom_time(p_v1=0.5, p0_v1=0.1, r_v1="invalid"))
})

test_that("dom_time matches with varfreq_freq_logistic", {
    expect_equal(dom_time(0.6, 0.1), Inf)
    expect_equal(dom_time(0.4, -0.1), Inf)
    
    dt_curve <- min(which(varfreq_freq_logistic(0:100, 1e-4, 0.11) >= 0.5)) - 1
    expect_equal(dom_time(1e-4, 0.11), dt_curve)
})
