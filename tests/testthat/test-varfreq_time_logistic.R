
test_that("varfreq_time_logistic input arguments are valid", {

    expect_error(varfreq_time_logistic(p_v1 = "invalid", p0_v1 = 0.1, r_v1 = 0.1, c_ratio = 1))

    expect_error(varfreq_time_logistic(p_v1 = 1, p0_v1 = 0.1, r_v1 = 0.1, c_ratio = 1))

    expect_error(varfreq_time_logistic(p_v1 = 0.5, p0_v1 = 1.1, r_v1 = 0.1, c_ratio = 1))

    expect_error(varfreq_time_logistic(p_v1 = 0.5, p0_v1 = "invalid", r_v1 = 0.1, c_ratio = 1))

    expect_error(varfreq_time_logistic(p_v1 = 0.5, p0_v1 = 0.1, r_v1 = "invalid", c_ratio = 1))

    expect_error(varfreq_time_logistic(p_v1 = 0.5, p0_v1 = 0.1, r_v1 = 0, c_ratio = 1))

    expect_error(varfreq_time_logistic(p_v1 = 0.5, p0_v1 = 0.1, r_v1 = 0.1, c_ratio = 0))
})

test_that("varfreq_time_logistic return object is valid double", {

    expect_type(varfreq_time_logistic(p_v1 = 0.5, p0_v1 = 0.1, r_v1 = 0.1, c_ratio = 1),
        "double")
})

test_that("varfreq_time_logistic is the inverse of varfreq_freq_logistic", {

    t <- varfreq_time_logistic(p_v1 = 0.03, p0_v1 = 1/10000, r_v1 = 0.12, c_ratio = 1.37)

    expect_equal(varfreq_freq_logistic(t = t, p0_v1 = 1/10000, r_v1 = 0.12, c_ratio = 1.37),
        0.03)

    expect_equal(varfreq_time_logistic(p_v1 = varfreq_freq_logistic(t = 14, p0_v1 = 0.6,
        r_v1 = -0.2), p0_v1 = 0.6, r_v1 = -0.2), 14)
})

test_that("varfreq_time_logistic is negative when prevalence was reached before time 0", {

    expect_lt(varfreq_time_logistic(p_v1 = 0.01, p0_v1 = 0.1, r_v1 = 0.1), 0)
})

test_that("varfreq_time_logistic variant reaches prevalence sooner when c_ratio is above 1", {

    expect_lt(varfreq_time_logistic(p_v1 = 0.01, p0_v1 = 1/10000, r_v1 = 0.1, c_ratio = 1.2),
        varfreq_time_logistic(p_v1 = 0.01, p0_v1 = 1/10000, r_v1 = 0.1, c_ratio = 1))
})
