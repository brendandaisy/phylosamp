
test_that("vartrack_samplesize_growth_power input arguments are valid", {

    expect_error(vartrack_samplesize_growth_power(power = 0.8, omega = 0.9, p0_v1 = 1/10000,
        r_v1 = 0.05), "must be provided")

    expect_error(vartrack_samplesize_growth_power(power = 0.8, t = 30, p_v1 = 0.01, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.05), "either")

    expect_error(vartrack_samplesize_growth_power(power = 1, t = 30, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.05))

    expect_error(vartrack_samplesize_growth_power(power = 0.8, t = 30, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.05, alpha = 0))

    expect_error(vartrack_samplesize_growth_power(power = 0.8, t = 30, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0))

    expect_error(vartrack_samplesize_growth_power(power = 0.8, t = 30, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.05, c_ratio = 0))

    expect_error(vartrack_samplesize_growth_power(power = 0.8, t = 3, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.05, tau = 3), "tau")
})

test_that("vartrack_samplesize_growth_power return object is valid double", {

    expect_type(vartrack_samplesize_growth_power(power = 0.8, t = 30, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.05), "double")
})

test_that("vartrack_samplesize_growth_power is the inverse of vartrack_power_growth", {

    n <- vartrack_samplesize_growth_power(power = 0.8, t = 30, omega = 0.7, p0_v1 = 1/10000,
        r_v1 = 0.05, c_ratio = 1.37, tau = 2)

    expect_equal(vartrack_power_growth(n = n, t = 30, omega = 0.7, p0_v1 = 1/10000,
        r_v1 = 0.05, c_ratio = 1.37, tau = 2), 0.8)
})

test_that("vartrack_samplesize_growth_power results match original implementation", {

    expect_equal(ceiling(vartrack_samplesize_growth_power(power = 0.8, t = 30, omega = 1,
        p0_v1 = 1/10000, r_v1 = 0.05)), 802)

    expect_equal(ceiling(vartrack_samplesize_growth_power(power = 0.9, t = 30, omega = 1,
        p0_v1 = 1/10000, r_v1 = 0.05, alpha = 0.01)), 1688)
})
