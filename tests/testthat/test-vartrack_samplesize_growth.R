
test_that("vartrack_samplesize_growth input arguments are valid", {

    expect_error(vartrack_samplesize_growth(prob = 0.95, precision = 0.15, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.84), "must be provided")

    expect_error(vartrack_samplesize_growth(prob = 0.95, precision = 0.15, t = 8, p_v1 = 0.03,
        omega = 0.9, p0_v1 = 1/10000, r_v1 = 0.84), "either")

    expect_error(vartrack_samplesize_growth(prob = 1.1, precision = 0.15, t = 8, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.84))

    expect_error(vartrack_samplesize_growth(prob = 0.95, precision = 0, t = 8, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.84))

    expect_error(vartrack_samplesize_growth(prob = 0.95, precision = 0.15, t = "invalid",
        omega = 0.9, p0_v1 = 1/10000, r_v1 = 0.84))

    expect_error(vartrack_samplesize_growth(prob = 0.95, precision = 0.15, p_v1 = 1.2,
        omega = 0.9, p0_v1 = 1/10000, r_v1 = 0.84))

    expect_error(vartrack_samplesize_growth(prob = 0.95, precision = 0.15, t = 8, omega = 1.1,
        p0_v1 = 1/10000, r_v1 = 0.84))

    expect_error(vartrack_samplesize_growth(prob = 0.95, precision = 0.15, t = 8, omega = 0.9,
        p0_v1 = 0, r_v1 = 0.84))

    expect_error(vartrack_samplesize_growth(prob = 0.95, precision = 0.15, t = 8, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0))

    expect_error(vartrack_samplesize_growth(prob = 0.95, precision = 0.15, t = 8, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.84, c_ratio = 0))

    expect_error(vartrack_samplesize_growth(prob = 0.95, precision = 0.15, t = 8, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.84, tau = -1))

    expect_error(vartrack_samplesize_growth(prob = 0.95, precision = 0.15, t = 3, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.84, tau = 3), "tau")
})

test_that("vartrack_samplesize_growth return object is valid double", {

    expect_type(vartrack_samplesize_growth(prob = 0.95, precision = 0.15, t = 8, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.84), "double")
})

test_that("vartrack_samplesize_growth sample size increases with sequencing failure, bias toward other variants, and turnaround time", {

    n <- vartrack_samplesize_growth(prob = 0.95, precision = 0.15, t = 8, omega = 1,
        p0_v1 = 1/10000, r_v1 = 0.84)

    expect_equal(vartrack_samplesize_growth(prob = 0.95, precision = 0.15, t = 8, omega = 0.5,
        p0_v1 = 1/10000, r_v1 = 0.84), 2 * n)
    expect_gt(vartrack_samplesize_growth(prob = 0.95, precision = 0.15, t = 8, omega = 1,
        p0_v1 = 1/10000, r_v1 = 0.84, c_ratio = 0.8), n)
    expect_gt(vartrack_samplesize_growth(prob = 0.95, precision = 0.15, t = 8, omega = 1,
        p0_v1 = 1/10000, r_v1 = 0.84, tau = 2), n)
})

test_that("vartrack_samplesize_growth prevalence threshold matches equivalent time step", {

    expect_equal(vartrack_samplesize_growth(prob = 0.95, precision = 0.15, p_v1 = 0.03,
        omega = 1, p0_v1 = 1/10000, r_v1 = 0.84),
        vartrack_samplesize_growth(prob = 0.95, precision = 0.15, t = 6, omega = 1,
        p0_v1 = 1/10000, r_v1 = 0.84))
})

test_that("vartrack_samplesize_growth results match original implementation", {

    expect_equal(ceiling(vartrack_samplesize_growth(prob = 0.95, precision = 0.15, t = 8,
        omega = 1, p0_v1 = 1/10000, r_v1 = 0.84)), 35)

    expect_equal(ceiling(vartrack_samplesize_growth(prob = 0.8, precision = 0.25, t = 8,
        omega = 1, p0_v1 = 1/10000, r_v1 = 0.84)), 6)
})
