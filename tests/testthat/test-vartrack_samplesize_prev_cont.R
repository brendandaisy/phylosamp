
test_that("vartrack_samplesize_prev_cont input arguments are valid", {

    expect_error(vartrack_samplesize_prev_cont(prob = 0.95, precision = 0.15, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.84), "must be provided")

    expect_error(vartrack_samplesize_prev_cont(prob = 0.95, precision = 0.15, t = 8, p_v1 = 0.03,
        omega = 0.9, p0_v1 = 1/10000, r_v1 = 0.84), "either")

    expect_error(vartrack_samplesize_prev_cont(prob = 0, precision = 0.15, t = 8, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.84))

    expect_error(vartrack_samplesize_prev_cont(prob = 0.95, precision = "invalid", t = 8,
        omega = 0.9, p0_v1 = 1/10000, r_v1 = 0.84))

    expect_error(vartrack_samplesize_prev_cont(prob = 0.95, precision = 0.15, p_v1 = 0,
        omega = 0.9, p0_v1 = 1/10000, r_v1 = 0.84))

    expect_error(vartrack_samplesize_prev_cont(prob = 0.95, precision = 0.15, t = 8, omega = 0,
        p0_v1 = 1/10000, r_v1 = 0.84))

    expect_error(vartrack_samplesize_prev_cont(prob = 0.95, precision = 0.15, t = 8, omega = 0.9,
        p0_v1 = 1.5, r_v1 = 0.84))

    expect_error(vartrack_samplesize_prev_cont(prob = 0.95, precision = 0.15, t = 8, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = "invalid"))

    expect_error(vartrack_samplesize_prev_cont(prob = 0.95, precision = 0.15, t = 8, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.84, c_ratio = "invalid"))

    expect_error(vartrack_samplesize_prev_cont(prob = 0.95, precision = 0.15, p_v1 = 0.03,
        omega = 0.9, p0_v1 = 1/10000, r_v1 = 0.84, tau = 6), "tau")
})

test_that("vartrack_samplesize_prev_cont return object is valid double", {

    expect_type(vartrack_samplesize_prev_cont(prob = 0.95, precision = 0.15, t = 8, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.84), "double")
})

test_that("vartrack_samplesize_prev_cont sample size increases with sequencing failure and turnaround time", {

    n <- vartrack_samplesize_prev_cont(prob = 0.95, precision = 0.15, t = 8, omega = 1,
        p0_v1 = 1/10000, r_v1 = 0.84)

    expect_equal(vartrack_samplesize_prev_cont(prob = 0.95, precision = 0.15, t = 8,
        omega = 0.5, p0_v1 = 1/10000, r_v1 = 0.84), 2 * n)
    expect_gt(vartrack_samplesize_prev_cont(prob = 0.95, precision = 0.15, t = 8, omega = 1,
        p0_v1 = 1/10000, r_v1 = 0.84, tau = 2), n)
})

test_that("vartrack_samplesize_prev_cont results match original implementation", {

    expect_equal(ceiling(vartrack_samplesize_prev_cont(prob = 0.95, precision = 0.15, t = 8,
        omega = 1, p0_v1 = 1/10000, r_v1 = 0.84)), 1347)

    expect_equal(ceiling(vartrack_samplesize_prev_cont(prob = 0.95, precision = 0.15,
        p_v1 = 0.03, omega = 1, p0_v1 = 1/10000, r_v1 = 0.84)), 9809)
})
