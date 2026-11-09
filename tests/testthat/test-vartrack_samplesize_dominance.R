
test_that("vartrack_samplesize_dominance input arguments are valid", {

    expect_error(vartrack_samplesize_dominance(prob = 0.95, precision = 0.15, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.84), "must be provided")

    expect_error(vartrack_samplesize_dominance(prob = 0.95, precision = 0.15, t = 8, p_v1 = 0.03,
        omega = 0.9, p0_v1 = 1/10000, r_v1 = 0.84), "either")

    expect_error(vartrack_samplesize_dominance(prob = -1, precision = 0.15, t = 8, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.84))

    expect_error(vartrack_samplesize_dominance(prob = 0.95, precision = 1, t = 8, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.84))

    expect_error(vartrack_samplesize_dominance(prob = 0.95, precision = 0.15, t = 8, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0))

    expect_error(vartrack_samplesize_dominance(prob = 0.95, precision = 0.15, t = 8, omega = 0.9,
        p0_v1 = 0.6, r_v1 = 0.84), "50%")

    expect_error(vartrack_samplesize_dominance(prob = 0.95, precision = 0.15, t = 8, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = -0.84), "50%")

    expect_error(vartrack_samplesize_dominance(prob = 0.95, precision = 0.15, t = 2, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.84, tau = 2), "tau")
})

test_that("vartrack_samplesize_dominance return object is valid double", {

    expect_type(vartrack_samplesize_dominance(prob = 0.95, precision = 0.15, t = 8, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.84), "double")
})

test_that("vartrack_samplesize_dominance supports declining variants", {

    expect_type(vartrack_samplesize_dominance(prob = 0.95, precision = 0.15, t = 8, omega = 0.9,
        p0_v1 = 0.9, r_v1 = -0.3), "double")
})

test_that("vartrack_samplesize_dominance sample size increases with sequencing failure and turnaround time", {

    n <- vartrack_samplesize_dominance(prob = 0.95, precision = 0.15, t = 8, omega = 1,
        p0_v1 = 1/10000, r_v1 = 0.84)

    expect_equal(vartrack_samplesize_dominance(prob = 0.95, precision = 0.15, t = 8,
        omega = 0.5, p0_v1 = 1/10000, r_v1 = 0.84), 2 * n)
    expect_gt(vartrack_samplesize_dominance(prob = 0.95, precision = 0.15, t = 8, omega = 1,
        p0_v1 = 1/10000, r_v1 = 0.84, tau = 2), n)
})

test_that("vartrack_samplesize_dominance results match original implementation", {

    expect_equal(ceiling(vartrack_samplesize_dominance(prob = 0.95, precision = 0.15, t = 8,
        omega = 1, p0_v1 = 1/10000, r_v1 = 0.84)), 731)

    expect_equal(ceiling(vartrack_samplesize_dominance(prob = 0.8, precision = 0.25, t = 8,
        omega = 1, p0_v1 = 1/10000, r_v1 = 0.84)), 110)
})
