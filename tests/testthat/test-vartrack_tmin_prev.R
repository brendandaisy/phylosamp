
test_that("vartrack_tmin_prev input arguments are valid", {

    expect_error(vartrack_tmin_prev(n = "invalid", prob = 0.95, precision = 0.15, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.12))

    expect_error(vartrack_tmin_prev(n = 100, prob = 0.95, precision = 0.15, omega = 0,
        p0_v1 = 1/10000, r_v1 = 0.12))

    expect_error(vartrack_tmin_prev(n = 100, prob = 0.95, precision = 0.15, omega = 0.9,
        p0_v1 = 1, r_v1 = 0.12))

    expect_error(vartrack_tmin_prev(n = 100, prob = 0.95, precision = 0.15, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.12, c_ratio = -1))

    expect_error(vartrack_tmin_prev(n = 100, prob = 0.95, precision = 0.15, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.12, tau = -1))
})

test_that("vartrack_tmin_prev return object is valid integer", {

    expect_type(vartrack_tmin_prev(n = 100, prob = 0.95, precision = 0.15, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.12), "integer")
})

test_that("vartrack_tmin_prev returns NA with a warning when precision is unreachable", {

    expect_warning(out <- vartrack_tmin_prev(n = 1, prob = 0.95, precision = 0.01,
        omega = 1, p0_v1 = 1/10000, r_v1 = 0.12, max_t = 50))
    expect_true(is.na(out))
})

test_that("vartrack_tmin_prev is the inverse of vartrack_samplesize_prev_cont", {

    t <- vartrack_tmin_prev(n = 100, prob = 0.95, precision = 0.15, omega = 1,
        p0_v1 = 1/10000, r_v1 = 0.12, tau = 2)

    expect_lte(vartrack_samplesize_prev_cont(prob = 0.95, precision = 0.15, t = t, omega = 1,
        p0_v1 = 1/10000, r_v1 = 0.12, tau = 2), 100)
    expect_gt(vartrack_samplesize_prev_cont(prob = 0.95, precision = 0.15, t = t - 1, omega = 1,
        p0_v1 = 1/10000, r_v1 = 0.12, tau = 2), 100)
})
