
test_that("vartrack_tmin_growth input arguments are valid", {

    expect_error(vartrack_tmin_growth(n = 0, prob = 0.95, precision = 0.15, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.12))

    expect_error(vartrack_tmin_growth(n = 100, prob = 1, precision = 0.15, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.12))

    expect_error(vartrack_tmin_growth(n = 100, prob = 0.95, precision = -0.15, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.12))

    expect_error(vartrack_tmin_growth(n = 100, prob = 0.95, precision = 0.15, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0))

    expect_error(vartrack_tmin_growth(n = 100, prob = 0.95, precision = 0.15, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.12, tau = 5, max_t = 5))
})

test_that("vartrack_tmin_growth return object is valid integer", {

    expect_type(vartrack_tmin_growth(n = 100, prob = 0.95, precision = 0.15, omega = 0.9,
        p0_v1 = 1/10000, r_v1 = 0.12), "integer")
})

test_that("vartrack_tmin_growth returns NA with a warning when precision is unreachable", {

    expect_warning(out <- vartrack_tmin_growth(n = 1, prob = 0.95, precision = 0.01,
        omega = 1, p0_v1 = 1/10000, r_v1 = 0.12, max_t = 50))
    expect_true(is.na(out))
})

test_that("vartrack_tmin_growth is the inverse of vartrack_samplesize_growth", {

    t <- vartrack_tmin_growth(n = 100, prob = 0.95, precision = 0.15, omega = 1,
        p0_v1 = 1/10000, r_v1 = 0.12, tau = 2)

    expect_lte(vartrack_samplesize_growth(prob = 0.95, precision = 0.15, t = t, omega = 1,
        p0_v1 = 1/10000, r_v1 = 0.12, tau = 2), 100)
    expect_gt(vartrack_samplesize_growth(prob = 0.95, precision = 0.15, t = t - 1, omega = 1,
        p0_v1 = 1/10000, r_v1 = 0.12, tau = 2), 100)
})
