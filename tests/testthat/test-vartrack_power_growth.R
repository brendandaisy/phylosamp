
test_that("vartrack_power_growth input arguments are valid", {

    expect_error(vartrack_power_growth(n = 100, omega = 0.9, p0_v1 = 1/10000, r_v1 = 0.05),
        "must be provided")

    expect_error(vartrack_power_growth(n = 0, t = 30, omega = 0.9, p0_v1 = 1/10000,
        r_v1 = 0.05))

    expect_error(vartrack_power_growth(n = 100, t = 30, omega = 0.9, p0_v1 = 1/10000,
        r_v1 = 0.05, alpha = 1))

    expect_error(vartrack_power_growth(n = 100, t = 30, omega = 0, p0_v1 = 1/10000,
        r_v1 = 0.05))

    expect_error(vartrack_power_growth(n = 100, t = 30, omega = 0.9, p0_v1 = 1/10000,
        r_v1 = 0.05, tau = -1))
})

test_that("vartrack_power_growth return object is valid probability", {

    p <- vartrack_power_growth(n = 100, t = 30, omega = 0.9, p0_v1 = 1/10000, r_v1 = 0.05)

    expect_type(p, "double")
    expect_gte(p, 0)
    expect_lte(p, 1)
})

test_that("vartrack_power_growth power increases with sample size and is symmetric in growth direction", {

    expect_gt(vartrack_power_growth(n = 200, t = 30, omega = 1, p0_v1 = 1/10000, r_v1 = 0.05),
        vartrack_power_growth(n = 100, t = 30, omega = 1, p0_v1 = 1/10000, r_v1 = 0.05))

    expect_equal(vartrack_power_growth(n = 20, t = 10, omega = 1, p0_v1 = 0.4, r_v1 = 0.2),
        vartrack_power_growth(n = 20, t = 10, omega = 1, p0_v1 = 0.6, r_v1 = -0.2))
})

test_that("vartrack_power_growth results match original implementation", {

    expect_equal(round(vartrack_power_growth(n = 2, t = 30, omega = 1, p0_v1 = 1/10000,
        r_v1 = 0.05), 4), 0.0642)
})
