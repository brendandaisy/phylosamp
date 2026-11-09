test_that("vartrack_samplesize_prev input arguments are valid", {

    expect_error(vartrack_samplesize_prev(p_v1 = 0.1, prob = 0.9, precision = 0.1, omega = 0.9, c_ratio = 1.5, sampling_freq = "fake"),
        "sampling frequency argument")
    expect_error(vartrack_samplesize_prev(p_v1 = 0.1, prob = 0.9, precision = 0.1, omega = 0.9, c_ratio = 1.5, sampling_freq = "cont"),
        "correct arguments for the cont method")
    expect_error(vartrack_samplesize_prev(p_v1 = 0.1, t = 8, prob = 0.9, precision = 0.1, omega = 0.9, p0_v1 = 1/10000, r_v1 = 0.12, c_ratio = 1.5, sampling_freq = "cont"),
        "correct arguments for the cont method")
})

test_that("vartrack_samplesize_prev return object is valid double", {

    expect_type(vartrack_samplesize_prev(p_v1 = 0.1, prob = 0.9, omega = 0.9, precision = 0.1, c_ratio = 1.5, sampling_freq = "xsect"), "double")
    expect_type(vartrack_samplesize_prev(p_v1 = 0.03, prob = 0.9, omega = 0.9, precision = 0.1, p0_v1 = 1/10000, r_v1 = 0.12, c_ratio = 1.5, sampling_freq = "cont"), "double")
})

