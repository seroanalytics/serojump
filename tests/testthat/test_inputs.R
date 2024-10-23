test_that("Check that data_sero has no single observations", {
    # these take a while, so don't run on CI
        data_sero <- data.frame(
            id = c(1, 1, 2, 2, 3, 3, 3, 4),
            time = c(20, 300, 31, 200, 30, 50, 100, 40),
            IgG = c(0, 4, 0, 4, 0, 2, 4, 0)
        )
        expect_error(check_sero_no_single_entries(data_sero))
    }
)

test_that("Check that data_sero has no small intervals", {
    # these take a while, so don't run on CI
        data_sero <- data.frame(
            id = c(1, 1, 2, 2, 3, 3, 3, 4),
            time = c(20, 24, 31, 200, 30, 34, 100, 40),
            IgG = c(0, 4, 0, 4, 0, 2, 4, 0)
        )
        expect_error(check_sero_timings(data_sero))
    }
)

test_that("Check that data_sero has no small intervals", {
    # these take a while, so don't run on CI
        data_sero <- data.frame(
            id = c(1, 1, 2, 2, 3, 3, 3),
            time = c(20, 300, 31, 200, 30, 34, 100),
            IgG = c(0, 4, 0, 4, 0, 2, 4)
        )

        known_exp <- data.frame(
            id = c(1,  3, 3),
            time = c(10,  32, 300),
            expsoure_type = c("vax", "inf", "vax")
        )

        expect_error(check_exposures_times(data_sero, known_exp))
    }
)
