test_that("calculate_reference_titre_expectation Check that is gives sensible predictions", {
    # these take a while, so don't run on CI
        model_summary <- readRDS(here::here("outputs", "fits", "transvir_data", "wave3_base", "model_summary.rds"))
        calculate_reference_titre_expectation(model_summary)

        expect_error(check_sero_no_single_entries(data_sero))
    }
)