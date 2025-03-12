library(testthat)
library(data.table)
library(dplyr)

modelhier_A <- readRDS(here::here("outputs", "fits", "fudan_e3", "base_hier", "model_summary.rds"))

# Sample test data
test_df <- data.table(
  id = c(1, 1, 2, 2),
  sample = c(1, 2, 1, 2),
  biomarker = c("X", "Y", "X", "Y"),
  time = c(1, 2, 1, 2)
)

# Test filter_and_sort_data
test_that("filter_and_sort_data works correctly", {
  result <- filter_and_sort_data(test_df, 1, 1, "X")
  expect_equal(nrow(result), 1)
  expect_equal(result$time, 1)
})

# Sample model outline for hierarchical parameters
test_model_outline <- list(
  abkineticsModel = list(
    "key1" = list(
      hierFlag = TRUE,
      dataHier = c(1, 2, 3)
    ),
    "key2" = list(
    )
  )
)

# Test extract_hierarchical_params
test_that("extract_hierarchical_params extracts correctly", {
    result <- extract_hierarchical_params_data(test_model_outline, "key1", 2)
    expect_equal(result, 2)

    result <- extract_hierarchical_params_data(test_model_outline, "key2", 2)
    expect_equal(result, 1)


    extract_hierarchical_params_data(modelhier_A, 2, 2)

})

# Sample model summary for adjusting parameters
test_model_summary <- list(
  fit = list(
    model = list(
      infoModel = list(
        logitBoundaries = data.frame(
          par_name = "param1",
          lb = 0,
          ub = 1
        )
      )
    )
  )
)

test_post_fit <- data.frame(param1 = 0.5)

# Test adjust_hierarchical_params
test_that("adjust_hierarchical_params modifies parameters correctly", {
  result <- adjust_hierarchical_params(test_post_fit, list("param1"), "param1", "param1", test_model_summary)
  expect_true("param1" %in% colnames(result))
})

# Dummy function for calculating trajectory
dummy_ab_func <- function(titre, t, params) {
  return(titre + t * params[1])
}

# Test calculate_trajectory
test_that("calculate_trajectory computes correctly", {
  result <- calculate_trajectory(10, 5, dummy_ab_func, c(2))
  expect_equal(length(result), 5)
  expect_equal(result[1], 12)
  expect_equal(result[5], 20)
})
