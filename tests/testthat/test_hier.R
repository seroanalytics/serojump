correctFunc <- function(titre, time, pars) {
  # A dummy body (for testing purposes)
  titre
}

test_that("errors when any argument is NULL", {
  expect_error(
    addAbkineticsModelHier(NULL, "bio", "exp", c("a", "b"), c("a"), 1:10, correctFunc),
    "One or more arguments of `addAbkineticsModelHier` are NULL."
  )
})

test_that("errors if dataHier is not a subset of pars", {
  expect_error(
    addAbkineticsModelHier("id", "bio", "exp", c("a", "b"), c("a", "c", "d"), 1:10, correctFunc),
    "All elements of parsHier must be a subset of pars."
  )
})

test_that("errors if dataHier is not a numeric vector", {
  expect_error(
    addAbkineticsModelHier("id", "bio", "exp", c("a", "b"), c("a"), letters[1:10], correctFunc),
    "dataHier must be a numeric vector."
  )
})

test_that("errors if funcForm is not a function", {
  expect_error(
    addAbkineticsModelHier("id", "bio", "exp", c("a", "b"), c("a"), 1:10, "not a function"),
    "funcForm must be a function."
  )
})


test_that("returns a proper list when all conditions are met", {
  result <- addAbkineticsModelHier("id", "bio", "exp", c("a", "b"), c("a"), c(rep(1, 5), rep(2, 5)), correctFunc)
  
  expect_type(result, "list")
  expect_equal(result$id, "id")
  expect_equal(result$biomarker, "bio")
  expect_equal(result$exposureType, "exp")
  expect_equal(result$pars, c("a", "z_a_1", "z_a_2", "sigma_a", "b"))
  expect_equal(result$parsHier, c("a"))
  expect_equal(result$dataHier, c(rep(1, 5), rep(2, 5)))
  expect_equal(result$funcForm, correctFunc)
})

test_that("returns a proper list when all conditions are met", {
  result <- addAbkineticsModelHier("id", "bio", "exp", c("a", "b"), c("a", "b"), c(rep(1, 5), rep(2, 5)), correctFunc)
  
  expect_type(result, "list")
  expect_equal(result$id, "id")
  expect_equal(result$biomarker, "bio")
  expect_equal(result$exposureType, "exp")
  expect_equal(result$pars, c("a", "z_a_1", "z_a_2", "sigma_a", "b", "z_b_1", "z_b_2", "sigma_b"))
  expect_equal(result$parsHier, c("a", "b"))
  expect_equal(result$dataHier, c(rep(1, 5), rep(2, 5)))
  expect_equal(result$funcForm, correctFunc)
})

test_that("runSeroJump produces reproducible output with hierarchical effects", {

    seed_i <- 123
    set.seed(seed_i)

    data_titre_model <- data.frame(
        id = rep(1:10, each = 10),
        time = rep(c(2, 5, 7, 12, 14, 15, 19, 30, 34, 40), 10),
        sVNT = c(100, rnorm(9, 200, 10), rnorm(40, 100, 10), 100, rnorm(9, 400, 10), rnorm(40, 100, 10)) / 100,
        age_group = rep(c("1", "2"), each = 50)
    )

    covar_key <- data_titre_model %>% select(id, age_group) %>% unique %>% pull(age_group)

    obsLogLikelihood <- function(titre_val, titre_est, pars) {
        ll <- dnorm(titre_val, titre_est, pars[1], log = TRUE)
    }

    infSerumKinetics <- function(titre_est, timeSince, pars) {
        a <- pars[1]
        titre_est <- titre_est + a
        titre_est
    }

    noInfSerumKinetics <- function(titre_est, timeSince, pars) {
        titre_est_log <- titre_est - pars[1] * (timeSince)
        titre_est_log <- max(0, titre_est_log)
        titre_est_log
    }

    # Define the biomarkers and exposure types in the model
    biomarkers <- c("sVNT")
    exposureTypes <- c("none", "inf")
    exposureFitted <- "inf"

    # Define the observational model
    observationalModel <- list(
        names = c("sVNT"),
        model = makeModel(
            addObservationalModel("sVNT", c("sigma"), obsLogLikelihood)
            ), # observational model,
        prior = bind_rows(
                addPrior("sigma", 0.0001, 100, "unif", 0.0001, 100)
        )
    )

    # Define the antibody kinetics model
    abkineticsModel <- list(
        model = makeModel(
                addAbkineticsModel("none", "sVNT", "none",  c("wane"), noInfSerumKinetics),
                addAbkineticsModelHier("inf", "sVNT", "inf", c("a"),  "a", as.numeric(covar_key), infSerumKinetics)
            ),
        prior = bind_rows(
            addPrior("wane", 0.0, 0.01, "unif", 0.0, 0.01), # observational model
            addPriorHier("a", 0.0, 6, "unif", 0.0, 6, "exp", 1, NA, 2) # ab kinetics
        )
    )

    model_cop <- createSeroJumpModel(
        data_sero = data_titre_model, 
        data_known = NULL,
        biomarkers = biomarkers,
        exposureTypes = exposureTypes,
        exposureFitted = exposureFitted,
        observationalModel = observationalModel,
        abkineticsModel = abkineticsModel,
        seed = seed_i)

    rj_settings <- list(
        numberChainRuns = 4, 
        iterations = 2000,
        consoleUpdates = 1000,
        burninPosterior = 1000,
        thin = 1
        
    )

    save_info_cop <- list(
        file_name = "test_data",
        model_name = "testA"
    )

    # Run these but take a while
    output_1 <- runSeroJump(model_cop, rj_settings, save_info = save_info_cop, seed = seed_i)
    set.seed(seed_i)
    output_2 <- runSeroJump(model_cop, rj_settings, save_info = save_info_cop, seed = seed_i)
   # output_3 <- readRDS(test_path("testdata", "dummy_output_123.rds"))

  #  output_1$fit$model$samplePriorDistribution
  #  output_3$fit$model$samplePriorDistribution


    expect_equal(output_1$post, output_2$post)
  #  expect_equal(output_1$post, output_3$post)

})