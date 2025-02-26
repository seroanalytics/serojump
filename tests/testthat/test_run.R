test_that("runSeroJump produces reproducible output", {

    seed_i <- 123
    set.seed(seed_i)

    data_titre_model <- data.frame(
        id = rep(1:10, each = 10),
        time = rep(c(2, 5, 7, 12, 14, 15, 19, 30, 34, 40), 10),
        sVNT = rnorm(100, 100, 10)
    )

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
                addPrior("sigma", 0.0001, 4, "unif", 0.0001, 4)
        )
    )

    # Define the antibody kinetics model
    abkineticsModel <- list(
        model = makeModel(
                addAbkineticsModel("none", "sVNT", "none",  c("wane"), noInfSerumKinetics),
                addAbkineticsModel("inf", "sVNT", "inf", c("a"), infSerumKinetics)
            ),
        prior = bind_rows(
            addPrior("wane", 0.0, 0.01, "unif", 0.0, 0.01), # observational model
            addPrior("a", 0, 6, "unif",  0, 6), # ab kinetics
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
    output_3 <- readRDS(test_path("testdata", "dummy_output_123.rds"))

  #  output_1$fit$model$samplePriorDistribution
  #  output_3$fit$model$samplePriorDistribution


    expect_equal(output_1$post, output_2$post)
    expect_equal(output_1$post, output_3$post)

})