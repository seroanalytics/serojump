library(devtools)

library(Rcpp)
devtools::load_all()

c("h1", "h3")
seroModel_full <- readRDS(here::here("hpc", "nih_2024", "nih_2024_model.RData"))
seroModel <- seroModel_full[[1]]

settings <-  list(
    numberChainRuns = 4,
    numberCores = 4,
    iterations = 5000,
    burninPosterior = 2500,
    thin = 10,
    consoleUpdates = 100,
    onAdaptiveCov = TRUE,
    updatesAdaptiveCov = 10,
    burninAdaptiveCov = 1000,
    covarInitVal = 1e-2, # make very small if struggling to sample to beginning
    covarInitValAdapt = 1e-2, # make very small if struggling to sample to beginning
    covarMaxVal = 1, # decrease if struggling toc sample in the middle
    runParallel = TRUE,
    noGibbsSteps = 1,
    onDebug = FALSE
)


ab_values <- c(0:10)
names_ab <- 2^c(0:10) * 5
names_ab[2] <- "<10"
names(ab_values) <- names_ab

runInfRJMCMC(seroModel, settings, "hpc/nih_2024_inf/prior1", "h1")
postprocessFigs("hpc/nih_2024/prior1", "h1", 4, ab_values)

