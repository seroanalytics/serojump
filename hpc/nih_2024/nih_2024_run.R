library(devtools)

library(Rcpp)
devtools::load_all()

seroModel <- readRDS(here::here("hpc", "nih_2024", "nih_2024_model.RData"))

settings <-  list(
    numberChainRuns = 4,
    numberCores = 4,
    iterations = 100000,
    burninPosterior =250000,
    thin = 100,
    consoleUpdates = 100,
    onAdaptiveCov = TRUE,
    updatesAdaptiveCov = 10,
    burninAdaptiveCov = 1000,
    covarInitVal = 1e-2, # make very small if struggling to sample to beginning
    covarInitValAdapt = 1e-2, # make very small if struggling to sample to beginning
    covarMaxVal = 1, # decrease if struggling toc sample in the middle
    runParallel = TRUE,
    noGibbsSteps = 5,
    onDebug = FALSE
)


runRJMCMC(seroModel, settings, "hpc/nih_2024", "h3")
postprocessFigs("hpc/nih_2024", "h3", 4)
