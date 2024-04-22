library(devtools)

library(Rcpp)
#Rcpp::compileAttributes()
devtools::load_all()

seroModel <- readRDS(here::here("hpc", "transvir_w2", "transvir_w2_model.RData"))

## check entriee
#seroModel$data$times_list # id, t
#seroModel$data$titre_list #id, bio, t

settings <-  list(
    numberChainRuns = 4,
    numberCores = 4,
    iterations = 400000,
    burninPosterior = 200000,
    thin = 1000,
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

runRJMCMC(seroModel, settings, "hpc/transvir", "wave2")
postprocessFigs("hpc/transvir", "wave2", 4)
