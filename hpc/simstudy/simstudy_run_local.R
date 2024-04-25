library(devtools)

library(Rcpp)
#Rcpp::compileAttributes()
devtools::load_all()
i <- 7

seroModel_full <- readRDS(here::here("hpc", "simstudy", "simstudy_model.RData"))

## check entriee
#seroModel$data$times_list # id, t
#seroModel$data$titre_list #id, bio, t

settings <-  list(
    numberChainRuns = 4,
    numberCores = 4,
    iterations = 10000,
    burninPosterior = 5000,
    thin = 10,
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

runRJMCMC(seroModel_full[[i]]$model, settings, seroModel_full[[i]]$filename, seroModel_full[[i]]$modelname)
postprocessFigs(seroModel_full[[i]]$filename, seroModel_full[[i]]$modelname, 4)
