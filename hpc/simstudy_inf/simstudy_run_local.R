library(devtools)

library(Rcpp)
#Rcpp::compileAttributes()
devtools::load_all()
i <- Sys.getenv("SLURM_ARRAY_TASK_ID")
i <- as.integer(i)

seroModel_full <- readRDS(here::here("hpc", "simstudy_inf", "simstudy_model.RData"))

## check entriee
#seroModel$data$times_list # id, t
#seroModel$data$titre_list #id, bio, t

settings <-  list(
    numberChainRuns = 4,
    numberCores = 4,
    iterations = 1000,
    burninPosterior = 500,
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

runInfRJMCMC(seroModel_full[[i]]$model, settings, seroModel_full[[i]]$filename, seroModel_full[[i]]$modelname)
postprocessFigsInf(seroModel_full[[i]]$filename, seroModel_full[[i]]$modelname, 4)
