library(devtools)
library(Rcpp)

devtools::load_all()
i <- 6
seroW2_full <- readRDS(here::here("hpc", "transvir_w2_inf", "transvir_w23_model.RData"))
prior_names <- c("p1", "p2", "p3", "p1", "p2", "p3")
prior_names_wave <- c("w2", "w2", "w2", "w3", "w3", "w3")

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

runInfRJMCMC(seroW2_full[[i]], settings, paste0("hpc/transvir_w2_inf/", prior_names[i]), prior_names_wave[i])
postprocessFigsInf(paste0("hpc/transvir_w2_inf/",  prior_names[i]), prior_names_wave[i], 4)
