library(devtools)
library(Rcpp)

devtools::load_all()
i <- 2
seroW2_full <- readRDS(here::here("hpc", "transvir_w2_inf", "transvir_w2_model.RData"))
prior_names <- c("p1", "p2", "p3")
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
    noGibbsSteps = 1,
    onDebug = FALSE
)

runInfRJMCMC(seroW2_full[[i]], settings, paste0("hpc/transvir_w2_inf/", prior_names[i]), "w2")
postprocessFigsInf(paste0("hpc/transvir_w2_inf/",  prior_names[i]), "w2", 4)
