library(devtools)
library(Rcpp)
library(data.table)

devtools::load_all()

i <- 6
seroModel_full <- readRDS(here::here("hpc", "nih_2024_inf", "nih_2024_model.RData"))
name_vec <- names(seroModel_full)
name_vec_spl <- str_split(name_vec, "_")

seroModel <- seroModel_full[[i]]
vec_names_i <- name_vec_spl[[i]]

settings <-  list(
    numberChainRuns = 4,
    numberCores = 4,
    iterations = 500,
    burninPosterior = 250,
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


runInfRJMCMC(seroModel, settings, paste0("hpc/nih_2024_inf/", vec_names_i[2]), vec_names_i[1])
postprocessFigsInf(paste0("hpc/nih_2024_inf/", vec_names_i[2]), vec_names_i[1], 4, ab_values)