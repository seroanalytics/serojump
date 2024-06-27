library(devtools)

library(Rcpp)
devtools::load_all()
i <- Sys.getenv("SLURM_ARRAY_TASK_ID")
i <- as.integer(i)

seroModel_full <- readRDS(here::here("hpc", "nih_2024", "nih_2024_model.RData"))

output_file <- c("h1", "h3")[i]


seroModel <- seroModel_full[[i]]

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

ab_values <- c(0:10)
names_ab <- 2^c(0:10) * 5
names_ab[2] <- "<10"
names(ab_values) <- names_ab

runRJMCMC(seroModel, settings, "hpc/nih_2024", output_file)
postprocessFigs("hpc/nih_2024", output_file, 4, ab_values)
