devtools::load_all()

library(lubridate)
library(patchwork)
library(tidybayes)
library(ggdist)

# Load data clean in the formet required for the model
data_known <- NULL # This is the empirical prior for the exposure time
exp_prior <- data.frame(
    lb = "1", # lower bound
    ub = "120", # upper bound
    func = "unif", # functional form
    par1 = "1", # parameter 1 of functional form
    par2 = "120" # parameter 2 of functional form
) # This is the empirical prior for the exposure time


## Define the functions required for the rjmcmc model

obsLogLikelihood = function(titre_val, titre_est, pars) {
    ll <- dnorm(titre_val, titre_est, pars[1], log = TRUE)
}

noInfSerumKinetics <- function(titre_est, timeSince, pars) {
    titre_est
}

infSerumKinetics <- function(titre_est, timeSince, pars) {
    a <- pars[1]
    b <- pars[2] 
    c <- pars[3]
    if (timeSince < 14) {
        titre_est <- titre_est + log(exp(a) + exp(c)) * (timeSince) / 14;
    } else {
        titre_est <- titre_est + log(exp(a) * exp(-b/10 * (timeSince - 14)) + exp(c));
    }
    titre_est
}




## Define the models for the rjmcmc model

# Define the biomarkers and exposure types in the model
biomakers <- c("IgG")
exposureTypes <- c("none", "delta")
exposureFitted <- c("delta")

# Define the observational model
observationalModel <- list(
    names = c("IgG"),
    model = makeModel(addObservationalModel("IgG", c("sigma"), obsLogLikelihood)),
    prior = add_par_df("sigma", 0.0001, 4, "unif", 0.0001, 4) # observational model,
)

# Define the antibody kinetics model
abkineticsModel <- list(
    model = makeModel(
            addAbkineticsModel("none", "IgG", "none", c("a"), noInfSerumKinetics),
            addAbkineticsModel("delta", "IgG", "delta", c("a", "b", "c"), infSerumKinetics)
        ),
    prior = bind_rows(
        add_par_df("a", -2, 2, "norm",  0, 1), # ab kinetics
        add_par_df("b", 0, 1, "norm",  0.3, 0.05), # ab kinetics
        add_par_df("c", 0, 4, "unif", 0,  4), # ab kinetics
    )
)

inf_prior_1 <- function(N, E, I, K) {
    0
}


inf_prior_2 <- function(N, E, I, K) {
    N_adj <- N - K
    E_adj <- E - K 
    logPriorExpInf <- lfactorial(E_adj) + lfactorial(N_adj - E_adj) - lfactorial(N_adj )
    logPriorExpInf
}

inf_prior_3 <- function(N, E, I, K) {
    N_adj <- N - K
    E_adj <- E - K 
    logPriorExpInf <- lfactorial(E_adj) + lfactorial(N_adj - E_adj) - lfactorial(N_adj ) + dbinom(E_adj, N_adj, 0.3, log = TRUE)
    logPriorExpInf
}

modeldefinition_p1 <- list(
    biomarkers = biomakers,
    exposureTypes = exposureTypes,
    exposureFitted = exposureFitted,
    observationalModel = observationalModel,
    abkineticsModel = abkineticsModel,
    expInfPrior = inf_prior_1,
    exposurePrior = exp_prior,
    exposurePriorType = "func",
    prior_name = "p1"
)

modeldefinition_p2 <- modeldefinition_p1
modeldefinition_p2$expInfPrior <- inf_prior_2
modeldefinition_p2$prior_name <- "p2"

modeldefinition_p3 <- modeldefinition_p1
modeldefinition_p3$expInfPrior <- inf_prior_3
modeldefinition_p3$prior_name <- "p3"


# Define the settings
settings <-  list(
    numberChainRuns = 4,
    numberCores = 4,
    iterations = 5000,
    burninPosterior = 2500,
    thin = 10,
    consoleUpdates = 100,
    onAdaptiveCov = TRUE,
    updatesAdaptiveCov = 100,
    burninAdaptiveCov = 1000,
    updatesAdaptiveTemp = 10,
    covarInitVal = 1e-2, # make very small if struggling to sample to beginning
    covarInitValAdapt = 1e-2, # make very small if struggling to sample to beginning
    covarMaxVal = 1, # decrease if struggling to sample in the middle
    runParallel = TRUE,
    onDebug = FALSE,
    noGibbsSteps = 1
)

modeldefinition_p1 <- modeldefinition_p1
modeldefinition_p2$expInfPrior <- inf_prior_2

modeldefinition_p3 <- modeldefinition_p1
modeldefinition_p3$expInfPrior <- inf_prior_3


sim_data_vec <- c("cesCOP_notd", "cesNoCOP_notd")
obs_vec <- c("0.5")
priors <- c("p1", "p2", "p3")
model_list <- vector(mode = "list", length = 6)
# unknown exposure
for (sim_data in sim_data_vec) {
    for (obs_er in obs_vec) {
        for (modeldefinition in list(modeldefinition_p1, modeldefinition_p2, modeldefinition_p3)) {
            data_titre_model <- clean_simulated_rjmcmc(sim_data, obs_er) %>% rename(IgG = titre)
            modelname <- paste0("obs_", obs_er)

            model_list[[i]] <- list(
                filename = paste0("hpc/simstudy_inf/", sim_data, "/", modeldefinition$prior_name),
                modelname = modelname,
                model = createSeroJumpModel(data_titre_model, data_known, modeldefinition)
            )
            i <- i + 1
        }
    }
}

saveRDS(model_list, file = here::here("hpc", "simstudy_inf", "simstudy_model.RData"))

#"local/sim_study/cesCOP_notd/InferExp"