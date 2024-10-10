library(devtools)
library(serojump)
library(tidyverse)


obsLogLikelihood <- function(titre_val, titre_est, pars) {
    ll <- dnorm(titre_val, titre_est, pars[1], log = TRUE)
}

noInfSerumKinetics <- function(titre_est, timeSince, pars) {
    titre_est_new <- titre_est - (pars[1] * titre_est) * (timeSince)
    titre_est_new <- max(titre_est_new, 0)
    titre_est_new
}

infTuenisPower2016 <- function(titre_est, timeSince, pars) {

    y1 <- pars[1]
    t1 <- pars[2]
    r <- pars[3]
    alpha <- pars[4]

    v <- 0.001
    mu <- 1 / t1 * y1

    if (timeSince < t1) {
        titre_est_boost <- exp(mu * timeSince)
    } else {
        titre_est_boost <- exp(y1) * (1 + (r - 1) * exp(y1)^{r - 1} * v * (timeSince - t1)) ^ {-1 / (r - 1)}
    }

    titre_est_log <- titre_est + log(titre_est_boost) * max(0, 1 - titre_est * alpha)
    titre_est_log 
}



data_titre_e1 <- read.csv( file = here::here("data", "fudan", "df_soro_model_e1.csv")) %>% select(!X) 
exp_prior_e1 <- read.csv(file = here::here("data", "fudan", "exp_prior_e1.csv")) %>% select(!X)

# Define the biomarkers and exposure types in the model
biomarkers_e1 <- c("PreF")
exposureTypes_e1 <- c("none", "inf")
exposureFitted_e1 <- c("inf")

# Define the observational model
observationalModel_e1 <- list(
    names = c("PreF"),
    model = makeModel(
        addObservationalModel("PreF", c("sigma"), obsLogLikelihood)
        ),
    prior = bind_rows(
        add_par_df("sigma", 0.0001, 4, "unif", 0.0001, 4)
        ) # observational model,
)


# Define the antibody kinetics model
abkineticsModel_e1 <- list(
    model = makeModel(
            addAbkineticsModel("none", "PreF", "none", c("wane"), noInfSerumKinetics),
            addAbkineticsModel("inf", "PreF", "inf", c("y1_h1", "t1_h1", "r_h1", "s"), infTuenisPower2016)
        ),
    prior = bind_rows(
        add_par_df("y1_h1", 1, 6, "unif",  1, 6), # ab kinetics
        add_par_df("t1_h1", 7, 24, "norm", 14, 3), # ab kinetics
        add_par_df("r_h1", 1, 5, "unif", 1, 5), # ab kinetics 
        add_par_df("s", 0, 1, "unif", 0, 1), # ab kinetics 
        add_par_df("wane", 0, 1, "unif", 0, 1) # ab kinetics
    )
)

inf_prior_1 <- function(N, E, I, K) {
    N_adj <- N - K
    E_adj <- E - K 
    logPriorExpInf <- lfactorial(E_adj) + lfactorial(N_adj - E_adj) - lfactorial(N_adj ) 
    logPriorExpInf
}

inf_prior_2 <- function(N, E, I, K) {
    N_adj <- N - K
    E_adj <- E - K
    logPriorExpInf <- lfactorial(E_adj) + lfactorial(N_adj - E_adj) - lfactorial(N_adj ) + dbinom(E_adj, N_adj, 0.05, log = TRUE)
    logPriorExpInf
}


modeldefinition_e1_p1 <- list(
    biomarkers = biomarkers_e1,
    exposureTypes = exposureTypes_e1,
    exposureFitted = exposureFitted_e1,
    observationalModel = observationalModel_e1,
    abkineticsModel = abkineticsModel_e1,
    expInfPrior = inf_prior_1,
    exposurePrior = exp_prior_e1,
    exposurePriorType = "empirical"
)

modeldefinition_e1_p2 <- modeldefinition_e1_p1
modeldefinition_e1_p2$expInfPrior <- inf_prior_2

seroModel_phirst_e1_p1 <- createSeroJumpModel(data_titre_e1, NULL, modeldefinition_e1_p1)
seroModel_phirst_e1_p2 <- createSeroJumpModel(data_titre_e1, NULL, modeldefinition_e1_p2)

outputs_modelA <- list(seroModel_phirst_e1_p1, seroModel_phirst_e1_p2)


save_info_list <- list(
    list(file_name = "fudan_e1", model_name = "e1_base"),
    list(file_name = "fudan_e1", model_name = "e1_inform")
)

model_run_info <- 
    list(model = outputs_modelA, save = save_info_list)
saveRDS(model_run_info, file = here::here("data", "fudan", "model_A.RDS"))
