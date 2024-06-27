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




obsLogLikelihood = function(titre_val, titre_est, pars) {
    ll <- dnorm(titre_val, titre_est, pars[1], log = TRUE)
}

copFuncForm <- function(inf_status, esttitreExp, params) {
    beta0 <- params[1]
    beta1 <- params[2]
    p <- 1.0 / (1.0 + exp(- (beta0 + (-beta0 / 4 - b1_sample) * esttitreExp) ) )
}

#b0_sample <- runif(1000, 0, 10)
#a_sample <- runif(1000, 0, 5)
#b1_sample <- -b0_sample/4 - a_sample


copLogLikelihood <- function(inf_status, esttitreExp, pars) {
    # COP parameters
    beta0 <- pars[1]
    beta1 <- pars[2]
    p <- 1.0 / (1.0 + exp(- (beta0 + (-beta0 / 4 - b1_sample) * esttitreExp) ) )
    ll <- inf_status * log(p) + (1 - inf_status) * log(1 - p)
    ll <- 0
    ll
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

# Define the COP model
copModel <- list( 
        names = c("IgG"),
        model = makeModel(addCopModel("IgG", "delta", c("beta0", "beta1"), copFuncForm, copLogLikelihood)),
        prior = bind_rows(
            add_par_df("beta0", -10, 10, "unif", 0, 10), # cop model (not used here)
            add_par_df("beta1", -10, 10, "unif", 0, 5)
        ) # cop model (not used here),
)


modeldefinition <- list(
    biomarkers = biomakers,
    exposureTypes = exposureTypes,
    exposureFitted = exposureFitted,
    observationalModel = observationalModel,
    abkineticsModel = abkineticsModel,
    copModel = copModel,
    exposurePrior = exp_prior,
    exposurePriorType = "func"
)

sim_data_vec <- c("cesCOP_notd", "cesNoCOP_notd")
obs_vec <- c("0.1", "0.3", "0.5")
model_list <- vector(mode = "list", length = 12)
i <- 1

#known exposure
for (sim_data in sim_data_vec) {
    for (obs_er in obs_vec) {
        data_titre_model <- clean_simulated_rjmcmc(sim_data, obs_er, 0) %>% rename(IgG = titre)
        modelname <- paste0("obs_", obs_er)

        modeli <- readRDS(here::here("outputs", "sim_data", sim_data, "inputs.RDS"))
         known_exp <- modeli$simpar$exp %>% as.data.frame %>% 
            mutate(i = 1:modeli$simpar$N) %>% pivot_longer(!i, names_to = "t", values_to = "exp") %>% mutate(t = as.numeric(substr(t, 2, 4))) %>% 
            filter(exp == 1) %>% complete(i = 1:modeli$simpar$N, fill = list(t = -1, exp = -1)) %>% pull(t)

        model_list[[i]] <- list(
            filename = paste0("hpc/simstudy_noCOP/", sim_data, "/KnownExp"),
            modelname = modelname,
            model = createSeroJumpModel(data_titre_model, data_known, modeldefinition, known_exp)
        )
        i <- i + 1
    }
}

# unknown exposure
# # 7–12
for (sim_data in sim_data_vec) {
    for (obs_er in obs_vec) {
        data_titre_model <- clean_simulated_rjmcmc(sim_data, obs_er, 0, known_exp = FALSE) %>% rename(IgG = titre)
        modelname <- paste0("obs_", obs_er)



        model_list[[i]] <- list(
            filename = paste0("hpc/simstudy_noCOP/", sim_data, "/InferExp"),
            modelname = modelname,
            model = createSeroJumpModel(data_titre_model, data_known, modeldefinition)
        )
        i <- i + 1
    }
}

saveRDS(model_list, file = here::here("hpc", "simstudy_noCOP", "simstudy_model.RData"))

#"local/sim_study/cesCOP_notd/InferExp"