devtools::load_all()

library(lubridate)
library(patchwork)
library(tidybayes)
library(ggdist)

# Load data clean in the formet required for the model
gambia_pvnt_w2 <- get_data_titre_model_wave2() %>% rename(sVNT = titre)# This is the empirical prior for the exposure time
gambia_exp_w2 <- get_exposures_wave2() # This is the empirical prior for the exposure time
exp_prior_w2 <- get_exp_prior_wave2() # This is the empirical prior for the exposure time


obsLogLikelihood = function(titre_val, titre_est, pars) {
    if (titre_val <= log10(40)) {
        ll <- pcauchy(log10(40), titre_est, pars[1], log.p = TRUE)
    } else {
        ll <- dcauchy(titre_val, titre_est, pars[1], log = TRUE)
    }
    ll
}

copFuncForm <- function(inf_status, esttitreExp, params) {
    beta0 <- params[1]
    beta1 <- params[2]
    p <- 1.0 / (1.0 + exp(- (beta0 + beta1 * esttitreExp) ) )
}

copLogLikelihood <- function(inf_status, esttitreExp, params) {
    # COP parameters
    beta0 <- params[1]
    beta1 <- params[2]
    p <- 1.0 / (1.0 + exp(- (beta0 + beta1 * esttitreExp) ) )
    ll <- inf_status * log(p) + (1 - inf_status) * log(1 - p)
    ll
}

noInfSerumKinetics <- function(titre_est, timeSince, pars) {
    titre_est <- titre_est - pars[1] * (timeSince)
    titre_est <- max(0, titre_est)
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
    titre_est <- max(0, titre_est)
}



# Define the biomarkers and exposure types in the model
biomarkers <- c("sVNT")
exposureTypes <- c("none", "delta", "vax", "predelta")
exposureFitted <- "delta"

# Define the observational model
observationalModel <- list(
    names = c("sVNT"),
    model = makeModel(addObservationalModel("sVNT", c("sigma"), obsLogLikelihood)),
    prior = add_par_df("sigma", 0.0001, 4, "unif", 0.0001, 4) # observational model,
)

# Define the antibody kinetics model
abkineticsModel <- list(
    model = makeModel(
            addAbkineticsModel("none", "sVNT", "none",  c("wane"), noInfSerumKinetics),
            addAbkineticsModel("delta", "sVNT", "delta", c("a_d", "b_d", "c_d"), infSerumKinetics),
            addAbkineticsModel("vax", "sVNT", "vax", c("a_vax", "b_vax", "c_vax"), infSerumKinetics),
            addAbkineticsModel("predelta", "sVNT", "predelta", c("a_pd", "b_pd", "c_pd"), infSerumKinetics)
        ),
    prior = bind_rows(
        add_par_df("a_vax", -2, 2, "norm",  0, 1), # ab kinetics
        add_par_df("b_vax", 0, 1, "norm",  0.3, 0.05), # ab kinetics
        add_par_df("c_vax", 0, 4, "unif", 0,  4), # ab kinetics 
        add_par_df("a_pd", -2, 2, "norm",  0, 1), # ab kinetics
        add_par_df("b_pd", 0, 1, "norm",  0.3, 0.05), # ab kinetics
        add_par_df("c_pd", 0, 4, "unif", 0,  4), # ab kinetics
        add_par_df("a_d", -2, 2, "norm",  0, 1), # ab kinetics
        add_par_df("b_d", 0, 1, "norm",  0.3, 0.05), # ab kinetics
        add_par_df("c_d", 0, 4, "unif", 0,  4), # ab kinetics
        add_par_df("wane", 0.0, 0.1, "unif", 0.0, 0.1) # observational model
    )
)

# Define the COP model
copModel <- list( 
        names = c("sVNT"),
        model = makeModel(addCopModel("sVNT", "delta", c("beta0", "beta1"), copFuncForm, copLogLikelihood)),
        prior = bind_rows(
            add_par_df("beta0", -10, 10, "unif", -10, 10), # cop model (not used here)
            add_par_df("beta1", -10, 10, "unif", -10, 10)
        ) # cop model (not used here),
)

modeldefinition <- list(
    biomarkers = biomarkers,
    exposureTypes = exposureTypes,
    exposureFitted = exposureFitted,
    observationalModel = observationalModel,
    abkineticsModel = abkineticsModel,
    copModel = copModel,
    exposurePrior = exp_prior_w2,
    exposurePriorType = "empirical"
)

modelW2 <- createSeroJumpModel(gambia_pvnt_w2, gambia_exp_w2, modeldefinition)
saveRDS(modelW2, file = here::here("hpc", "transvir_w2", "transvir_w2_model.RData"))
