devtools::load_all()

library(lubridate)
library(patchwork)
library(tidybayes)
library(ggdist)

# Load data clean in the formet required for the model
data_exp <- get_data_titre_nih_2023()
data_titre <- data_exp[[1]] #%>% check_titre
known_exposure <- data_exp[[2]] %>% unique # must be unique
exp_prior <- data_exp[[3]]


obsLogLikelihood <- function(titre_val, titre_est, pars) {
    if (titre_val < 0) {
        ll <- pnorm(0, titre_est, pars[1], log.p = TRUE)
    } else if (titre_val >= 8){
        ll <- 1 - pnorm(8, titre_est, pars[1], log.p = TRUE)
    } else {
        ll <- log(pnorm(titre_val + 1, titre_est, pars[1]) - pnorm(titre_val , titre_est, pars[1]))
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
    titre_est <- max(titre_est, 0)
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
biomarkers <- c("A/Darwin/06/2021", "A/Darwin/09/2021e")
exposureTypes <- c("vax", "h3_2023")
exposureFitted <- c("h3_2023")

# Define the observational model
observationalModel <- list(
    model = makeModel(
        addObservationalModel("A/Darwin/06/2021", c("sigma"), obsLogLikelihood),
        addObservationalModel("A/Darwin/09/2021e", c("sigma_e"), obsLogLikelihood)
        ),
    prior = bind_rows(
        add_par_df("sigma", 0.0001, 4, "unif", 0.0001, 4),
        add_par_df("sigma_e", 0.0001, 4, "unif", 0.0001, 4),
        ) # observational model,
)

# Define the antibody kinetics model
abkineticsModel <- list(
    model = makeModel(
            addAbkineticsModel("vax", "A/Darwin/06/2021", "vax", c("a_vax", "b_vax", "c_vax"), infSerumKinetics),
            addAbkineticsModel("h3_2023","A/Darwin/06/2021", "h3_2023", c("a_v", "b_v", "c_v"), infSerumKinetics),
            addAbkineticsModel("vax_e", "A/Darwin/09/2021e", "vax", c("a_vax_e", "b_vax_e", "c_vax_e"), infSerumKinetics),
            addAbkineticsModel("h3_2023_e", "A/Darwin/09/2021e", "h3_2023", c("a_v_e", "b_v_e", "c_v_e"), infSerumKinetics)
        ),
    prior = bind_rows(
        add_par_df("a_vax", -2, 2, "norm",  0, 1), # ab kinetics
        add_par_df("b_vax", 0, 1, "norm",  0.3, 0.05), # ab kinetics
        add_par_df("c_vax", 0, 4, "unif", 0,  4), # ab kinetics 
        add_par_df("a_v", -2, 2, "norm",  0, 1), # ab kinetics
        add_par_df("b_v", 0, 1, "norm",  0.3, 0.05), # ab kinetics
        add_par_df("c_v", 0, 4, "unif", 0,  4), # ab kinetics
        add_par_df("a_vax_e", -2, 2, "norm",  0, 1), # ab kinetics
        add_par_df("b_vax_e", 0, 1, "norm",  0.3, 0.05), # ab kinetics
        add_par_df("c_vax_e", 0, 4, "unif", 0,  4), # ab kinetics 
        add_par_df("a_v_e", -2, 2, "norm",  0, 1), # ab kinetics
        add_par_df("b_v_e", 0, 1, "norm",  0.3, 0.05), # ab kinetics
        add_par_df("c_v_e", 0, 4, "unif", 0,  4) # ab kinetics
    )
)

# Define the COP model
copModel <- list( 
        model = makeModel(
            addCopModel("A/Darwin/06/2021", "h3_2023", c("beta0", "beta1"), copFuncForm, copLogLikelihood),
            addCopModel("A/Darwin/09/2021e", "h3_2023", c("beta0_e", "beta1_e"), copFuncForm, copLogLikelihood)
        ),
        prior = bind_rows(
            add_par_df("beta0", -10, 10, "unif", -10, 10), # cop model (not used here)
            add_par_df("beta1", -10, 10, "unif", -10, 10),
            add_par_df("beta0_e", -10, 10, "unif", -10, 10), # cop model (not used here)
            add_par_df("beta1_e", -10, 10, "unif", -10, 10)
        ) # cop model (not used here),
)

modeldefinition <- list(
    biomarkers = biomarkers,
    exposureTypes = exposureTypes,
    exposureFitted = exposureFitted,
    observationalModel = observationalModel,
    abkineticsModel = abkineticsModel,
    copModel = copModel,
    exposurePrior = exp_prior,
    exposurePriorType = "empirical"
)

seroModel <- createSeroJumpModel(data_titre, known_exposure, modeldefinition)

saveRDS(seroModel, file = here::here("hpc", "nih_2024", "nih_2024_model.RData"))
