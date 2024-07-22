devtools::load_all()

library(lubridate)
library(patchwork)
library(tidybayes)
library(ggdist)

obsLogLikelihood <- function(titre_val, titre_est, pars) {
    if (titre_val < 0) {
        ll <- pnorm(0, titre_est, pars[1], log.p = TRUE)
    } else if (titre_val >= 12){
        ll <- 1 - pnorm(12, titre_est, pars[1], log.p = TRUE)
    } else {
        ll <- log(pnorm(titre_val + 1, titre_est, pars[1]) - pnorm(titre_val , titre_est, pars[1]))
    }
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


infSerumKinetics_titredep <- function(titre_est, timeSince, pars) {
    a <- pars[1]
    b <- pars[2]
    c <- pars[3]
    s <- pars[4]

    if (timeSince < 14) {
        titre_est <- titre_est + max(1 - s * titre_est, 0) * log(exp(a) + exp(c)) * (timeSince) / 14;
    } else {
        titre_est <- titre_est + max(1 - s * titre_est, 0) * (log(exp(a) * exp(-b/10 * (timeSince - 14)) + exp(c)));
    }
    titre_est
}

########################################################################
######## H1 ############
########################################################################

data_exp <- get_data_titre_nih_2023_h1()
data_titre_h1 <- data_exp[[1]] #%>% check_titre
known_exposure_h1 <- data_exp[[2]] %>% unique # must be unique
exp_prior <- data_exp[[3]]


# Get the number of four fold rises between point 2 and 3

known_exposure_h1_add <- known_exposure_h1 %>% filter(exposure_type == "h1_2023") %>% select(pid, id) %>% mutate(type = "known")
sero_confirmed <- left_join(
    data_titre_h1 %>% group_by(id) %>% mutate(r = row_number(), r_max = max(r)) %>% filter(r_max >=3) %>% 
        filter(r == 2 ) %>% rename(titre_2 = `A/Sydney/5/2021`, titre_2e = `A/Sydney/5/2021e` ) %>% 
        select(pid, id, titre_2, titre_2e),

    data_titre_h1 %>% group_by(id) %>% mutate(r = row_number(), r_max = max(r)) %>% filter(r_max >=3) %>% 
        filter(r == 3 ) %>% rename(titre_3 = `A/Sydney/5/2021`, titre_3e = `A/Sydney/5/2021e` ) %>% 
        select(pid, id, titre_3, titre_3e)
    ) %>% mutate(titre_diff = titre_3 - titre_2, titre_diffe = titre_3e - titre_2e ) %>%
    mutate(diff = if (titre_diff >= 2 & titre_diffe >= 2) {four_fold = 1} else {four_fold = 0}) %>% filter(diff == 1)

sero_confirmed_trim <- sero_confirmed %>% left_join(known_exposure_h1_add) %>% filter(is.na(type))
rate_inf_h1 <- (sero_confirmed_trim %>% nrow) / (data_titre_h1$id %>% unique %>% length)

# Define the biomarkers and exposure types in the model
biomarkers_h1 <- c("A/Sydney/5/2021", "A/Sydney/5/2021e")
exposureTypes_h1 <- c("vax", "h1_2023")
exposureFitted_h1 <- c("h1_2023")

# Define the observational model
observationalModel_h1 <- list(
    model = makeModel(
        addObservationalModel("A/Sydney/5/2021", c("sigma"), obsLogLikelihood),
        addObservationalModel("A/Sydney/5/2021e", c("sigma_e"), obsLogLikelihood)
        ),
    prior = bind_rows(
        add_par_df("sigma", 0.0001, 4, "unif", 0.0001, 4),
        add_par_df("sigma_e", 0.0001, 4, "unif", 0.0001, 4),
        ) # observational model,
)

# Define the antibody kinetics model
abkineticsModel_h1 <- list(
    model = makeModel(
            addAbkineticsModel("vax", "A/Sydney/5/2021", "vax", c("a_vax", "b_vax", "c_vax", "s"), infSerumKinetics_titredep),
            addAbkineticsModel("h1_2023","A/Sydney/5/2021", "h1_2023", c("a_v", "b_v", "c_v", "s"), infSerumKinetics_titredep),
            addAbkineticsModel("vax_e", "A/Sydney/5/2021e", "vax", c("a_vax_e", "b_vax_e", "c_vax_e", "s_e"), infSerumKinetics_titredep),
            addAbkineticsModel("h1_2023_e", "A/Sydney/5/2021e", "h1_2023", c("a_v_e", "b_v_e", "c_v_e", "s_e"), infSerumKinetics_titredep)
        ),
    prior = bind_rows(
        add_par_df("a_vax", -5, 6, "norm",  1, 2), # ab kinetics
        add_par_df("b_vax", 0, 1, "norm",  0.3, 0.05), # ab kinetics
        add_par_df("c_vax", 0, 4, "unif", 0, 4), # ab kinetics 
        add_par_df("a_v", -5, 6, "norm",  1, 2), # ab kinetics
        add_par_df("b_v", 0, 1, "norm",  0.3, 0.05), # ab kinetics
        add_par_df("c_v", 1, 4, "unif", 1, 4), # ab kinetics
        add_par_df("s", 0, 1, "unif", 0, 1), # ab kinetics 
        add_par_df("a_vax_e", -5, 6, "norm",  1, 2), # ab kinetics
        add_par_df("b_vax_e", 0, 1, "norm",  0.3, 0.05), # ab kinetics
        add_par_df("c_vax_e", 0, 4, "unif", 0, 4), # ab kinetics 
        add_par_df("a_v_e", -5, 6, "norm",  1, 2), # ab kinetics
        add_par_df("b_v_e", 0, 1, "norm",  0.3, 0.05), # ab kinetics
        add_par_df("c_v_e", 1, 4, "unif", 1, 4), # ab kinetics
        add_par_df("s_e", 0, 1, "unif", 0, 1) # ab kinetics 
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
    logPriorExpInf <- lfactorial(E_adj) + lfactorial(N_adj - E_adj) - lfactorial(N_adj ) + dbinom(E_adj, N_adj, 0.008482564, log = TRUE)
    logPriorExpInf
}


# Define the COP model


modeldefinition_h1_p1 <- list(
    biomarkers = biomarkers_h1,
    exposureTypes = exposureTypes_h1,
    exposureFitted = exposureFitted_h1,
    observationalModel = observationalModel_h1,
    abkineticsModel = abkineticsModel_h1,
    expInfPrior = inf_prior_1,
    exposurePrior = exp_prior,
    exposurePriorType = "empirical"
)

modeldefinition_h1_p2 <- modeldefinition_h1_p1
modeldefinition_h1_p2$expInfPrior <- inf_prior_2

modeldefinition_h1_p3 <- modeldefinition_h1_p1
modeldefinition_h1_p3$expInfPrior <- inf_prior_3


seroModel_nih_h1_p1 <- createSeroJumpModel(data_titre_h1, known_exposure_h1, modeldefinition_h1_p1, TRUE)
seroModel_nih_h1_p2 <- createSeroJumpModel(data_titre_h1, known_exposure_h1, modeldefinition_h1_p2)
seroModel_nih_h1_p3 <- createSeroJumpModel(data_titre_h1, known_exposure_h1, modeldefinition_h1_p3)



########################################################################
######## H3 ############
########################################################################


data_exp_h3 <- get_data_titre_nih_2023_h3()
data_titre_h3 <- data_exp_h3[[1]] #%>% check_titre
known_exposure_h3 <- data_exp_h3[[2]] %>% unique # must be unique


known_exposure_h3_add <- known_exposure_h3 %>% filter(exposure_type == "h3_2023") %>% select(pid, id) %>% mutate(type = "known")
sero_confirmed <- left_join(
    data_titre_h3 %>% group_by(id) %>% mutate(r = row_number(), r_max = max(r)) %>% filter(r_max >=3) %>% 
        filter(r == 2 ) %>% rename(titre_2 = `A/Darwin/06/2021`, titre_2e = `A/Darwin/09/2021e` ) %>% 
        select(pid, id, titre_2, titre_2e),

    data_titre_h3 %>% group_by(id) %>% mutate(r = row_number(), r_max = max(r)) %>% filter(r_max >=3) %>% 
        filter(r == 3 ) %>% rename(titre_3 = `A/Darwin/06/2021`, titre_3e = `A/Darwin/09/2021e` ) %>% 
        select(pid, id, titre_3, titre_3e)
    ) %>% mutate(titre_diff = titre_3 - titre_2, titre_diffe = titre_3e - titre_2e ) %>%
    mutate(diff = if (titre_diff >= 2 & titre_diffe >= 2) {four_fold = 1} else {four_fold = 0}) %>% filter(diff == 1)

sero_confirmed_trim <- sero_confirmed %>% left_join(known_exposure_h3_add) %>% filter(is.na(type))
rate_inf_h3 <- (sero_confirmed_trim %>% nrow) / (data_titre_h3$id %>% unique %>% length)


# Define the biomarkers and exposure types in the model
biomarkers_h3 <- c("A/Darwin/06/2021", "A/Darwin/09/2021e")
exposureTypes_h3 <- c("vax", "h3_2023")
exposureFitted_h3 <- c("h3_2023")

# Define the observational model
observationalModel_h3 <- list(
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
abkineticsModel_h3 <- list(
    model = makeModel(
            addAbkineticsModel("vax", "A/Darwin/06/2021", "vax", c("a_vax", "b_vax", "c_vax", "s"), infSerumKinetics_titredep),
            addAbkineticsModel("h3_2023","A/Darwin/06/2021", "h3_2023", c("a_v", "b_v", "c_v", "s"), infSerumKinetics_titredep),
            addAbkineticsModel("vax_e", "A/Darwin/09/2021e", "vax", c("a_vax_e", "b_vax_e", "c_vax_e", "s_e"), infSerumKinetics_titredep),
            addAbkineticsModel("h3_2023_e", "A/Darwin/09/2021e", "h3_2023", c("a_v_e", "b_v_e", "c_v_e", "s_e"), infSerumKinetics_titredep)
        ),
 prior = bind_rows(
        add_par_df("a_vax", -5, 6, "norm",  1, 2), # ab kinetics
        add_par_df("b_vax", 0, 1, "norm",  0.3, 0.05), # ab kinetics
        add_par_df("c_vax", 0, 4, "unif", 0, 4), # ab kinetics 
        add_par_df("a_v", -5, 6, "norm",  1, 2), # ab kinetics
        add_par_df("b_v", 0, 1, "norm",  0.3, 0.05), # ab kinetics
        add_par_df("c_v", 1, 4, "unif", 1, 4), # ab kinetics
        add_par_df("s", 0, 1, "unif", 0, 1), # ab kinetics 
        add_par_df("a_vax_e", -5, 6, "norm",  1, 2), # ab kinetics
        add_par_df("b_vax_e", 0, 1, "norm",  0.3, 0.05), # ab kinetics
        add_par_df("c_vax_e", 0, 4, "unif", 0, 4), # ab kinetics 
        add_par_df("a_v_e", -5, 6, "norm",  1, 2), # ab kinetics
        add_par_df("b_v_e", 0, 1, "norm",  0.3, 0.05), # ab kinetics
        add_par_df("c_v_e", 1, 4, "unif", 1, 4), # ab kinetics
        add_par_df("s_e", 0, 1, "unif", 0, 1) # ab kinetics 
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
    logPriorExpInf <- lfactorial(E_adj) + lfactorial(N_adj - E_adj) - lfactorial(N_adj ) + dbinom(E_adj, N_adj,  0.005655042, log = TRUE)
    logPriorExpInf
}
modeldefinition_h3_p1 <- list(
    biomarkers = biomarkers_h3,
    exposureTypes = exposureTypes_h3,
    exposureFitted = exposureFitted_h3,
    observationalModel = observationalModel_h3,
    abkineticsModel = abkineticsModel_h3,
    expInfPrior = inf_prior_1,
    exposurePrior = exp_prior,
    exposurePriorType = "empirical"
)

modeldefinition_h3_p2 <- modeldefinition_h3_p1
modeldefinition_h3_p2$expInfPrior <- inf_prior_2

modeldefinition_h3_p3 <- modeldefinition_h3_p1
modeldefinition_h3_p3$expInfPrior <- inf_prior_3


seroModel_nih_h3_p1 <- createSeroJumpModel(data_titre_h3, known_exposure_h3, modeldefinition_h3_p1, TRUE)
seroModel_nih_h3_p2 <- createSeroJumpModel(data_titre_h3, known_exposure_h3, modeldefinition_h3_p2)
seroModel_nih_h3_p3 <- createSeroJumpModel(data_titre_h3, known_exposure_h3, modeldefinition_h3_p3)

seroModel_nih_h3_p1$data$knownExpVec %>% names
##########################################
############ Make models
##########################################

seroModel <- list(h1_p1 = seroModel_nih_h1_p1, h1_p2 = seroModel_nih_h1_p2, h1_p3 = seroModel_nih_h1_p3, h3_p1 = seroModel_nih_h3_p1, h3_p2 = seroModel_nih_h3_p2, h3_p3 = seroModel_nih_h3_p3)
saveRDS(seroModel, file = here::here("hpc", "nih_2024_inf", "nih_2024_model.RData"))
