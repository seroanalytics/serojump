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

copFuncForm <- function(inf_status, esttitreExp, params, maxtitre) {
    beta0 <- params[1]
    beta1 <- params[2]
    mu <- params[3]

    p <- mu / (1.0 + exp(- (beta0 + beta1 * esttitreExp) ) )
}


copFuncFormInformed <- function(inf_status, esttitreExp, params, maxtitre) {
    ep <- params[1]
    beta1 <- params[2]
    mu <- params[3]

    beta0 <- log(0.1) - beta1 * maxtitre - ep
    p <- mu / (1.0 + exp(- (beta0 + beta1 * esttitreExp) ) )
}


copLogLikelihood <- function(inf_status, esttitreExp, params, maxtitre) {
    # COP parameters
    ep <- params[1]
    beta1 <- params[2]
    mu <- params[3]

    beta0 <- log(0.1) - beta1 * maxtitre - ep
    p <- mu / (1.0 + exp(- (beta0 + beta1 * esttitreExp) ) )

    ll <- inf_status * log(p) + (1 - inf_status) * log(1 - p)
    ll
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
            addAbkineticsModel("vax", "A/Sydney/5/2021", "vax", c("y1_vax", "t1_vax", "r_vax", "s"), infTuenisPower2016),
            addAbkineticsModel("h1_2023","A/Sydney/5/2021", "h1_2023", c("y1_v", "t1_v", "r_v", "s"), infTuenisPower2016),
            addAbkineticsModel("vax_e", "A/Sydney/5/2021e", "vax", c("y1_vax_e", "t1_vax_e", "r_vax_e", "s_e"), infTuenisPower2016),
            addAbkineticsModel("h1_2023_e", "A/Sydney/5/2021e", "h1_2023", c("y1_v_e", "t1_v_e", "r_v_e", "s_e"), infTuenisPower2016)
        ),
    prior = bind_rows(
        add_par_df("y1_vax", 1, 6, "unif",  1, 6), # ab kinetics
        add_par_df("t1_vax", 7, 24, "norm", 14, 3), # ab kinetics
        add_par_df("r_vax", 1, 5, "unif", 1, 5), # ab kinetics 
        add_par_df("y1_v", 1, 9, "unif",  1, 9), # ab kinetics
        add_par_df("t1_v", 7, 24, "norm", 14, 3), # ab kinetics
        add_par_df("r_v", 1, 5, "unif", 1, 5), # ab kinetics 
        add_par_df("s", 0, 1, "unif", 0, 1), # ab kinetics 
        add_par_df("y1_vax_e", 1, 6, "unif",  1, 6), # ab kinetics
        add_par_df("t1_vax_e", 7, 24, "norm", 14, 3), # ab kinetics
        add_par_df("r_vax_e", 1, 5, "unif", 1, 5), # ab kinetics 
        add_par_df("y1_v_e", 1, 6, "unif",  1, 6), # ab kinetics
        add_par_df("t1_v_e", 7, 24, "norm", 14, 3), # ab kinetics
        add_par_df("r_v_e", 1, 5, "unif", 1, 5), # ab kinetics 
        add_par_df("s_e", 0, 1, "unif", 0, 1) # ab kinetics 
    )
)

copModel_h1 <- list( 
        model = makeModel(
            addCopModel("A/Sydney/5/2021", "h1_2023", c("beta0", "beta1", "mu"), copFuncForm, copLogLikelihood),
            addCopModel("A/Sydney/5/2021e", "h1_2023", c("beta0_e", "beta1_e", "mu_e"), copFuncForm, copLogLikelihood)
        ),
        prior = bind_rows(
            add_par_df("beta0", 0, 10, "norm", 4.5, 2.5), # cop model (not used here)
            add_par_df("beta1", -5, 0, "unif", -5, 0),
            add_par_df("mu", 0.01, 0.5, "unif", 0.01, 0.5),
            add_par_df("beta0_e", 0, 10, "norm",  4.5, 2.5), # cop model (not used here)
            add_par_df("beta1_e", -5, 0, "unif", -5, 0),
            add_par_df("mu_e", 0.01, 0.5, "unif", 0.01, 0.5),
        )
)

          #  add_par_df("beta0", 0, 10, "norm", 5, 2.5), # cop model (not used here)
           # add_par_df("beta1", -5, 0, "unif", -5, 0),
            #add_par_df("mu",  0.01, 0.5, "unif", 0.01, 0.5),

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
    copModel = NULL,
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
            addAbkineticsModel("vax", "A/Darwin/06/2021", "vax", c("y1_vax", "t1_vax", "r_vax", "s"), infTuenisPower2016),
            addAbkineticsModel("h3_2023","A/Darwin/06/2021", "h3_2023", c("y1_v", "t1_v", "r_v", "s"), infTuenisPower2016),
            addAbkineticsModel("vax_e", "A/Darwin/09/2021e", "vax", c("y1_vax_e", "t1_vax_e", "r_vax_e", "s_e"), infTuenisPower2016),
            addAbkineticsModel("h3_2023_e", "A/Darwin/09/2021e", "h3_2023", c("y1_v_e", "t1_v_e", "r_v_e", "s_e"), infTuenisPower2016)
        ),
    prior = bind_rows(
        add_par_df("y1_vax", 1, 6, "unif",  1, 6), # ab kinetics
        add_par_df("t1_vax", 7, 24, "norm", 14, 3), # ab kinetics
        add_par_df("r_vax", 1, 5, "unif", 1, 5), # ab kinetics 
        add_par_df("y1_v", 1, 9, "unif",  1, 9), # ab kinetics
        add_par_df("t1_v", 7, 24, "norm", 14, 3), # ab kinetics
        add_par_df("r_v", 1, 5, "unif", 1, 5), # ab kinetics 
        add_par_df("s", 0, 1, "unif", 0, 1), # ab kinetics 
        add_par_df("y1_vax_e", 1, 6, "unif",  1, 6), # ab kinetics
        add_par_df("t1_vax_e", 7, 24, "norm", 14, 3), # ab kinetics
        add_par_df("r_vax_e", 1, 5, "unif", 1, 5), # ab kinetics 
        add_par_df("y1_v_e", 1, 6, "unif",  1, 6), # ab kinetics
        add_par_df("t1_v_e", 7, 24, "norm", 14, 3), # ab kinetics
        add_par_df("r_v_e", 1, 5, "unif", 1, 5), # ab kinetics 
        add_par_df("s_e", 0, 1, "unif", 0, 1) # ab kinetics 
    )
)

copModel_h3 <- list( 
        model = makeModel(
            addCopModel("A/Darwin/06/2021", "h3_2023", c("beta0", "beta1", "mu"), copFuncForm, copLogLikelihood),
            addCopModel("A/Darwin/09/2021e", "h3_2023", c("beta0_e", "beta1_e", "mu_e"), copFuncForm, copLogLikelihood)
        ),
        prior = bind_rows(
            add_par_df("beta0", 0, 10, "norm", 4, 2.5), # cop model (not used here)
            add_par_df("beta1", -5, 0, "unif", -5, 0),
            add_par_df("mu", 0.01, 0.5, "unif", 0.01, 0.5),
            add_par_df("beta0_e", 0, 10, "norm",  5, 2.5), # cop model (not used here)
            add_par_df("beta1_e", -5, 0, "unif", -5, 0),
            add_par_df("mu_e", 0.01, 0.5, "unif", 0.01, 0.5),
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
    copModel = NULL,
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

seroModel_nih_h3_p1$data$max_titre

##########################################
############ Make models
##########################################

seroModel <- list(h1_p1 = seroModel_nih_h1_p1, h1_p2 = seroModel_nih_h1_p2, h1_p3 = seroModel_nih_h1_p3, h3_p1 = seroModel_nih_h3_p1, h3_p2 = seroModel_nih_h3_p2, h3_p3 = seroModel_nih_h3_p3)
saveRDS(seroModel, file = here::here("hpc", "nih_2024_inf", "nih_2024_model.RData"))
