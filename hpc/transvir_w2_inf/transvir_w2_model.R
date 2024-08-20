devtools::load_all()

library(lubridate)
library(patchwork)
library(tidybayes)
library(ggdist)

# Load data clean in the formet required for the model
gambia_pvnt_w2 <- get_data_titre_model_wave2_pvnt_iga()# This is the empirical prior for the exposure time
gambia_exp_w2 <- get_exposures_wave2() # This is the empirical prior for the exposure time
exp_prior_w2 <- get_exp_prior_wave2() # This is the empirical prior for the exposure time


known_exposure_gambia_exp_w2_add <- gambia_exp_w2 %>% mutate(type = "known")
sero_confirmed <-  gambia_pvnt_w2 %>% group_by(id) %>% mutate(r = row_number(), r_max = max(r)) %>% filter(r_max >=2) %>%
        select(!c(IgA, time )) %>% pivot_wider(names_from = r, values_from = sVNT) %>% 
        mutate(diff = `1` - `2`) %>% filter(diff < -0.903) %>% left_join(known_exposure_gambia_exp_w2_add) %>% 
        filter(is.na(type))


rate_inf_w2 <- (sero_confirmed %>% nrow) / (gambia_pvnt_w2$id %>% unique %>% length)


obsLogLikelihoodSerum = function(titre_val, titre_est, pars) {
    if (titre_val <= log10(40)) {
        ll <- pnorm(log10(40), titre_est, pars[1], log.p = TRUE)
    } else {
        ll <- dnorm(titre_val, titre_est, pars[1], log = TRUE)
    }
    ll
}

obsLogLikelihoodIgA = function(titre_val, titre_est, pars) {

    ll <- dnorm(titre_val, titre_est, pars[1], log = TRUE)

}


noInfSerumKinetics <- function(titre_est, timeSince, pars) {
    titre_est <- titre_est - pars[1] * (timeSince)
    titre_est <- max(0, titre_est)
}


#copFuncForm <- function(inf_status, esttitreExp, params, maxtitre) {
#    beta0 <- params[1]
#    mu <- params[3]
#    beta1 <- params[2]

#    p <- mu / (1.0 + exp(- (beta0 + beta1 * esttitreExp) ) )
#}

#copFuncFormInformed <- function(inf_status, esttitreExp, params, maxtitre) {
#    ep <- params[1]
#    beta1 <- params[2]
#    mu <- params[3]

#    beta0 <- log(0.1) - beta1 * maxtitre - ep
#    p <- mu / (1.0 + exp(- (beta0 + beta1 * esttitreExp) ) )
#}


#copLogLikelihood <- function(inf_status, esttitreExp, params, maxtitre) {
    # COP parameters
#    ep <- params[1]
#    beta1 <- params[2]
#    mu <- params[3]

#    beta0 <- log(0.1) - beta1 * maxtitre - ep
 #   p <- mu / (1.0 + exp(- (beta0 + beta1 * esttitreExp) ) )

 #   ll <- inf_status * log(p) + (1 - inf_status) * log(1 - p)
 #   ll
#}

copFuncFormGeneralised <- function(inf_status, esttitreExp, params, maxtitre ) {

    gamma <- params[1]
    k <- params[2]
    beta_mean <- params[3]
    beta_z <- params[4]
    beta_sigma <- params[5]
    alpha_mean <- params[6]
    alpha_z <- params[7]
    alpha_sigma <- params[8]

    beta <- inv.logit(beta_mean + beta_z * beta_sigma) * -3
    alpha <- inv.logit(alpha_mean + alpha_z * alpha_sigma) * (maxtitre * 0.75)

    r <- beta * (esttitreExp - alpha)
    p <- gamma * ((r / (1 + abs(r)^k)^(1 / k)) * 0.5 + 0.5)
    p
}

copLogLikelihood <- function(inf_status, esttitreExp, params, maxtitre ) {

    gamma <- params[1]
    k <- params[2]
    beta_mean <- params[3]
    beta_z <- params[4]
    beta_sigma <- params[5]
    alpha_mean <- params[6]
    alpha_z <- params[7]
    alpha_sigma <- params[8]

    beta <- inv.logit(beta_mean + beta_z * beta_sigma) * - 3
    alpha <- inv.logit(alpha_mean + alpha_z * alpha_sigma) * (maxtitre * 0.75)

    r <- beta * (esttitreExp - alpha)
    p <- gamma * ((r / (1 + abs(r)^k)^(1 / k)) * 0.5 + 0.5)
    ll <- inf_status * log(p) + (1 - inf_status) * log(1 - p)
    ll
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

    titre_est_log <- titre_est + log(titre_est_boost)  * max(0, 1 - titre_est * alpha)
    titre_est_log
}




# Define the biomarkers and exposure types in the model
biomarkers <- c("sVNT", "IgA")
exposureTypes <- c("none", "delta", "vax", "predelta")
exposureFitted <- "delta"

# Define the observational model
observationalModel <- list(
    names = c("sVNT", "IgA"),
    model = makeModel(
        addObservationalModel("sVNT", c("sigma"), obsLogLikelihoodSerum),
        addObservationalModel("IgA", c("sigma_a"), obsLogLikelihoodIgA)),
    prior = bind_rows(
        add_par_df("sigma", 0.0001, 4, "unif", 0.0001, 4),
        add_par_df("sigma_a", 0.0001, 4, "unif", 0.0001, 4),
        ) # observational model,
)

# Define the antibody kinetics model
abkineticsModel <- list(
    model = makeModel(
            addAbkineticsModel("none", "sVNT", "none",  c("wane"), noInfSerumKinetics),
            addAbkineticsModel("delta", "sVNT", "delta", c("y1_d", "t1_d", "r_d", "s"), infTuenisPower2016),
            addAbkineticsModel("vax", "sVNT", "vax", c("y1_vax", "t1_vax", "r_vax", "s"), infTuenisPower2016),
            addAbkineticsModel("predelta", "sVNT", "predelta", c("y1_pd", "t1_pd", "r_pd", "s"), infTuenisPower2016),
            addAbkineticsModel("none_a", "IgA", "none",  c("wane_a"), noInfSerumKinetics),
            addAbkineticsModel("delta_a", "IgA", "delta", c("y1_d_a", "t1_d_a", "r_d_a", "s_a"), infTuenisPower2016),
            addAbkineticsModel("vax_a", "IgA", "vax", c("y1_vax_a", "t1_vax_a", "r_vax_a", "s_a"), infTuenisPower2016),
            addAbkineticsModel("predelta_a", "IgA", "predelta", c("y1_pd_a", "t1_pd_a", "r_pd_a", "s_a"), infTuenisPower2016)
        ),
    prior = bind_rows(
        add_par_df("y1_d", 1, 6, "unif",  1, 6), # ab kinetics
        add_par_df("t1_d", 10, 25, "unif",  10, 25), # ab kinetics
        add_par_df("r_d", 1, 5, "unif", 1, 5), # ab kinetics 
        add_par_df("y1_vax", 1, 6, "unif",  1, 6), # ab kinetics
        add_par_df("t1_vax", 10, 25, "unif",  10, 25), # ab kinetics
        add_par_df("r_vax", 1, 5, "unif", 1, 5), # ab kinetics 
        add_par_df("y1_pd", 1, 9, "unif",  1, 9), # ab kinetics
        add_par_df("t1_pd", 10, 25, "unif",  10, 25), # ab kinetics
        add_par_df("r_pd", 1, 7, "unif", 1, 7), # ab kinetics 
        add_par_df("s",  0, 1, "unif", 0, 1), # ab kinetics 
        add_par_df("wane", 0.0, 0.01, "unif", 0.0, 0.01), # observational model
        add_par_df("y1_d_a", 1, 6, "unif",  1, 6), # ab kinetics
        add_par_df("t1_d_a", 25, 75, "unif",  25, 75), # ab kinetics
        add_par_df("r_d_a", 1, 5, "unif", 1, 5), # ab kinetics 
        add_par_df("y1_vax_a", 1, 6, "unif",  1, 6), # ab kinetics
        add_par_df("t1_vax_a", 25, 75, "unif",  25, 75), # ab kinetics
        add_par_df("r_vax_a", 1, 5, "unif", 1, 5), # ab kinetics 
        add_par_df("y1_pd_a", 1, 6, "unif",  1, 6), # ab kinetics
        add_par_df("t1_pd_a", 25, 75, "unif",  25, 75), # ab kinetics
        add_par_df("r_pd_a", 1, 7, "unif", 1, 7), # ab kinetics 
        add_par_df("s_a", 0, 1, "unif", 0, 1), # ab kinetics 
        add_par_df("wane_a", 0.0, 0.01, "unif", 0.0, 0.01) # observational model
    )
)


# modelW2_p1$data$max_titre
max_titre <- 3.72153 / 2

copModel <- list( 
        model = makeModel(
            addCopModel("sVNT", "delta", c("gamma", "k", "beta", "beta_z_1", "beta_sigma", "alpha", "alpha_z_1", "alpha_sigma"), copFuncFormGeneralised, copLogLikelihood),
            addCopModel("IgA", "delta", c("gamma", "k", "beta", "beta_z_2", "beta_sigma", "alpha", "alpha_z_2", "alpha_sigma"), copFuncFormGeneralised, copLogLikelihood)
        ),
        prior = bind_rows(
            add_par_df("gamma", 0, 1, "unif", 0, 1),
            add_par_df("k", 1, 3, "unif", 1, 3),
            add_par_df("beta", 0, 10, "norm", 0, 1.5),
            add_par_pool_df_non_centered("beta_z", 2, 0, "beta_sigma", "exp", 2, NA),
            add_par_df("alpha", -10, 10, "norm", 0, 1),
            add_par_pool_df_non_centered("alpha_z", 2, 0, "alpha_sigma", "exp", 2, NA)
        )
)

#copModel_old <- list( 
#        model = makeModel(
#            addCopModel("sVNT", "delta", c("beta0", "beta1", "mu"), copFuncForm, copLogLikelihood),
#            addCopModel("IgA", "delta", c("beta0_a", "beta1_a", "mu_a"), copFuncForm, copLogLikelihood)
#        ),
#        prior = bind_rows(
#            add_par_df("beta0", 0, 4, "norm", 1.76852, 1), # cop model (not used here)
#            add_par_df("beta1", -5, 0, "unif", -5, 0),
#            add_par_df("mu", 0.01, 1, "unif", 0.01, 1),
#            add_par_df("beta0_a", 0, 4, "norm", 1.76852, 1), # cop model (not used here)
#            add_par_df("beta1_a", -5, 0, "unif", -5, 0),
#            add_par_df("mu_a", 0.01, 1, "unif", 0.01, 1),
#        )
#)




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
    logPriorExpInf <- lfactorial(E_adj) + lfactorial(N_adj - E_adj) - lfactorial(N_adj ) + dbinom(E_adj, N_adj, 0.137931, log = TRUE)
    logPriorExpInf
}

modeldefinition_w2_p1 <- list(
    biomarkers = biomarkers,
    exposureTypes = exposureTypes,
    exposureFitted = exposureFitted,
    observationalModel = observationalModel,
    abkineticsModel = abkineticsModel,
    copModel = copModel,
    exposurePrior = exp_prior_w2,
    exposurePriorType = "empirical",
    expInfPrior = inf_prior_1
)




modeldefinition_w2_p2 <- modeldefinition_w2_p1
modeldefinition_w2_p2$expInfPrior <- inf_prior_2

modeldefinition_w2_p3 <- modeldefinition_w2_p1
modeldefinition_w2_p3$expInfPrior <- inf_prior_3


modelW2_p1 <- createSeroJumpModel(gambia_pvnt_w2, gambia_exp_w2, modeldefinition_w2_p1, TRUE)
modelW2_p2 <- createSeroJumpModel(gambia_pvnt_w2, gambia_exp_w2, modeldefinition_w2_p2)
modelW2_p3 <- createSeroJumpModel(gambia_pvnt_w2, gambia_exp_w2, modeldefinition_w2_p3)

#seroW2 <- list(p1 = modelW2_p1, p2 = modelW2_p2, p3 = modelW2_p3)
#saveRDS(seroW2, file = here::here("hpc", "transvir_w2_inf", "transvir_w2_model.RData"))



gambia_pvnt_w3 <- get_data_titre_model_wave3_pvnt_iga()# This is the empirical prior for the exposure time
gambia_exp_w3 <- get_exposures_wave3() # This is the empirical prior for the exposure time
exp_prior_w3 <- get_exp_prior_wave3() # This is the empirical prior for the exposure time


known_exposure_gambia_exp_w3_add <- gambia_exp_w3 %>% mutate(type = "known")
sero_confirmed <-  gambia_pvnt_w3 %>% group_by(id) %>% mutate(r = row_number(), r_max = max(r)) %>% filter(r_max >=2) %>%
        select(!c(IgA, time )) %>% pivot_wider(names_from = r, values_from = sVNT) %>% 
        mutate(diff = `1` - `2`) %>% filter(diff < -0.903) %>% left_join(known_exposure_gambia_exp_w3_add) %>% 
        filter(is.na(type))

rate_inf_w3 <- (sero_confirmed %>% nrow) / (gambia_pvnt_w3$id %>% unique %>% length)


# Define the biomarkers and exposure types in the model
biomarkers <- c("sVNT", "IgA")
exposureTypes <- c("none", "omicron", "vax")
exposureFitted <- "omicron"

# Define the observational model
observationalModel <- list(
    names = c("sVNT", "IgA"),
    model = makeModel(
        addObservationalModel("sVNT", c("sigma"), obsLogLikelihoodSerum),
        addObservationalModel("IgA", c("sigma_a"), obsLogLikelihoodIgA)),
    prior = bind_rows(
        add_par_df("sigma", 0.0001, 4, "unif", 0.0001, 4),
        add_par_df("sigma_a", 0.0001, 4, "unif", 0.0001, 4)
        ) # observational model,
)

# Define the antibody kinetics model
abkineticsModel <- list(
    model = makeModel(
            addAbkineticsModel("none", "sVNT", "none",  c("wane"), noInfSerumKinetics),
            addAbkineticsModel("omicron", "sVNT", "omicron", c("y1_o", "t1_o", "r_o", "s"), infTuenisPower2016),
            addAbkineticsModel("vax", "sVNT", "vax", c("y1_vax", "t1_vax", "r_vax", "s"), infTuenisPower2016),
            addAbkineticsModel("none_a", "IgA", "none",  c("wane_a"), noInfSerumKinetics),
            addAbkineticsModel("omicron_a", "IgA", "omicron", c("y1_o_a", "t1_o_a", "r_o_a", "s_a"), infTuenisPower2016),
            addAbkineticsModel("vax_a", "IgA", "vax", c("y1_vax_a", "t1_vax_a", "r_vax_a", "s_a"), infTuenisPower2016)
        ),
    prior = bind_rows(
        add_par_df("y1_o", 1, 6, "unif",  1, 6), # ab kinetics
        add_par_df("t1_o", 10, 25, "unif",  10, 25), # ab kinetics
        add_par_df("r_o", 1, 5, "unif", 1, 5), # ab kinetics 
        add_par_df("y1_vax", 1, 9, "unif",  1, 9), # ab kinetics
        add_par_df("t1_vax", 10, 25, "unif",  10, 25), # ab kinetics
        add_par_df("r_vax", 1, 5, "unif", 1, 5), # ab kinetics 
        add_par_df("s",  0, 1, "unif", 0, 1), # ab kinetics 
        add_par_df("wane", 0.0, 0.01, "unif", 0.0, 0.01), # observational model
        add_par_df("y1_o_a", 1, 6, "unif",  1, 6), # ab kinetics
        add_par_df("t1_o_a", 25, 75, "unif",  25, 75), # ab kinetics
        add_par_df("r_o_a", 1, 5, "unif", 1, 5), # ab kinetics 
        add_par_df("y1_vax_a", 1, 6, "unif",  1, 6), # ab kinetics
        add_par_df("t1_vax_a", 25, 75, "unif",  25, 75), # ab kinetics
        add_par_df("r_vax_a", 1, 5, "unif", 1, 5), # ab kinetics 
        add_par_df("s_a", 0, 1, "unif", 0, 1), # ab kinetics 
        add_par_df("wane_a", 0.0, 0.01, "unif", 0.0, 0.01) # observational model
    )
)




copModel <- list( 
        model = makeModel(
            addCopModel("sVNT", "delta", c("gamma", "k", "beta", "beta_z_1", "beta_sigma", "alpha", "alpha_z_1", "alpha_sigma"), copFuncFormGeneralised, copLogLikelihood),
            addCopModel("IgA", "delta", c("gamma", "k", "beta", "beta_z_2", "beta_sigma", "alpha", "alpha_z_2", "alpha_sigma"), copFuncFormGeneralised, copLogLikelihood)
        ),
        prior = bind_rows(
            add_par_df("gamma", 0, 1, "unif", 0, 1),
            add_par_df("k", 1, 3, "unif", 1, 3),
            add_par_df("beta", 0, 10, "norm", 0, 1.5),
            add_par_pool_df_non_centered("beta_z", 2, 0, "beta_sigma", "exp", 2, NA),
            add_par_df("alpha", -10, 10, "norm", 0, 1),
            add_par_pool_df_non_centered("alpha_z", 2, 0, "alpha_sigma", "exp", 2, NA)
        )
)
#copModel <- list( 
 #       model = makeModel(
  #          addCopModel("sVNT", "omicron", c("beta0", "beta1", "mu"), copFuncFormGeneralised, copLogLikelihood),
  #          addCopModel("IgA", "omicron", c("beta0_a", "beta1_a", "mu_a"), copFuncFormGeneralised, copLogLikelihood)
  #      ),
#        prior = bind_rows(
 #           add_par_df("beta0", 0, 4, "norm", 1.76852, 1), # cop model (not used here)
  #          add_par_df("beta1", -5, 0, "unif", -5, 0),
  #          add_par_df("mu", 0.01, 1, "unif", 0.01, 1),
   #         add_par_df("beta0_a", 0, 4, "norm", 1.76852, 1), # cop model (not used here)
#            add_par_df("beta1_a", -5, 0, "unif", -5, 0),
 #           add_par_df("mu_a", 0.01, 1, "unif", 0.01, 1),
  #      )
#)



   # add_par_df("d_vax_a", -1, 50, "unif",  -1, 50), # ab kinetics
   # add_par_df("e_vax_a", 1, 5, "unif",  1, 5), # ab kinetics
   # add_par_df("s_vax_a", 2, 5, "unif", 2, 5), # ab kinetics 
   # add_par_df("h_vax_a", 0.5, 5, "unif", 0.5,  3), # ab kinetics 
   # add_par_df("d_pd_a", -1, 50, "unif",  -1, 50), # ab kinetics
   # add_par_df("e_pd_a", 1, 5, "unif",  1, 5), # ab kinetics
   # add_par_df("s_pd_a", 2, 5, "unif", 2, 5), # ab kinetics 
   # add_par_df("h_pd_a", 0.5, 5, "unif", 0.5,  3), # ab kinetics 
   # add_par_df("d_d_a", -1, 50, "unif",  -1, 50), # ab kinetics
   # add_par_df("e_d_a", 1, 5, "unif",  1, 5), # ab kinetics
   # add_par_df("s_d_a", 2, 5, "unif", 2, 5), # ab kinetics 
   # add_par_df("h_d_a", 0.5, 5, "unif", 0.5,  3), # ab kinetics 


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
    logPriorExpInf <- lfactorial(E_adj) + lfactorial(N_adj - E_adj) - lfactorial(N_adj ) + dbinom(E_adj, N_adj, 0.1395349, log = TRUE)
    logPriorExpInf
}


modeldefinition_w3_p1 <- list(
    biomarkers = biomarkers,
    exposureTypes = exposureTypes,
    exposureFitted = exposureFitted,
    observationalModel = observationalModel,
    abkineticsModel = abkineticsModel,
    copModel = copModel,
    exposurePrior = exp_prior_w3,
    exposurePriorType = "empirical",
    expInfPrior = inf_prior_1
)

modeldefinition_w3_p2 <- modeldefinition_w3_p1
modeldefinition_w3_p2$expInfPrior <- inf_prior_2

modeldefinition_w3_p3 <- modeldefinition_w3_p1
modeldefinition_w3_p3$expInfPrior <- inf_prior_3


modelW3_p1 <- createSeroJumpModel(gambia_pvnt_w3, gambia_exp_w3, modeldefinition_w3_p1, TRUE)
modelW3_p2 <- createSeroJumpModel(gambia_pvnt_w3, gambia_exp_w3, modeldefinition_w3_p2)
modelW3_p3 <- createSeroJumpModel(gambia_pvnt_w3, gambia_exp_w3, modeldefinition_w3_p3)

seroW2W3 <- list(w2_p1 = modelW2_p1, w2_p2 = modelW2_p2, w2_p3 = modelW2_p3, w3_p1 = modelW3_p1, w3_p2 = modelW3_p2, w3_p3 = modelW3_p3)
saveRDS(seroW2W3, file = here::here("hpc", "transvir_w2_inf", "transvir_w23_model.RData"))
