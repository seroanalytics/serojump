devtools::load_all()

library(lubridate)
library(patchwork)
library(tidybayes)
library(ggdist)

# Load data clean in the formet required for the model
gambia_pvnt_w2 <- get_data_titre_model_wave2_pvnt_iga()# This is the empirical prior for the exposure time
gambia_exp_w2 <- get_exposures_wave2() # This is the empirical prior for the exposure time
exp_prior_w2 <- get_exp_prior_wave2() # This is the empirical prior for the exposure time


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

a <- 4
b <- 0.03
c <- 2
data.frame(
    t = 0:100,
    x = log(exp(a) * exp(-b/10 * (0:100 - 14)) + exp(c))
) %>% 
    ggplot() + geom_line(aes(t, x))


infIgAinetics <- function(titre_est, timeSince, pars) {
    d <- pars[1]
    e <- pars[2]
    s <- pars[3]
    h <- pars[4]
    if (timeSince > 249) {
        timeSince <- 250 - 1
    }

    if (timeSince == 0) {
        boost <- 0
    } else {
        x <- timeSince / 250

        sub <- ((1 - x^s) / (d*x^s + 1))^e

        num <- (d + 1) * e * s * x^(s - 1) * sub
        denom <- (1 - x^s) * (d*x^s + 1)
        boost <- num / denom * h
        if (boost < 0) {
            boost <- 0
        } else if (boost > 8){
            boost <- 8
        }
    }
    titre_est + boost
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
            addAbkineticsModel("delta", "sVNT", "delta", c("a_d", "b_d", "c_d"), infSerumKinetics),
            addAbkineticsModel("vax", "sVNT", "vax", c("a_vax", "b_vax", "c_vax"), infSerumKinetics),
            addAbkineticsModel("predelta", "sVNT", "predelta", c("a_pd", "b_pd", "c_pd"), infSerumKinetics),
            addAbkineticsModel("none_a", "IgA", "none",  c("wane_a"), noInfSerumKinetics),
            addAbkineticsModel("delta_a", "IgA", "delta", c("d_d_a", "e_d_a", "s_d_a", "h_d_a"), infIgAinetics),
            addAbkineticsModel("vax_a", "IgA", "vax", c("d_vax_a", "e_vax_a", "s_vax_a", "h_vax_a"), infIgAinetics),
            addAbkineticsModel("predelta_a", "IgA", "predelta", c("d_pd_a", "e_pd_a", "s_pd_a", "h_pd_a"), infIgAinetics)
        ),
    prior = bind_rows(
        add_par_df("a_vax", -6, 6, "norm",  0, 2), # ab kinetics
        add_par_df("b_vax", 0, 1, "norm",  0.3, 0.05), # ab kinetics
        add_par_df("c_vax", 0, 4, "unif", 0,  4), # ab kinetics 
        add_par_df("a_pd", -6, 6, "norm",  0, 2), # ab kinetics
        add_par_df("b_pd", 0, 1, "norm",  0.3, 0.05), # ab kinetics
        add_par_df("c_pd", 0, 4, "unif", 0,  4), # ab kinetics
        add_par_df("a_d", -6, 6, "norm",  0, 2), # ab kinetics
        add_par_df("b_d", 0, 1, "norm",  0.3, 0.05), # ab kinetics
        add_par_df("c_d", 0, 4, "unif", 0,  4), # ab kinetics
        add_par_df("wane", 0.0, 0.01, "unif", 0.0, 0.01), # observational model
        add_par_df("d_vax_a", -1, 50, "unif",  -1, 50), # ab kinetics
        add_par_df("e_vax_a", 1, 5, "unif",  1, 5), # ab kinetics
        add_par_df("s_vax_a", 2, 5, "unif", 2, 5), # ab kinetics 
        add_par_df("h_vax_a", 0.5, 5, "unif", 0.5,  3), # ab kinetics 
        add_par_df("d_pd_a", -1, 50, "unif",  -1, 50), # ab kinetics
        add_par_df("e_pd_a", 1, 5, "unif",  1, 5), # ab kinetics
        add_par_df("s_pd_a", 2, 5, "unif", 2, 5), # ab kinetics 
        add_par_df("h_pd_a", 0.5, 5, "unif", 0.5,  3), # ab kinetics 
        add_par_df("d_d_a", -1, 50, "unif",  -1, 50), # ab kinetics
        add_par_df("e_d_a", 1, 5, "unif",  1, 5), # ab kinetics
        add_par_df("s_d_a", 2, 5, "unif", 2, 5), # ab kinetics 
        add_par_df("h_d_a", 0.5, 5, "unif", 0.5,  3), # ab kinetics 
        add_par_df("wane_a", 0.0, 0.01, "unif", 0.0, 0.01) # observational model
    )
)


inf_prior_1 <- function(N, E, I, K) {
    0
}

inf_prior_2 <- function(N, E, I, K) {
    N_adj <- N - K
    E_adj <- E - K 
    logPriorExpInf <- lfactorial(E_adj) + lfactorial(N_adj - E_adj) - lfactorial(N_adj ) #+ log(1 / N_adj)
    logPriorExpInf
}

inf_prior_3 <- function(N, E, I, K) {
    N_adj <- N - K
    E_adj <- E - K
    logPriorExpInf <- lfactorial(E_adj) + lfactorial(N_adj - E_adj) - lfactorial(N_adj ) + dbinom(E_adj, N_adj, 0.01, log = TRUE)
    logPriorExpInf
}

modeldefinition_p1 <- list(
    biomarkers = biomarkers,
    exposureTypes = exposureTypes,
    exposureFitted = exposureFitted,
    observationalModel = observationalModel,
    abkineticsModel = abkineticsModel,
    exposurePrior = exp_prior_w2,
    exposurePriorType = "empirical",
    expInfPrior = inf_prior_1
)

modeldefinition_p2 <- modeldefinition_p1
modeldefinition_p2$expInfPrior <- inf_prior_2

modeldefinition_p3 <- modeldefinition_p1
modeldefinition_p3$expInfPrior <- inf_prior_3


modelW2_p1 <- createSeroJumpModel(gambia_pvnt_w2, gambia_exp_w2, modeldefinition_p1)
modelW2_p2 <- createSeroJumpModel(gambia_pvnt_w2, gambia_exp_w2, modeldefinition_p2)
modelW2_p3 <- createSeroJumpModel(gambia_pvnt_w2, gambia_exp_w2, modeldefinition_p3)

seroW2 <- list(p1 = modelW2_p1, p2 = modelW2_p2, p3 = modelW2_p3)
saveRDS(seroW2, file = here::here("hpc", "transvir_w2_inf", "transvir_w2_model.RData"))
