library(serosim)

## Load additional packages required 
library(tidyverse)
library(data.table)
library(ggplot2)
library(patchwork)
library(reshape2)
library(extraDistr)
####################################
####### Titre-dep boost and no wane #######
####################################


library(finalsize)
r0 <- 1.5
uk_pop <- 67 * 1e6
contact_matrix <- matrix(1.0) / uk_pop
susceptibility_full <- matrix(1)
susceptibility <- matrix(0.5)
p_susceptibility <- matrix(1)
# calculate final size
final_size_data <- final_size(
  r0 = r0,
  contact_matrix = contact_matrix,
  demography_vector = uk_pop,
  susceptibility = susceptibility_full,
  p_susceptibility = p_susceptibility
)

# view the output data frame
final_size_data

biomarker_protection <- function(biomarker_quantity, biomarker_prot_midpoint, biomarker_prot_width) {
    risk <- 1 - 1/(1 + exp(biomarker_prot_width * (biomarker_quantity - biomarker_prot_midpoint)))
    return(risk)
}

get_foe_par <- function(modeli) {
  biomarker_map_original <- tibble(exposure_id=c("covid_ifxn"),biomarker_id=c("covid_svnt"))
  modeli$simpar$biomarker_map <- reformat_biomarker_map(biomarker_map_original)

  modeli$simpar$exp <- array(0, dim=c(modeli$simpar$N, max(modeli$simpar$times), n_distinct(modeli$simpar$biomarker_map$exposure_id)))
  modeli$simpar$foe_pars <- array(0, dim=c(modeli$simpar$N, max(modeli$simpar$times), n_distinct(modeli$simpar$biomarker_map$exposure_id)))
  for (i in 1:modeli$simpar$N) {
    cat(i, "\n")
    if (runif(1, 0, 1) < 0.6) {
      d <- max(min(round(rnorm(1, 60, 20), 0), 120 - 21), 7)
      modeli$simpar$exp[i, d, 1] <- 1
      if (runif(1, 0, 1) > biomarker_protection(modeli$simpar$start_titre[i], 2, modeli$protection_curve)) {
        modeli$simpar$foe_pars[i, d, 1] <- 1
      }
    }
  }
  modeli$simpar$N_exp <- modeli$simpar$exp %>% sum
  modeli$simpar$N_inf <- modeli$simpar$foe_pars %>% sum

  modeli
}

get_start_titre <- function(simpar) {
  simpar$start_titre <- pmin(runif(simpar$N, 0, 4), 4)
  simpar
}

T <- 120

source("R/sim_data_utils.R")

run_uncert_model_bi <- function(modeli, uncert_vec) {

  uncert_vec %>% 
    map( 
    function(uncert) {
    N <- modeli$simpar$N
    # 1. Simulation settings
    ## Specify the number of time periods to simulate 
    ## Set simulation settings
    dir.create(here::here("outputs", "sim_data", modeli$name), showWarnings = FALSE)
    simulation_settings <- list("t_start"=1,"t_end"=max(modeli$simpar$times))

    # 2. Population demography
    ## Specify the number of individuals in the simulation 
    ## Generate the population demography tibble
    demography <- generate_pop_demography(N, modeli$simpar$times, age_min=0, birth_times = 0, prob_removal=0)
    # 3. Exposure to biomarker mapping

    ## Specify a simple exposure model which calculates the probability of exposure directly from the force of exposure at that time step
    exposure_model<-exposure_model_indiv_fixed
    start_titre <- modeli$simpar$start_titre
    foe_pars <- modeli$simpar$foe_par
    biomarker_map <- modeli$simpar$biomarker_map
    immune_histories_fixed <- foe_pars

    ## Examine the probability of exposure over time for the specified exposure model
    #plot_exposure_model(exposure_model=exposure_model_simple_FOE, times=times, n_groups = 1, n_exposures = 1, foe_pars=foe_pars)

    # 5. Immunity model
    immunity_model<-immunity_model_ifxn_biomarker_prot
    # COP
    #plot_biomarker_mediated_protection(0:4, biomarker_prot_midpoint=2, biomarker_prot_width=1)

    # 6.  Antibody model, antibody kinetics parameters, and `draw_parameters
    antibody_model<-antibody_model_exp_asymp

    a = modeli$simpar$a
    b = modeli$simpar$b
    c = modeli$simpar$c
    wane_long = modeli$simpar$wane_long

    biomarker_ceiling_gradient = modeli$simpar$biomarker_ceiling_gradient
    biomarker_ceiling_threshold = modeli$simpar$biomarker_ceiling_threshold
    uncert = uncert

    model_pars_original_1 <- data.frame(
      exposure_id = rep("covid_ifxn", 4),
      biomarker_id = rep("covid_svnt", 4),
      name = c("a", "b", "c", "wane_short"),
      mean = c(a, b, c, 0),
      sd = c(a * 0.01, b * 0.01, c * 0.01,  0),
      distribution = c("normal", "normal", "normal", "normal")
    )

    model_pars_original_2 <- data.frame(
        exposure_id = c(rep("covid_ifxn", 4), NA),
        biomarker_id = rep("covid_svnt", 5),
        name = c("biomarker_ceiling_threshold", "biomarker_ceiling_gradient", 
            "biomarker_prot_midpoint", "biomarker_prot_width", "obs_sd"),
        mean = c(biomarker_ceiling_threshold, biomarker_ceiling_gradient, 2, 1, NA),
        sd = c(rep(NA, 4), uncert),
        distribution = c(rep("", 4), "normal")
    )

    model_pars_original <- bind_rows(model_pars_original_1, model_pars_original_2)

    ## Reformat model_pars for runserosim
    model_pars<-reformat_biomarker_map(model_pars_original)
    draw_parameters<-draw_parameters_random_fx
    #plot_biomarker_dependent_boosting(start=0, end=4, by=0.1, biomarker_ceiling_threshold=biomarker_ceiling_threshold, biomarker_ceiling_gradient=biomarker_ceiling_gradient)

    # 7. Observation model and `observation_times`
    bounds_1<-tibble(biomarker_id=c(1),name=c("lower_bound", "upper_bound"),value=c(0,Inf))
    bounds <- bind_rows(bounds_1)
    observation_model<-observation_model_continuous_bounded_noise

    ## Specify assay sensitivity and specificity needed for the observation model
    sensitivity<-1
    specificity<-1
    observation_times <- modeli$observation_times

    VERBOSE<-100

    ## Run the core simulation and save outputs in "res"
    res <- runserosim_edit(
      simulation_settings,
      demography,
      observation_times,
      foe_pars,
      biomarker_map,
      model_pars,
      exposure_model,
      immunity_model,
      antibody_model,
      observation_model,
      draw_parameters,

      ## Specify other arguments needed
      VERBOSE=VERBOSE,
      bounds=bounds,
      max_events=1,
      
      immune_histories_fixed = immune_histories_fixed,
      start_titre = start_titre
    )

    dir.create(here::here("outputs", "sim_data", modeli$name), showWarnings = FALSE)
    saveRDS(res, file = here::here("outputs", "sim_data", modeli$name, paste0("sim_data_", round(uncert, 2), ".rds")))

    }
  )
}


plot_sim <- function(modeli) {

  uncert_vec %>% 
      map( 
      function(uncert) {
      res <- readRDS(file = here::here("outputs", "sim_data", modeli$name, paste0("sim_data_", round(uncert, 1), ".rds")))

  #plot_subset_individuals_history(res$biomarker_states, res$immune_histories_long, subset=30, demography)
  #plot_biomarker_quantity(res$biomarker_states) 
  #plot_obs_biomarkers_paired_sample(res$observed_biomarker_states)

    N <- modeli$simpar$N
    exp_times_sim <- modeli$simpar$exp %>% as.data.frame %>% 
      mutate(i = 1:N) %>% pivot_longer(!i, names_to = "t", values_to = "exp") %>% mutate(t = as.numeric(substr(t, 2, 4))) %>% 
      filter(exp == 1) %>% complete(i = 1:N, fill = list(t = NA)) %>% arrange(t)
    ids_order <- exp_times_sim %>% pull(i)
    exp_times_sim <- exp_times_sim %>% mutate(i = factor(i, levels = rev(ids_order))) 
    inf_times_sim <- res$immune_histories_long %>% filter(value == 1) %>% complete(i = 1:N, fill = list(t = NA)) %>% arrange(t) %>% 
      mutate(i = as.character(i)) %>% rename(inf = value)
    epi_times_sims <- left_join(exp_times_sim, inf_times_sim) %>% mutate(
      exp_type = 
      case_when(
        (exp == 1 & is.na(inf)) ~ "Exposed and not infected",
        (exp == 1 & inf == 1)~"Exposed and infected"
      )
    ) %>% filter(!is.na(exp_type))

    N_exp = (epi_times_sims$exp == 1) %>% sum
    N_inf =(!is.na(epi_times_sims$inf == 1)) %>% sum


    start_bleed <- res$observed_biomarker_states %>% as.data.frame %>% group_by(i) %>% filter(t == min(t)) %>%
      mutate(i = factor(i, levels = rev(ids_order)))
    end_bleed <- res$observed_biomarker_states %>% as.data.frame %>% group_by(i) %>% filter(t == max(t)) %>%
      mutate(t = t - 7) %>%
      mutate(i = factor(i, levels = rev(ids_order)))

    p1 <- res$observed_biomarker_states %>% mutate(i = factor(i, levels = rev(ids_order))) %>%
      ggplot() + 
        geom_tile(data = start_bleed, aes(y = i, x = t / 2, width = t ), alpha = 0.9, fill = "gray90") +
        geom_tile(data = end_bleed, aes(y = i, x = t + (T - t) / 2, width = (T - t) ), alpha = 0.9, fill = "gray90") +
        geom_point(aes(x = t, y = i), shape = "|", size = 2.5) + 
        geom_point(data = epi_times_sims, aes(x = t, y = i, color = exp_type) ,  shape = 4) + theme_bw() + 
        labs(x = "Time in study (days)", y = "Id of individual", color = "Exposure type") + 
        geom_hline(yintercept = seq(0.5, N + 0.5, 1), color = "gray90") +
        theme(panel.grid.major.y = element_blank())

    p2 <- data.frame(start_model = res$start_titre) %>% 
      pivot_longer(everything(), names_to = "type", values_to = "value") %>%
      ggplot() + geom_density(aes(value), fill = "blue", alpha = 0.2) + 
      labs(x = "Titre at start of study", y = "Density") + theme_bw() 

    p1 / p2 + plot_layout(heights = c(3, 1))

    dir.create(here::here("outputs", "sim_data", modeli$name, "figs"))
    dir.create(here::here("outputs", "sim_data", modeli$name, "figs", paste0("obs_", round(uncert, 1))))
    ggsave(here::here("outputs", "sim_data", modeli$name, "figs", paste0("obs_", round(uncert, 1)), "simulated_summary.png"), height = 15, width = 15)
    }
  )
}

# Run simulations for a continuous epidemic surveillance study (CES)
# Run simulations for a Pre- and Post-epidemic serosurvey (PPES)

uncert_vec <- c(0.01, seq(0.05, 1, 0.05))


simpar <- list(
  N = 200,
  T = T,
  S = 5,
  a = 1.5, # boosting after infection
  b = 0.3, # boosting after infection
  c = 1,
  wane_long = 0.005,
  # Titre-dep boosting
  biomarker_ceiling_gradient = 0, # decline per sVNT unit
  biomarker_ceiling_threshold = 4, # fixed
  times = seq(1, T, by=1)
)
simpar <- simpar %>% get_start_titre  # add foe_pars, biomarker_map, start_titre


observation_times_ces <-  
  1:simpar$N %>% 
    map_df( ~tibble(i=.x, t = floor(
      c(runif(1, 1, 7), runif(simpar$S - 2, 8, simpar$T - 15), runif(1, 106, 120))
      ), b = 1) )

observation_times_ppes <-  
  1:simpar$N %>% 
    map_df( ~tibble(i=.x, t = floor(
      c(runif(1, 1, 7),  runif(1, 106, 120))
      ), b = 1) )

modelcesNoCOP <- list(
  name = "cesNoCOP_notd",
  observation_times = observation_times_ces,
  protection_curve = 0,
  simpar = simpar,
  uncert_vec = uncert_vec
) %>% get_foe_par

modelcesCOP <- list(
  name = "cesCOP_notd",
  observation_times = observation_times_ces,
  protection_curve = 2,
  simpar = simpar,
  uncert_vec = uncert_vec
) %>% get_foe_par

modelppesNoCOP <- list(
  name = "ppesNoCOP_notd",
  observation_times = observation_times_ppes,
  protection_curve = 0,
  simpar = simpar,
  uncert_vec = uncert_vec
) %>% get_foe_par

modelppesCOP <- list(
  name = "ppesCOP_notd",
  observation_times = observation_times_ppes,
  protection_curve = 2,
  simpar = simpar,
  uncert_vec = uncert_vec
) %>% get_foe_par



# Run both simulated data
run_uncert_model_bi(modelcesNoCOP, uncert_vec)
run_uncert_model_bi(modelcesCOP, uncert_vec)
# Run both simulated data
run_uncert_model_bi(modelppesNoCOP, uncert_vec)
run_uncert_model_bi(modelppesCOP, uncert_vec) 

saveRDS(modelcesNoCOP, here::here("outputs", "sim_data", modelcesNoCOP$name, "inputs.RDS"))
saveRDS(modelcesCOP, here::here("outputs", "sim_data", modelcesCOP$name, "inputs.RDS"))
saveRDS(modelppesNoCOP, here::here("outputs", "sim_data", modelppesNoCOP$name, "inputs.RDS"))
saveRDS(modelppesCOP, here::here("outputs", "sim_data", modelppesCOP$name, "inputs.RDS"))