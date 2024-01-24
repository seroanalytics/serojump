library(serosim)

## Load additional packages required 
library(tidyverse)
library(data.table)
library(ggplot2)
library(patchwork)
library(reshape2)
library(extraDistr)
require(tidybayes)

####################################
####### Titre-dep boost and no wane #######
####################################
uncert_vec <- c(0.1, 0.3, 0.5)


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

    df_titres <- data.frame(i = as.character(1:N), start_model = res$start_titre)
    epi_times_sims <- left_join(epi_times_sims, df_titres)
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

    p2 <- epi_times_sims %>% 
      ggplot() + geom_density(aes(start_model), fill = "blue", alpha = 0.2) + 
      labs(x = "Titre at start of study", y = "Density") + theme_bw() 

    p3 <- epi_times_sims %>% replace_na(list(inf = 0)) %>% 
      ggplot() + geom_point(aes(start_model, inf), fill = "red", alpha = 0.2) + 
      geom_smooth(aes(start_model, inf)) +
      labs(x = "Titre at start of study", y = "Infection status") + theme_bw() 

    p1 / (p2 + p3) + plot_layout(heights = c(3, 1))

    dir.create(here::here("outputs", "sim_data", modeli$name, "figs"))
    dir.create(here::here("outputs", "sim_data", modeli$name, "figs", paste0("obs_", round(uncert, 1))))
    ggsave(here::here("outputs", "sim_data", modeli$name, "figs", paste0("obs_", round(uncert, 1)), "simulated_summary.png"), height = 15, width = 15)
    }
  )
}


modelcesNoCOP <- readRDS(here::here("outputs", "sim_data", "cesNoCOP", "inputs.RDS"))
modelcesCOP <- readRDS(here::here("outputs", "sim_data", "cesCOP", "inputs.RDS"))
modelppesNoCOP <- readRDS(here::here("outputs", "sim_data", "ppesNoCOP", "inputs.RDS"))
modelppesCOP <- readRDS(here::here("outputs", "sim_data", "ppesCOP", "inputs.RDS"))

# Plot sim everything
plot_sim(modelcesNoCOP)
plot_sim(modelcesCOP)
plot_sim(modelppesNoCOP)
plot_sim(modelppesCOP)


################################################################################
# THIS SECTION PRODUCES FIGURES FOR THE SIMULATED DATA VERSION OF MANUSCRIPT #
################################################################################


# COLORS!
## #882255, Exposed but not infected
## #DDCC77, Infected
## #888888, Model A COP
## #661100, Model B COP

## PLOT A: Epidemiology and survery design
modeli <- modelcesNoCOP
res <- readRDS(file = here::here("outputs", "sim_data", modeli$name, paste0("sim_data_", round(0.1, 1), ".rds")))

res$kinetics_parameters

N <- modeli$simpar$N
exp_times_sim <- modeli$simpar$exp %>% as.data.frame %>% 
  mutate(i = 1:N) %>% pivot_longer(!i, names_to = "t", values_to = "exp") %>% mutate(t = as.numeric(substr(t, 2, 4))) %>% 
  filter(exp == 1) %>% complete(i = 1:N, fill = list(t = NA)) %>% arrange(t)
ids_order <- exp_times_sim %>% pull(i)
exp_times_sim <- exp_times_sim %>% mutate(i = factor(i, levels = rev(ids_order))) %>% arrange(t)
inf_times_sim <- res$immune_histories_long %>% filter(value == 1) %>% complete(i = 1:N, fill = list(t = NA)) %>% arrange(t) %>% 
  mutate(i = as.character(i)) %>% rename(inf = value)
epi_times_sims <- left_join(exp_times_sim, inf_times_sim) %>% mutate(
  exp_type = 
  case_when(
    (exp == 1 & is.na(inf)) ~ "Exposed and not infected",
    (exp == 1 & inf == 1)~"Exposed and infected"
  )
) %>% filter(!is.na(exp_type))

N_exp <- (epi_times_sims$exp == 1) %>% sum
N_inf <- (!is.na(epi_times_sims$inf == 1)) %>% sum


start_bleed <- res$observed_biomarker_states %>% as.data.frame %>% group_by(i) %>% filter(t == min(t)) %>%
  mutate(i = factor(i, levels = rev(ids_order)))
end_bleed <- res$observed_biomarker_states %>% as.data.frame %>% group_by(i) %>% filter(t == max(t)) %>%
  mutate(t = t - 7) %>%
  mutate(i = factor(i, levels = rev(ids_order)))

p1 <- res$observed_biomarker_states %>% mutate(i = factor(i, levels = rev(ids_order))) %>%
  ggplot() + 
    geom_tile(data = start_bleed, aes(y = i, x = t / 2, width = t ), alpha = 0.9, fill = "gray90") +
    geom_tile(data = end_bleed, aes(y = i, x = t + (T - t) / 2, width = (T - t) ), alpha = 0.9, fill = "gray90") +
    geom_point(aes(x = t, y = i), shape = "|", size = 2) + 
    geom_point(data = epi_times_sims, aes(x = t, y = i, color = exp_type),  shape = 15) + theme_bw() + 
    labs(x = "Time in study (days)", y = "Id of individual", color = "Exposure type") + 
    scale_color_manual(values = c("#882255", "#DDCC77")) + 
    geom_hline(yintercept = seq(0.5, N + 0.5, 1), color = "gray90") +
    theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 5)) + 
    ggtitle("Overview of serological survey and underlying epidemic for simulated CES data")


## PLOT B: Plot initial starting titres

p2 <- data.frame(start_model = res$start_titre) %>% 
  pivot_longer(everything(), names_to = "type", values_to = "value") %>%
  ggplot() + geom_density(aes(value), fill = "gray", alpha = 0.2) + 
  labs(x = "Titre at start of study", y = "Density") + theme_bw() +
  ggtitle("Distribution of initial titre values, z")


biomarker_protection <- function(biomarker_quantity, biomarker_prot_midpoint, biomarker_prot_width) {
    risk <- 1 - 1/(1 + exp(biomarker_prot_width * (biomarker_quantity - biomarker_prot_midpoint)))
    return(risk)
}

## PLOT C: Plot initial starting titres

df_cop_show <- data.frame(
    titre = seq(0, 4, 0.1),
    `Model A COP` = (1 - biomarker_protection(seq(0, 4, 0.1), 2, modelcesNoCOP$protection_curve)),
    `Model B COP` = (1 - biomarker_protection(seq(0, 4, 0.1), 2,  modelcesCOP$protection_curve))
) %>% pivot_longer(!titre, names_to = "COPtype", values_to = "prot")


p3 <- df_cop_show %>% ggplot() + 
    geom_line(aes(x = titre, y = prot, color = COPtype), size = 3) + 
    theme_bw() + 
    scale_color_manual(values = c("#888888", "#661100")) +
    labs(x = "Titre at exposure, z", y = "Probability of infection given exposure", color = "Curve type") +
    ggtitle("COP curves for models A + B", expression(f[COP](z)))


## PLOT D: Antibody kinetics model (obtain shematc from somewhere elese)


sim_ab_par = c(modelcesCOP$simpar$a, modelcesCOP$simpar$b, modelcesCOP$simpar$c, modelcesCOP$simpar$biomarker_ceiling_gradient)
names(sim_ab_par) <- c("a", "b", "c", "alpha")

ab_function <- function(a, b, c, alpha, T) {
    0:T %>% map( 
        function(t) {
            if (t < 14) {
                titre_init <- 0 + max(0, (1 - alpha * 0)) * (log(exp(a) +  exp(c)) * (t) / 14);
            } else {
                titre_init <- 0 + max(0, (1 - alpha * 0)) * (log(exp(a) * exp(-b/10 * ((t) - 14)) + exp(c)));
            }
        }
    ) %>% unlist
}

traj_true <- data.frame(
    time = 0:T,
    traj = ab_function(modelcesCOP$simpar$a, modelcesCOP$simpar$b, modelcesCOP$simpar$c, modelcesCOP$simpar$biomarker_ceiling_gradient, T)
)

p4 <- traj_true %>% 
  ggplot() + 
    geom_line(aes(x = time, y = traj), size = 3, alpha = 0.8) + 
    theme_bw() + 
    geom_vline(xintercept = 14, linetype = "dotted") + 
    geom_hline(yintercept = log(exp(sim_ab_par[["a"]]) + exp(sim_ab_par[["c"]])), linetype = "dashed") + 
    geom_hline(yintercept = sim_ab_par[["c"]], linetype = "dashed") + 
   # geom_segment(x = 7, xend = 70, y = 2, yend = 1.15, color = "blue") + 
    geom_segment(aes(x = 7, y = 0, xend = 7, yend = 0 + log(exp(sim_ab_par[["a"]]) + exp(sim_ab_par[["c"]]))),
                  arrow = arrow(length = unit(0.5, "cm"), ends = "both"), size = 2, color = "darkgreen") +
    geom_segment(aes(x = 80, y = sim_ab_par[["c"]], xend = 80, yend = log(exp(sim_ab_par[["a"]]) + exp(sim_ab_par[["c"]]))),
                  arrow = arrow(length = unit(0.5, "cm"), ends = "both"), size = 2, color = "darkgreen") +
    geom_segment(aes(x = 110, y = 0, xend = 110, yend = 0 + sim_ab_par[["c"]]),
                  arrow = arrow(length = unit(0.5, "cm"), ends = "both"), size = 2, color = "darkgreen") +
    geom_label(x = 7, y = 1.2, label = "Peak boost") +
    geom_label(x = 80, y = 1.5, label = "Wane") +
    geom_label(x = 110, y = 0.5, label = "Set point") +

    scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 4.5, 0.5), labels = round(seq(0, 4.5, 0.5), 1)) + 
    labs(x = "Time post-infection, s", y = "Titre increase") + 
    ggtitle("Antibody kinetics post-infection", expression(F[ab](s)))


## PLOT E: Titre dependent boosting

td_true <- data.frame(
    time = seq(0, 4, 0.1),
    traj = pmax(1 - 0.3 * seq(0, 4, 0.1), 0)
)

p5 <- td_true %>% 
  ggplot() + 
    geom_line(aes(x = time, y = traj), size = 3, alpha = 0.8)  + 
    theme_bw() + 
    labs(x = "Log titre at infection", y = "Relative boost in magnitude of kinetics") + 
    ggtitle("Titre-dependent boosting", expression(f[td](z)))


p1 | (p2 / p3 / p4 / p5) + plot_layout(width = c(1, 1)) + plot_annotation(tag_levels = "A")
ggsave(here::here("outputs", "sim_data", "summary_fig_A_CES.pdf"))

traj_post_all <- c(0.1, 0.3, 0.5) %>% map_df(
  function(sigma) {
    a_dist <- rnorm(1000, sim_ab_par[["a"]], sim_ab_par[["a"]] * sigma)
    b_dist <- rnorm(1000, sim_ab_par[["b"]], sim_ab_par[["b"]] * sigma)
    c_dist <- rnorm(1000, sim_ab_par[["c"]], sim_ab_par[["c"]] * sigma)
    alpha_dist <- rnorm(1000, sim_ab_par[["alpha"]], sim_ab_par[["alpha"]] * sigma)

    traj_post <- 1:1000 %>% map_df(
        ~data.frame(
            time = 0:T,
            sigma = sigma,
            value = ab_function(a_dist[.x], b_dist[.x], c_dist[.x], alpha_dist[.x], T)
        )
    )
    traj_post
  }
)

traj_post_all %>% ggplot() + stat_lineribbon(aes(x = time, y = value), 
  alpha = 0.7, fill = "gray", .width = c(0.5, 0.95)) + theme_bw() + 
  labs(x = "Time post-infection, s", y = "Titre increase") + 
  facet_grid(rows = vars(sigma)) + 
  ggtitle("Models with different levels on individual variability in antibody kinetics")
ggsave(here::here("outputs", "sim_data", "summary_fig_B_CES.pdf"))