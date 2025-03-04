# WAVE 2: PCR 

devtools::load_all()
require(cmdstanr)
library(pracma)  # for numerical integration

model_summary_w2 <- readRDS(file = "./outputs/fits/transvir_data/wave2_base/model_summary.RDS")

save_info_w2 <- list(
    file_name = "transvir_data",
    model_name = "wave2_base"
)


# Rerun these if needed

#
#plotMCMCDiagnosis(model_summary_w2, save_info = save_info_w2)
# takes ages to run these so comment out!
#plotPostFigs(model_summary_w2, save_info = save_info_w2) 


# Make Figure 4
## Get antibody kinetics

abkin_df <- readRDS(file = "./outputs/fits/transvir_data/wave2_base/figs/post/plt_data/ab_kinetics_recov.RDS")

recode_exposure <- c("delta" = "Infection with Delta variant", "predelta" = "Infection prior to delta variant", "vax" = "Vaccination")

## get AUC 
abkin_df[[1]] %>% group_by(exposure_type, biomarker) %>% summarise(auc = trapz(time, value))

abkin_df[[1]] %>% group_by(exposure_type, biomarker) %>% mutate(above = value > 2) %>% filter(above) %>% as.data.frame



pA <- abkin_df[[1]] %>%  filter(exposure_type != "none") %>% mutate(exposure_type = recode(exposure_type, !!!recode_exposure)) %>%
    ggplot() + 
    geom_ribbon(aes(ymin = .lower, ymax = .upper, x = time, fill = biomarker), size = 3, alpha = 0.3) + 
    geom_line(aes(y = value, x = time, color = biomarker), size = 2) + 
    scale_color_manual(values = c("#003741", "#bf702d")) +
    scale_fill_manual(values = c("#003741", "#bf702d")) +
    theme_minimal() + 
    labs(x = "Time post-infection/since bleed, s", y = "Fold-rise in titre", fill = "Biomarker", color = "Biomarker") + 
    ggtitle("A. Fitted antibody kinetic trajectories") +
    facet_grid(cols = vars(exposure_type)) + 
    theme(legend.position = "bottom")


file_path <- here::here("outputs", "fits", save_info_w2$file_name, save_info_w2$model_name,  "figs", "post")
#plot_exp_times_rec(model_summary_w2, file_path)
#plot_cop_rec(model_summary_w2, file_path)

exptime_df <- readRDS(file = "./outputs/fits/transvir_data/wave2_base/figs/post/plt_data/exposure_time.RDS")

require(lubridate)

exptime_df[[1]] <- exptime_df[[1]] %>% mutate(date = as.Date("2021-03-02") + days(time))
exptime_df[[2]] <- exptime_df[[2]] %>% map(~ .x %>% mutate(date = as.Date("2021-03-02") + days(round(bin, 0))))



pB <- exptime_df[[1]] %>% ggplot() +
      geom_histogram(aes(x = date, y = after_stat(count), fill = "gray40"), size = 2, alpha = 0.5, color = "black") + 
    #geom_line(data = bind_rows(hist_data), aes(x = bin, y = count), alpha = 0.1, color = "blue") +
    geom_line(data = data.frame(bin = exptime_df[[2]][[1]]$date, count = exptime_df[[3]]), aes(x = bin, y = count), color = "#cf725e", size = 1.5) +
    geom_ribbon(data = data.frame(bin = exptime_df[[2]][[1]]$date, 
                                    lower = exptime_df[[3]] - exptime_df[[4]], 
                                    upper = exptime_df[[3]] + exptime_df[[4]]), 
                aes(x = bin, ymin = lower, ymax = upper, fill = "red"), 
                 alpha = 0.5) +
    scale_fill_manual(
         values =c('gray40'='#62656e','red'='#cf725e'),
         labels = c('PCR-confirmed', 'Inferred total infections')) +
    theme_minimal() +
            labs(x = "Time in study (2021)", y = "Number of people infected", fill = "") + 
            ggtitle("B. Inferred epidemic wave") + 
    theme(legend.position = "bottom")


known_delta <- nrow(exptime_df[[1]])


df_plot <- data.frame(map(1:length(exptime_df[[2]]), ~(exptime_df[[2]][[.x]]$count %>% sum)) %>% unlist %>% quantile(., c(0.025, 0.5, 0.975)) ) %>% t
colnames(df_plot) <- c("lb", "mean", "ub")
N <- model_summary_w2$fit$data_t$N




pC <- df_plot %>% 
    ggplot() +
  geom_bar(aes(x = "Delta", y = df_plot[2] / N, fill = "Inferred total infections"), size = 2, alpha = 0.5, stat = "identity", position = position_dodge(width = 0.0), width = 0.5) +
  geom_bar(aes(x = "Delta", y = known_delta / N, fill = "PCR-confirmed"), stat = "identity", position = position_dodge(width = 0.0), width = 0.5, alpha = 1) +
  geom_errorbar( aes(x = "Delta", ymin = lb /N, ymax = ub/ N), width = 0.4, size = 2) + 
  scale_fill_manual(values = c("Inferred total infections" = "#cf725e", "PCR-confirmed" = "#62656e")) +
  labs(title = "C. Inferred attack rates",
       x = "SARS-CoV-2 wave",
       y = "Attack rate", fill = "")  + theme_minimal() + theme(legend.position = "bottom")


pBC <- pB + pC + plot_layout(widths = c(2, 1))

 
### COP work 


model_summary <- model_summary_w2

df_rt_summary <- calculate_reference_titre_expectation(model_summary) %>% 
    group_by(biomarker) %>% mutate(titre_val_s = (titre_val - min(titre_val))/ (max(titre_val) - min(titre_val)))  

biomarkers <- c("spike", "NCP")
df_outputs_cop <- map_df(biomarkers,
    function(b) {
        df_rt_summary_s <- df_rt_summary %>% filter(biomarker == b)

        data_list_A <- list(
            N = nrow(df_rt_summary_s),
            x = df_rt_summary_s$titre_val_s,
            y = df_rt_summary_s$prop
        )

        log_curve_fit <- cmdstan_model( here::here("src", "stan", "cops", "logit_curve_fit.stan") )
        post_A <- log_curve_fit$sample(data = data_list_A, parallel_chains = 4, iter_warmup = 1000, iter_sampling = 1000)

        y_hat_new_post_full <- post_A$draws() %>% spread_draws(y_hat_rel[i], y_prot[i], y_hat_new[i], x_new[i]) %>% 
            mutate(biomarker = b)
    }
)

df_rt_summary_spike <- df_rt_summary %>% filter(biomarker == "spike")
labels_x_spike <- seq(min(df_rt_summary_spike$titre_val), max(df_rt_summary_spike$titre_val), length.out = 10) %>% round(2)

p1A <- df_outputs_cop %>% filter(biomarker == "spike") %>%
    ggplot() +
    stat_lineribbon(aes(x = x_new, y = y_hat_new, group = biomarker, color = biomarker, fill = biomarker), .width = 0.95,alpha = 0.3) + 
    geom_point(data = df_rt_summary %>% filter(biomarker == "spike"), aes(x = titre_val_s, y = prop, color = biomarker), alpha = 0.6, size = 2) +
    scale_color_manual(values = c("#003741", "#bf702d")) +
    scale_fill_manual(values = c("#003741", "#bf702d")) +
    facet_wrap(vars(biomarker)) + theme_minimal() +
    guides(color = "none", fill = "none") +
    scale_x_continuous(breaks = seq(0, 1, length.out = 10), labels = labels_x_spike) + 
    labs(x = "Log titre", y = "Posterior probability of infection", color = "Biomarker", fill = "Biomarker") + 
    ggtitle("D. Fitted curves for infection risk")

df_rt_summary_spike <- df_rt_summary %>% filter(biomarker == "NCP")
labels_x_ncp <- seq(min(df_rt_summary_spike$titre_val), max(df_rt_summary_spike$titre_val), length.out = 10) %>% round(2)

p1B <- df_outputs_cop %>% filter(biomarker == "NCP") %>%
    ggplot() +
    stat_lineribbon(aes(x = x_new, y = y_hat_new, group = biomarker, color = biomarker, fill = biomarker), .width = 0.95,alpha = 0.3) + 
    geom_point(data = df_rt_summary %>% filter(biomarker == "NCP"), aes(x = titre_val_s, y = prop, color = biomarker), alpha = 0.6, size = 2) +
    scale_color_manual(values = c("#bf702d")) +
    scale_fill_manual(values = c("#bf702d")) +
    facet_wrap(vars(biomarker)) + theme_minimal() +
    guides(color = "none", fill = "none") +
    scale_x_continuous(breaks = seq(0, 1, length.out = 10), labels = labels_x_ncp) + 
    labs(x = "Log titre", y = "Posterior probability of infection", color = "Biomarker", fill = "Biomarker")

p1 <- p1A + p1B

p2 <- df_outputs_cop %>% 
    ggplot() + stat_lineribbon(aes(x = x_new, y = y_prot, group = biomarker, color = biomarker, fill = biomarker), .width = 0.01, alpha = 1, size = 2) + theme_bw() + 
    scale_color_manual(values = c("#003741", "#bf702d")) + scale_fill_manual(values = c("#003741", "#bf702d"))  + ylim(0, 1) + 
    labs(x = "Normalised titre", y = "Protection probability", color = "Biomarker", fill = "Biomarker") + 
    scale_x_continuous(
        name = "log10 Titre value (spike)", 
        breaks = seq(0, 1, length.out = 10), labels = labels_x_spike,
        sec.axis = sec_axis(~ . * 1, 
        breaks = seq(0, 1, length.out = 10), labels = labels_x_ncp, name = "log10 Titre value (NCP)")
    ) + theme_minimal() + 
    ggtitle("E. Fitted curves for absolute COP")

p3 <- df_outputs_cop %>% 
    ggplot() + stat_lineribbon(aes(x = x_new, y = y_hat_rel, group = biomarker, color = biomarker, fill = biomarker), .width = 0.01,alpha = 1, size = 2) + theme_bw()+ 
    scale_color_manual(values = c("#003741", "#bf702d")) + scale_fill_manual(values = c("#003741", "#bf702d"))  + ylim(0, 1) + 
    labs(x = "Normalised titre", y = "Relative risk of infection", color = "Biomarker", fill = "Biomarker") + 
    scale_x_continuous(
        name = "log10 Titre value (spike)", 
        breaks = seq(0, 1, length.out = 10), labels = labels_x_spike,
        sec.axis = sec_axis(~ . * 1, 
        breaks = seq(0, 1, length.out = 10), labels = labels_x_ncp, name = "log10 Titre value (NCP)")
    ) + theme_minimal() + 
    ggtitle("F. Fitted curves for relative COP")


pCOP <- p1 / (p2 + p3) + plot_layout(guides = "collect")


(pA / pBC) / pCOP + plot_layout(height = c(1, 1, 3)) & theme(text = element_text(size = 12, color = "black")) 
ggsave(here::here("outputs", "figs", "fig4.png"), width = 12, height = 15)













####################################
######### WAVE 2: NO PCR  #########
####################################

devtools::load_all()
model_summary_w2 <- readRDS(file = "./outputs/fits/transvir_data/wave2_no_pcr/model_summary.RDS")

save_info_w2 <- list(
   file_name = "transvir_data",
    model_name = "wave2_no_pcr"
)


# Rerun these if needed

#
#plotMCMCDiagnosis(model_summary_w2, save_info = save_info_w2)
# takes ages to run these so comment out!
#plotPostFigs(model_summary_w2, save_info = save_info_w2) 


# Make Figure 4
## Get antibody kinetics

abkin_df <- readRDS(file = "./outputs/fits/transvir_data/wave2_no_pcr/figs/post/plt_data/ab_kinetics_recov.RDS")

recode_exposure <- c("delta" = "Infection with Delta variant", "predelta" = "Infection prior to delta variant", "vax" = "Vaccination")

pA <- abkin_df[[1]] %>%  filter(exposure_type != "none") %>% mutate(exposure_type = recode(exposure_type, !!!recode_exposure)) %>%
    ggplot() + 
    geom_ribbon(aes(ymin = .lower, ymax = .upper, x = time, fill = biomarker), size = 3, alpha = 0.3) + 
    geom_line(aes(y = value, x = time, color = biomarker), size = 2) + 
    scale_color_manual(values = c("#003741", "#bf702d")) +
    scale_fill_manual(values = c("#003741", "#bf702d")) +
    theme_minimal() + 
    labs(x = "Time post-infection/since bleed, s", y = "Fold-rise in titre", fill = "Biomarker", color = "Biomarker") + 
    ggtitle("A. Fitted antibody kinetic trajectories") +
    facet_grid(cols = vars(exposure_type)) + 
    theme(legend.position = "bottom")


file_path <- here::here("outputs", "fits", save_info_w2$file_name, save_info_w2$model_name,  "figs", "post")
#plot_exp_times_rec(model_summary_w2, file_path)
#plot_cop_rec(model_summary_w2, file_path)

exptime_df <- readRDS(file = "./outputs/fits/transvir_data/wave2_no_pcr/figs/post/plt_data/exposure_time.RDS")

#exptime_df[[1]] <- exptime_df[[1]] %>% mutate(date = as.Date("2021-03-02") + days(time))
exptime_df[[2]] <- exptime_df[[2]] %>% map(~ .x %>% mutate(date = as.Date("2021-03-02") + days(round(bin, 0))))


pB <- exptime_df[[1]] %>% ggplot() +
      geom_histogram(aes(x = time, y = after_stat(count), fill = "gray40"), size = 2, alpha = 0.5, color = "black") + 
    #geom_line(data = bind_rows(hist_data), aes(x = bin, y = count), alpha = 0.1, color = "blue") +
    geom_line(data = data.frame(bin = exptime_df[[2]][[1]]$date, count = exptime_df[[3]]), aes(x = bin, y = count), color = "#cf725e", size = 1.5) +
    geom_ribbon(data = data.frame(bin = exptime_df[[2]][[1]]$date, 
                                    lower = exptime_df[[3]] - exptime_df[[4]], 
                                    upper = exptime_df[[3]] + exptime_df[[4]]), 
                aes(x = bin, ymin = lower, ymax = upper, fill = "red"), 
                 alpha = 0.5) +
    scale_fill_manual(
         values =c('gray40'='#62656e','red'='#cf725e'),
         labels = c('PCR-confirmed', 'Inferred total infections')) +
    theme_minimal() +
            labs(x = "Time in study (2021)", y = "Number of people infected", fill = "") + 
            ggtitle("B. Inferred epidemic wave") + 
    theme(legend.position = "bottom")


known_delta <- nrow(exptime_df[[1]])


df_plot <- data.frame(map(1:length(exptime_df[[2]]), ~(exptime_df[[2]][[.x]]$count %>% sum)) %>% unlist %>% quantile(., c(0.025, 0.5, 0.975)) ) %>% t
colnames(df_plot) <- c("lb", "mean", "ub")
N <- model_summary_w2$fit$data_t$N
pC <- df_plot %>% 
    ggplot() +
  geom_bar(aes(x = "Delta", y = df_plot[2] / N, fill = "Inferred total infections"), size = 2, alpha = 0.5, stat = "identity", position = position_dodge(width = 0.0), width = 0.5) +
  geom_bar(aes(x = "Delta", y = known_delta / N, fill = "PCR-confirmed"), stat = "identity", position = position_dodge(width = 0.0), width = 0.5, alpha = 1) +
  geom_errorbar( aes(x = "Delta", ymin = lb /N, ymax = ub/ N), width = 0.4, size = 2) + 
  scale_fill_manual(values = c("Inferred total infections" = "#cf725e", "PCR-confirmed" = "#62656e")) +
  labs(title = "C. Inferred attack rates",
       x = "SARS-CoV-2 wave",
       y = "Attack rate", fill = "")  + theme_minimal() + theme(legend.position = "bottom")


pBC <- pB + pC + plot_layout(widths = c(2, 1))

 
### COP work 


model_summary <- model_summary_w2

df_rt_summary <- calculate_reference_titre_expectation(model_summary) %>% 
    group_by(biomarker) %>% mutate(titre_val_s = (titre_val - min(titre_val))/ (max(titre_val) - min(titre_val)))  

biomarkers <- c("spike", "NCP")
df_outputs_cop <- map_df(biomarkers,
    function(b) {
        df_rt_summary_s <- df_rt_summary %>% filter(biomarker == b)

        data_list_A <- list(
            N = nrow(df_rt_summary_s),
            x = df_rt_summary_s$titre_val_s,
            y = df_rt_summary_s$prop
        )

        log_curve_fit <- cmdstan_model( here::here("src", "stan", "cops", "logit_curve_fit.stan") )
        post_A <- log_curve_fit$sample(data = data_list_A, parallel_chains = 4, iter_warmup = 1000, iter_sampling = 1000)

        y_hat_new_post_full <- post_A$draws() %>% spread_draws(y_hat_rel[i], y_prot[i], y_hat_new[i], x_new[i]) %>% 
            mutate(biomarker = b)
    }
)

df_rt_summary_spike <- df_rt_summary %>% filter(biomarker == "spike")
labels_x_spike <- seq(min(df_rt_summary_spike$titre_val), max(df_rt_summary_spike$titre_val), length.out = 10) %>% round(2)

p1A <- df_outputs_cop %>% filter(biomarker == "spike") %>%
    ggplot() +
    stat_lineribbon(aes(x = x_new, y = y_hat_new, group = biomarker, color = biomarker, fill = biomarker), .width = 0.95,alpha = 0.3) + 
    geom_point(data = df_rt_summary %>% filter(biomarker == "spike"), aes(x = titre_val_s, y = prop, color = biomarker), alpha = 0.6, size = 2) +
    scale_color_manual(values = c("#003741", "#bf702d")) +
    scale_fill_manual(values = c("#003741", "#bf702d")) +
    facet_wrap(vars(biomarker)) + theme_minimal() +
    guides(color = "none", fill = "none") +
    scale_x_continuous(breaks = seq(0, 1, length.out = 10), labels = labels_x_spike) + 
    labs(x = "Log titre", y = "Posterior probability of infection", color = "Biomarker", fill = "Biomarker") + 
    ggtitle("D. Fitted curves for infection risk")

df_rt_summary_spike <- df_rt_summary %>% filter(biomarker == "NCP")
labels_x_ncp <- seq(min(df_rt_summary_spike$titre_val), max(df_rt_summary_spike$titre_val), length.out = 10) %>% round(2)

p1B <- df_outputs_cop %>% filter(biomarker == "NCP") %>%
    ggplot() +
    stat_lineribbon(aes(x = x_new, y = y_hat_new, group = biomarker, color = biomarker, fill = biomarker), .width = 0.95,alpha = 0.3) + 
    geom_point(data = df_rt_summary %>% filter(biomarker == "NCP"), aes(x = titre_val_s, y = prop, color = biomarker), alpha = 0.6, size = 2) +
    scale_color_manual(values = c("#bf702d")) +
    scale_fill_manual(values = c("#bf702d")) +
    facet_wrap(vars(biomarker)) + theme_minimal() +
    guides(color = "none", fill = "none") +
    scale_x_continuous(breaks = seq(0, 1, length.out = 10), labels = labels_x_ncp) + 
    labs(x = "Log titre", y = "Posterior probability of infection", color = "Biomarker", fill = "Biomarker")

p1 <- p1A + p1B

p2 <- df_outputs_cop %>% 
    ggplot() + stat_lineribbon(aes(x = x_new, y = y_prot, group = biomarker, color = biomarker, fill = biomarker), .width = 0.01, alpha = 1, size = 2) + theme_bw() + 
    scale_color_manual(values = c("#003741", "#bf702d")) + scale_fill_manual(values = c("#003741", "#bf702d"))  + ylim(0, 1) + 
    labs(x = "Normalised titre", y = "Protection probability", color = "Biomarker", fill = "Biomarker") + 
    scale_x_continuous(
        name = "log10 Titre value (spike)", 
        breaks = seq(0, 1, length.out = 10), labels = labels_x_spike,
        sec.axis = sec_axis(~ . * 1, 
        breaks = seq(0, 1, length.out = 10), labels = labels_x_ncp, name = "log10 Titre value (NCP)")
    ) + theme_minimal() + 
    ggtitle("E. Fitted curves for absolute COP")

p3 <- df_outputs_cop %>% 
    ggplot() + stat_lineribbon(aes(x = x_new, y = y_hat_rel, group = biomarker, color = biomarker, fill = biomarker), .width = 0.01,alpha = 1, size = 2) + theme_bw()+ 
    scale_color_manual(values = c("#003741", "#bf702d")) + scale_fill_manual(values = c("#003741", "#bf702d"))  + ylim(0, 1) + 
    labs(x = "Normalised titre", y = "Relative risk of infection", color = "Biomarker", fill = "Biomarker") + 
    scale_x_continuous(
        name = "log10 Titre value (spike)", 
        breaks = seq(0, 1, length.out = 10), labels = labels_x_spike,
        sec.axis = sec_axis(~ . * 1, 
        breaks = seq(0, 1, length.out = 10), labels = labels_x_ncp, name = "log10 Titre value (NCP)")
    ) + theme_minimal() + 
    ggtitle("F. Fitted curves for relative COP")


pCOP <- p1 / (p2 + p3) + plot_layout(guides = "collect")

(pA / pBC) / pCOP + plot_layout(height = c(1, 1, 3)) & theme(text = element_text(size = 12, color = "black")) 
ggsave(here::here("outputs", "figs", "supp", "fig4_no_pcr.png"), width = 12, height = 15)













## HOLD ON LEMME COOK ON THIS:

model_summary <- model_summary_w2

fitfull <- model_summary$fit    
outputfull <- model_summary$post

fit_states <- outputfull$fit_states
model_outline <- fitfull$model

df_full_info <- fit_states %>% mutate(titre_type = case_when(
    inf_time == -1 ~ "No inf",
    inf_time != -1 ~ "Inf")
)

P <- 10

infection_risk <- map_df(1:length(model_outline$observationalModel), 
    function(i) {
    biomarker <- model_outline$observationalModel[[i]]$biomarker

    df_full_info %>% mutate(titre_cut = cut(!!sym(biomarker), breaks = P)) %>% 
        group_by(titre_cut) %>% summarise(inf = sum(inf_ind), n = n()) %>% mutate(inf_risk = inf/n) %>% 
        complete(titre_cut, fill = list(inf = 0, n = 0, inf_risk = 0)) %>%
        mutate(biomarker = biomarker, decile = 1:P)  %>% mutate(dens = n / sum(n))
        }
    )
    


# FT-Style Theme
theme_ft <- function() {
  theme_minimal(base_family = "Arial") +  # Simple sans-serif font
    theme(
      panel.grid.major.y = element_line(color = "grey80", size = 0.3),  # Light grid lines
      panel.grid.major.x = element_blank(),  # No vertical grid lines
      panel.grid.minor = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 14, color = "grey40", margin = margin(b = 15)),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      plot.margin = margin(20, 20, 20, 20)
    )
}


# This is a total plot of the Infection risk per decile of titre 
p1 <- infection_risk %>% 
    ggplot() + 
    geom_col(aes(x = decile, y = dens, fill = biomarker), position = "dodge", size = 2, alpha = 0.2) +
    geom_point(aes(x = decile, y = inf_risk, fill = biomarker, size = n), color = "black", shape = 21) +
    geom_line(aes(x = decile, y = inf_risk, color = biomarker), size = 0.5, alpha = 0.5) + theme_minimal() + 
    scale_x_continuous(breaks = 1:P) + 
    guides(size = "none") + 
    scale_size(range = c(0, 10 )) +
    scale_color_manual(values = c("#003741", "#bf702d")) +
    scale_fill_manual(values = c("#003741", "#bf702d")) +
    labs(x = "Decile of titre", y = "Infection risk (points)\nProportion of density (bars)", color = "Biomarker", fill = "Biomarker") + 
    theme_ft() + theme(text = element_text(size = 15)) 


# Plot the relative risks of infection
infection_risk %>% group_by(biomarker) %>% mutate(r = row_number()) %>% 
    mutate(rr = inf_risk / inf_risk[1]) %>% 
    ggplot() + 
        geom_point(aes(x = decile, y = rr, fill = biomarker, size = n), color = "black", shape = 21) +
        geom_line(aes(x = decile, y = rr, color = biomarker), size = 0.5, alpha = 0.5) + theme_minimal() + 
        guides(size = "none") + 
        scale_size(range = c(0, 10 )) +
        scale_color_manual(values = c("#003741", "#bf702d")) +
        scale_fill_manual(values = c("#003741", "#bf702d")) +
        labs(x = "Decile of titre", y = "Relative risk of infection\n(relative to lowest tire)", color = "Biomarker", fill = "Biomarker") + 
        theme_ft() + theme(text = element_text(size = 15)) + 
        scale_x_continuous(breaks = 1:10) + scale_y_continuous(breaks = seq(0, 1, 0.2))

# I would now look into fitting the relative risk 
library(dplyr)
library(stringr)


infection_risk <- map_df(1:length(model_outline$observationalModel), 
    function(i) {
    biomarker <- model_outline$observationalModel[[i]]$biomarker

    df_full_info %>% mutate(titre_cut = cut(!!sym(biomarker), breaks = P)) %>% 
        group_by(titre_cut) %>% summarise(inf = sum(inf_ind), n = n(), .groups = "drop") %>% mutate(inf_risk = inf/n) %>% 
        complete(titre_cut, fill = list(inf = 0, n = 0, inf_risk = 0)) %>%
        mutate(biomarker = biomarker, decile = 1:P)  %>% mutate(dens = n / sum(n)) %>%
        mutate(    titre_cut = stringr::str_remove_all(titre_cut, "[\\(\\)\\[\\]]"),  # Remove brackets
                lower = as.numeric(stringr::str_extract(titre_cut, "^[^,]+")),     # Extract lower bound
                upper = as.numeric(stringr::str_extract(titre_cut, "[^,]+$")),     # Extract upper bound
                midpoint = (lower + upper) / 2                    # Compute midpoint                
        )
        }
    )


biomarker_i <- model_outline$observationalModel[[1]]$biomarker
infection_risk_bio <- infection_risk %>% filter(biomarker == biomarker_i)
 infection_risk_bio_n <- infection_risk_bio 


## use stan model
require(cmdstanr)
data_list <- list(
    N = nrow(infection_risk_bio_n),
    quantiles = infection_risk_bio_n$decile,
    adj_total = infection_risk_bio_n$n,
    count_1 = infection_risk_bio_n$inf
)

weighted_log_reg <- cmdstan_model( here::here("src", "stan", "weighted_log_reg.stan") )
post <- weighted_log_reg$sample(data = data_list, parallel_chains = 4, cores = 4, iter_warmup = 1000, iter_sampling = 1000)

post$draws() %>% gather_draws(p[i], k, x_0) %>% mutate(p_E_j = p_E[j]) %>% mutate(biomarker = biomarker_i)



# Now look into calculating the COP for a given calculation 
# Fit something under tha assumption everyone is exposure and thus
# P(Z = 1 || E = 1) = logistic regression 
# Example

library(dplyr)
library(stringr)



infection_risk <- map_df(1:length(model_outline$observationalModel), 
    function(i) {
    biomarker <- model_outline$observationalModel[[i]]$biomarker

    df_full_info %>% mutate(titre_cut = cut(!!sym(biomarker), breaks = P)) %>% 
        group_by(titre_cut) %>% summarise(inf = sum(inf_ind), n = n(), .groups = "drop") %>% mutate(inf_risk = inf/n) %>% 
        complete(titre_cut, fill = list(inf = 0, n = 0, inf_risk = 0)) %>%
        mutate(biomarker = biomarker, decile = 1:P)  %>% mutate(dens = n / sum(n)) %>%
        mutate(    titre_cut = stringr::str_remove_all(titre_cut, "[\\(\\)\\[\\]]"),  # Remove brackets
                lower = as.numeric(stringr::str_extract(titre_cut, "^[^,]+")),     # Extract lower bound
                upper = as.numeric(stringr::str_extract(titre_cut, "[^,]+$")),     # Extract upper bound
                midpoint = (lower + upper) / 2                    # Compute midpoint                
        )
        }
    )
    
# Example Data
set.seed(42)

#df_full_info %>% pull()

# Extract breakpoints
#breaks <- as.numeric(levels(infection_risk$titre_cut))  # Get breakpoints
#midpoints <- (head(breaks, -1) + tail(breaks, -1)) / 2  # Compute midpoints


#infection_risk_bio <- infection_risk %>% filter(biomarker == "NCP")
#glm_model <- glm(inf_risk ~ decile, weights = n, data = infection_risk %>% filter(biomarker == "NCP"), family = binomial)

# Generate a sequence of predictor values for plotting
#s_seq <- data.frame(decile = seq(1, 20, 1))

# Get predicted probabilities
#infection_risk_bio$predicted_prob <- predict(glm_model, newdata = s_seq, type = "response")


#ggplot(df_agg, aes(x = s_star, y = prop_1)) +
#  geom_point(alpha = 0.6) +
#  geom_smooth(method = "glm", method.args = list(family = binomial), se = FALSE) +
#  labs(title = "Logistic Regression on Aggregated Data", y = "Infection Probability", x = "s_star") +
#  theme_minimal()

p_E <- seq(0.0, 1, 0.1)

full_exp_explore <- map_df(1:length(model_outline$observationalModel), 
    function(b) {
        biomarker_i <- model_outline$observationalModel[[b]]$biomarker
        infection_risk_bio <- infection_risk %>% filter(biomarker == biomarker_i)
        map_df(1:11, 
            function(j) {
                infection_risk_bio_n <- infection_risk_bio %>% mutate(adj_n = round(inf + p_E[j] * (n - inf), 0))


                ## use stan model
                require(cmdstanr)
                data_list <- list(
                    N = nrow(infection_risk_bio_n),
                    quantiles = infection_risk_bio_n$decile,
                    adj_total = infection_risk_bio_n$adj_n,
                    count_1 = infection_risk_bio_n$inf
                )

                weighted_log_reg <- cmdstan_model( here::here("src", "stan", "weighted_log_reg.stan") )

                post <- weighted_log_reg$sample(data = data_list, parallel_chains = 4, cores = 4, iter_warmup = 1000, iter_sampling = 1000)

                post$draws() %>% gather_draws(p[i], k, x_0) %>% mutate(p_E_j = p_E[j]) %>% mutate(biomarker = biomarker_i)

            }

        )
    }
 )

full_exp_explore_p2 <- full_exp_explore %>% left_join(infection_risk %>% select(i = decile, midpoint, biomarker))
p2A <- full_exp_explore_p2 %>% filter(biomarker == "spike") %>%
    ggplot() + 
        stat_summary(aes(x = i, .value, color = as.character(p_E_j)), geom = "line") + 
        facet_wrap(vars(biomarker) ) + theme_bw() + 
        scale_x_continuous(breaks = c(1:P)[seq(1, P, 2)], labels = infection_risk %>% filter(biomarker == "spike") %>% pull(midpoint) %>% .[seq(1, P, 2)]) +  
        labs(x = "Titre (log2)", y = "P(I | E = 1)", color = "NIER", fill = "") + theme_ft() + theme(text = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust = 1)) 

p2B <- full_exp_explore_p2 %>% filter(biomarker == "NCP") %>%
    ggplot() + 
        stat_summary(aes(x = i, .value, color = as.character(p_E_j)), geom = "line") + 
        facet_wrap(vars(biomarker) ) + theme_bw() + 
        scale_x_continuous(breaks = c(1:P)[seq(1, P, 2)], labels = infection_risk %>% filter(biomarker == "NCP") %>% pull(midpoint) %>% .[seq(1, P, 2)]) + 
        labs(x = "Titre (log2)", y = "P(I | E = 1)", color = "NEIR", fill = "") + theme_ft() + theme(text = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust = 1)) 

p2 <- p2A + p2B + plot_layout(guides = "collect")

p3A <- full_exp_explore %>% filter(.variable == "k") %>% filter(p_E_j > 0) %>% 
    mutate() %>% 
    ggplot() + 
    stat_summary(aes(x = as.character(p_E_j), y = .value, color = biomarker, group = biomarker ), size = 1) +
    stat_summary(aes(x = as.character(p_E_j), y = .value, color = biomarker, group = biomarker ), geom = "line", size =1 ) +
    scale_color_manual(values = c("#003741", "#bf702d")) +
    scale_fill_manual(values = c("#003741", "#bf702d")) +
    labs(x = "Proportion of non-infected individuals who are exposed", 
        y = "Graident of logistic curve", title = "Posterior distribution of parameters across different NIER") + 
    theme_bw() + theme_ft() + theme(text = element_text(size = 15)) 

p3B <- full_exp_explore %>% filter(.variable == "x_0") %>% filter(p_E_j > 0) %>% 
    mutate() %>% 
    ggplot() + 
    stat_summary(aes(x = as.character(p_E_j), y = .value, color = biomarker, group = biomarker ), size = 1) +
    stat_summary(aes(x = as.character(p_E_j), y = .value, color = biomarker, group = biomarker ), geom = "line", size =1 ) +
    scale_color_manual(values = c("#003741", "#bf702d")) +
    scale_fill_manual(values = c("#003741", "#bf702d")) +
    labs(x = "Proportion of non-infected individuals who are exposed", 
        y = "Graident of logistic curve", title = "Posterior distribution of parameters across different NIER", color = "Biomarker") + 
    theme_bw() + theme(text = element_text(size = 15)) 

p1 / (p2) / p3A + plot_annotation(tag_levels = "A")  




#cop_df <- readRDS(file = "./outputs/fits/transvir_data/wave2_no_pcr/figs/post/plt_data/cop_data.RDS")
#pD <- cop_df[[1]] %>% 
#        ggplot() + geom_point(aes(x = titre_val, y = prop, color = biomarker), alpha = 0.6, size = 2) + theme_bw() + 
#            geom_smooth(method = "loess", span = 1, aes(x = titre_val, y = prop , color = biomarker), alpha = 0.1, size = 3) +
#            ylim(0, 1) +    scale_color_manual(values = c("#003741", "#bf702d")) +
            #facet_wrap(vars(biomarker)) + 
#            labs(x = expression("Titre at exposure"), y = expression("Posterior probability of infection"), color = "Biomarker") + 
#             theme_minimal() + theme(legend.position = "bottom") + 
#    ggtitle("D. Estimated correlate of risk") 



#(pA / pBC) | pD + plot_annotation(tag_levels = "A") + plot_layout(width = c(2, 1, 0)) & theme(text = element_text(size = 15, color = "black")) 
#ggsave(here::here("outputs", "figs", "supp", "fig4_no_pcr.png"), width = 18, height = 10)
