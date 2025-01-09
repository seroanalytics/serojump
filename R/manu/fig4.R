# WAVE 2: PCR 

devtools::load_all()
model_summary_w2 <- readRDS(file = "./outputs/fits/transvir_data/wave2_base/model_summary.RDS")

save_info_w2 <- list(
    file_name = "transvir_data",
    model_name = "wave2_base"
)


# Rerun these if needed

#
plotMCMCDiagnosis(model_summary_w2, save_info = save_info_w2)
# takes ages to run these so comment out!
plotPostFigs(model_summary_w2, save_info = save_info_w2) 


# Make Figure 4
## Get antibody kinetics

abkin_df <- readRDS(file = "./outputs/fits/transvir_data/wave2_base/figs/post/plt_data/ab_kinetics_recov.RDS")

recode_exposure <- c("delta" = "Infection with Delta variant", "predelta" = "Infection prior to delta vavriant", "vax" = "Vaccination")

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

pB <- exptime_df[[1]] %>% ggplot() +
      geom_histogram(aes(x = time, y = after_stat(count), fill = "gray40"), size = 2, alpha = 0.5, color = "black") + 
    #geom_line(data = bind_rows(hist_data), aes(x = bin, y = count), alpha = 0.1, color = "blue") +
    geom_line(data = data.frame(bin = exptime_df[[2]][[1]]$bin, count = exptime_df[[3]]), aes(x = bin, y = count), color = "#cf725e", size = 1.5) +
    geom_ribbon(data = data.frame(bin = exptime_df[[2]][[1]]$bin, 
                                    lower = exptime_df[[3]] - exptime_df[[4]], 
                                    upper = exptime_df[[3]] + exptime_df[[4]]), 
                aes(x = bin, ymin = lower, ymax = upper, fill = "red"), 
                 alpha = 0.5) +
    scale_fill_manual(
         values =c('gray40'='#62656e','red'='#cf725e'),
         labels = c('PCR-confirmed', 'Inferred total infections')) +
    theme_minimal() +
            labs(x = "Time in study", y = "Number of people infected", fill = "") + 
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


cop_df <- readRDS(file = "./outputs/fits/transvir_data/wave2_base/figs/post/plt_data/cop_data.RDS")

pD <- cop_df[[1]] %>% 
        ggplot() + geom_point(aes(x = titre_val, y = prop, color = biomarker), alpha = 0.6, size = 2) + theme_bw() + 
            geom_smooth(method = "loess", span = 1, aes(x = titre_val, y = prop , color = biomarker), alpha = 0.1, size = 3) +
            ylim(0, 1) +    scale_color_manual(values = c("#003741", "#bf702d")) +
            #facet_wrap(vars(biomarker)) + 
            labs(x = expression("Titre at exposure"), y = expression("Posterior probability of infection"), color = "Biomarker") + 
             theme_minimal() + theme(legend.position = "bottom") + 
    ggtitle("D. Estimated correlate of risk") 



(pA / pBC) | pD + plot_annotation(tag_levels = "A") + plot_layout(width = c(2, 1, 0)) & theme(text = element_text(size = 15, color = "black")) 


ggsave(here::here("outputs", "figs", "fig4.png"), width = 18, height = 10)



# WAVE 2: NO PCR 


devtools::load_all()
model_summary_w2 <- readRDS(file = "./outputs/fits/transvir_data/wave2_no_pcr/model_summary.RDS")

save_info_w2 <- list(
    file_name = "transvir_data",
    model_name = "wave2_no_pcr"
)


# Rerun these if needed

#
plotMCMCDiagnosis(model_summary_w2, save_info = save_info_w2)
# takes ages to run these so comment out!
plotPostFigs(model_summary_w2, save_info = save_info_w2) 


# Make Figure 4
## Get antibody kinetics

abkin_df <- readRDS(file = "./outputs/fits/transvir_data/wave2_no_pcr/figs/post/plt_data/ab_kinetics_recov.RDS")

recode_exposure <- c("delta" = "Infection with Delta variant", "predelta" = "Infection prior to delta vavriant", "vax" = "Vaccination")

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

pB <- exptime_df[[1]] %>% ggplot() +
      geom_histogram(aes(x = time, y = after_stat(count), fill = "gray40"), size = 2, alpha = 0.5, color = "black") + 
    #geom_line(data = bind_rows(hist_data), aes(x = bin, y = count), alpha = 0.1, color = "blue") +
    geom_line(data = data.frame(bin = exptime_df[[2]][[1]]$bin, count = exptime_df[[3]]), aes(x = bin, y = count), color = "#cf725e", size = 1.5) +
    geom_ribbon(data = data.frame(bin = exptime_df[[2]][[1]]$bin, 
                                    lower = exptime_df[[3]] - exptime_df[[4]], 
                                    upper = exptime_df[[3]] + exptime_df[[4]]), 
                aes(x = bin, ymin = lower, ymax = upper, fill = "red"), 
                 alpha = 0.5) +
    scale_fill_manual(
         values =c('gray40'='#62656e','red'='#cf725e'),
         labels = c('PCR-confirmed', 'Inferred total infections')) +
    theme_minimal() +
            labs(x = "Time in study", y = "Number of people infected", fill = "") + 
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


cop_df <- readRDS(file = "./outputs/fits/transvir_data/wave2_no_pcr/figs/post/plt_data/cop_data.RDS")

pD <- cop_df[[1]] %>% 
        ggplot() + geom_point(aes(x = titre_val, y = prop, color = biomarker), alpha = 0.6, size = 2) + theme_bw() + 
            geom_smooth(method = "loess", span = 1, aes(x = titre_val, y = prop , color = biomarker), alpha = 0.1, size = 3) +
            ylim(0, 1) +    scale_color_manual(values = c("#003741", "#bf702d")) +
            #facet_wrap(vars(biomarker)) + 
            labs(x = expression("Titre at exposure"), y = expression("Posterior probability of infection"), color = "Biomarker") + 
             theme_minimal() + theme(legend.position = "bottom") + 
    ggtitle("D. Estimated correlate of risk") 



(pA / pBC) | pD + plot_annotation(tag_levels = "A") + plot_layout(width = c(2, 1, 0)) & theme(text = element_text(size = 15, color = "black")) 
ggsave(here::here("outputs", "figs", "supp", "fig4_no_pcr.png"), width = 18, height = 10)
