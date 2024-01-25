
obs_er <- "0.1"
modelname <- "inferExp"
modelname_sim <-  "cesNoCOP"
outputfull_1 <- readRDS(file = here::here("outputs", "fits", modelname_sim, modelname, paste0("pp_obs_", obs_er, ".RDS")))
modelA_1 <- readRDS(here::here("outputs", "sim_data", modelname_sim, "inputs.RDS"))

obs_er <- "0.1"
modelname <- "inferExp"
modelname_sim <-  "cesCOP"
outputfull_2 <- readRDS(file = here::here("outputs", "fits", modelname_sim, modelname, paste0("pp_obs_", obs_er, ".RDS")))
modelA_2 <- readRDS(here::here("outputs", "sim_data", modelname_sim, "inputs.RDS"))

biomarker_protection <- function(biomarker_quantity, biomarker_prot_midpoint, biomarker_prot_width) {
    risk <- 1 - 1/(1 + exp(biomarker_prot_width * (biomarker_quantity - biomarker_prot_midpoint)))
    return(risk)
}

require(tidybayes)
df_here <- data.frame(
    titre = seq(0, 4, 0.1),
    true_data_A = (1 - biomarker_protection(seq(0, 4, 0.1), 2, modelA_1$protection_curve)),
    true_data_B = (1 - biomarker_protection(seq(0, 4, 0.1), 2, modelA_2$protection_curve))
)

p1 <- outputfull_1$sim_states %>% filter(exp_ind_sim == 1) %>% mutate(inf_state = case_when(inf_ind_sim==0~"Not infected", inf_ind_sim==1~"Infected")) %>%
    ggplot() + geom_point(aes(x = init_titre, y = inf_ind_sim, fill =  inf_state), alpha = 0.8, size = 2, shape = 21) +
        geom_line(data = df_here, aes(x = titre, y = true_data_A, color = "red"), size = 2) + 
        geom_smooth(aes(x = init_titre, y = inf_ind_sim , color = "black")) + theme_bw() + 
        scale_color_manual(name = "Line type", values = c("red" = "red", "black" = "black"), labels = c("Realised simualted curve", "Functional simulated curve")) + 
        labs(x = expression("Pre-infection titre, Z"[j]^0), y = "Probability of infection", fill = "Simulated infection state") + 
        ggtitle("Correlate of protection in simulated data for Model A (No COP)")

p2 <- outputfull_2$sim_states %>% filter(exp_ind_sim == 1)  %>% mutate(inf_state = case_when(inf_ind_sim==0~"Not infected", inf_ind_sim==1~"Infected")) %>%
    ggplot() + geom_point(aes(x = init_titre, y = inf_ind_sim, fill = inf_state), alpha = 0.8, size = 2, shape = 21) +
        geom_line(data = df_here, aes(x = titre, y = true_data_B, color = "red"), size = 2) + 
        geom_smooth(aes(x = init_titre, y = inf_ind_sim, color = "black" )) + theme_bw() + 
        scale_color_manual(name = "Line type", values = c("red" = "red", "black" = "black"), labels = c("Realised simualted curve", "Functional simulated curve")) + 
        labs(x = expression("Pre-infection titre, Z"[j]^0), y = "Probability of infection", fill = "Simulated infection state") + 
        ggtitle("Correlate of protection in simulated data for Model B (Logistic COP)")

p1 / p2 + plot_annotation(tag_levels = "A") + plot_layout(guides = "collect")
ggsave(here::here("outputs", "sim_data", "summary_fig_C_CES.pdf"))