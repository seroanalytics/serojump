# Detection recovery of `serojump`
devtools::load_all()

wave_full <- readRDS(file = "./vignettes/wave2_export.RDS")
sero_data_w2 <- wave_full$sero_data

pid_linked <- sero_data_w2 %>% select(pid) %>% unique %>% mutate(id = 1:nrow(.))
known_delta <- wave_full$known_inf %>% filter(exposure_type == "delta") 


## Extract the fourfold rise sensitivtiy 

inferred_fourfold_s <- sero_data_w2 %>% group_by(id) %>% mutate(inferred_4_fold = (spike - lag(spike)) > log10(4)) %>%
    filter(!is.na(inferred_4_fold))
recover_fourfold_s <- inferred_fourfold_s %>%
    left_join(known_delta %>% select(!time)) %>% filter(!is.na(exposure_type)) %>% mutate(inferred_4_fold = as.numeric(inferred_4_fold)) %>%
    mutate(true_pos = inferred_4_fold > 0.5) 

inferred_fourfold_n <- sero_data_w2 %>% group_by(id) %>% mutate(inferred_4_fold = (NCP - lag(NCP)) > log10(4)) %>%
    filter(!is.na(inferred_4_fold))
recover_fourfold_n <- inferred_fourfold_n %>%
    left_join(known_delta %>% select(!time)) %>% filter(!is.na(exposure_type)) %>% mutate(inferred_4_fold = as.numeric(inferred_4_fold)) %>%
    mutate(true_pos = inferred_4_fold > 0.5) 


## Extract the change is seroneg -> seropos sensitivity

recover_seropos_s <- sero_data_w2 %>% group_by(id) %>% mutate(inferred_pos = spike > log10(7)) %>%
    mutate(inferred_pos_prev = lag(inferred_pos)) %>% mutate(sero_pos = (inferred_pos == 1) & (inferred_pos_prev == 0)) %>%
    filter(!is.na(inferred_pos_prev)) %>% select(!c( inferred_pos, inferred_pos_prev )) %>%
    left_join(known_delta %>% select(!time)) %>% filter(!is.na(exposure_type)) %>% mutate(sero_pos = as.numeric(sero_pos)) %>%
    mutate(true_pos = sero_pos > 0.5) 

recover_seropos_n <- sero_data_w2 %>% group_by(id) %>% mutate(inferred_pos = NCP > log10(2.920)) %>%
    mutate(inferred_pos_prev = lag(inferred_pos)) %>% mutate(sero_pos = (inferred_pos == 1) & (inferred_pos_prev == 0)) %>%
    filter(!is.na(inferred_pos_prev)) %>% select(!c( inferred_pos, inferred_pos_prev )) %>%
    left_join(known_delta %>% select(!time)) %>% filter(!is.na(exposure_type)) %>% mutate(sero_pos = as.numeric(sero_pos)) %>%
    mutate(true_pos = sero_pos > 0.5) 


model_summary <- readRDS(file = "./outputs/fits/transvir_data/wave2_no_pcr/model_summary.RDS")

# model predict sensitivity of serojump
S <- (model_summary$post$fit_states$sample_c %>% unique %>% length) * (model_summary$post$fit_states$chain_no %>% unique %>% length)
N <- model_summary$fit$data$N
recover_serojump <- model_summary$post$fit_states %>% group_by(id) %>% filter(inf_time > -1) %>% summarise(n = n() / S) %>% 
    complete(id = 1:N, fill = list(n = 0)) %>% left_join(known_delta, by = join_by(id)) %>% as.data.frame() %>% filter(!is.na(pid)) %>% 
    mutate(true_pos = n > 0.5) 

sens_4_s <- mean(recover_fourfold_s$true_pos)
p1 <- recover_fourfold_s %>%
    ggplot() + geom_point(aes(x = id, inferred_4_fold, fill = true_pos), size = 5, shape = 21, alpha = 0.5) + 
    theme_bw() + theme(text = element_text(size = 20)) + 
    geom_hline(yintercept = sens_4_s, size = 2, color = "#E8B88D") + 
    geom_text(x = 100, y = sens_4_s + 0.05, label = paste0("Sensitivity: ", round(sens_4_s, 2)), size= 8) + 
    scale_fill_manual(values = c( "#E74C3C", "#6DBE45")) + 
    labs(x = "Id of individual", y = "Inferred infection status", title = "Four-fold rise (spike)", 
        fill = "Inferred positive") + theme(legend.position = "top")


sens_4_n <- mean(recover_fourfold_n$true_pos)
p2 <- recover_fourfold_n %>%
    ggplot() + geom_point(aes(x = id, inferred_4_fold, fill = true_pos), size = 5, shape = 21, alpha = 0.5) + 
    theme_bw() + theme(text = element_text(size = 20)) + 
    geom_hline(yintercept = sens_4_n, size = 2, color = "#E8B88D") + 
    geom_text(x = 100, y = sens_4_n + 0.05, label = paste0("Sensitivity: ", round(sens_4_n, 2)), size= 8) + 
    scale_fill_manual(values = c( "#E74C3C", "#6DBE45")) + 
    labs(x = "Id of individual", y = "Inferred infection status", title = "Four-fold rise (NCP)", 
        fill = "Inferred positive") + theme(legend.position = "top")

sens_p_s <- mean(recover_seropos_s$true_pos)

p3 <- recover_seropos_s %>%
    ggplot() + geom_point(aes(x = id, sero_pos, fill = true_pos), size = 5, shape = 21, alpha = 0.5) + 
    theme_bw() + theme(text = element_text(size = 20)) + 
    geom_hline(yintercept = sens_p_s, size = 2, color = "#E8B88D") + 
    scale_fill_manual(values = c( "#E74C3C", "#6DBE45")) + 
    geom_text(x = 100, y = sens_p_s + 0.05, label = paste0("Sensitivity: ", round(sens_p_s, 2)), size= 8) + 
    labs(x = "Id of individual", y = "Inferred infection status", title = "Seropos threshold (spike)", 
        fill = "Inferred positive") + theme(legend.position = "top")

sens_p_n <- mean(recover_seropos_n$true_pos)

p4 <- recover_seropos_n %>%
    ggplot() + geom_point(aes(x = id, sero_pos, fill = true_pos), size = 5, shape = 21, alpha = 0.5) + 
    theme_bw() + theme(text = element_text(size = 20)) + 
    geom_hline(yintercept = sens_p_n, size = 2, color = "#E8B88D") + 
    scale_fill_manual(values = c( "#E74C3C", "#6DBE45")) + 
    geom_text(x = 100, y = sens_p_n + 0.05, label = paste0("Sensitivity: ", round(sens_p_n, 2)), size= 8) + 
    labs(x = "Id of individual", y = "Inferred infection status", title = "Seropos threshold (NCP)", 
        fill = "Inferred positive") + theme(legend.position = "top")


sens_sj <- mean(recover_serojump$true_pos)

p5 <- recover_serojump %>%
    ggplot() + geom_point(aes(x = id, n, fill = true_pos), size = 5, shape = 21, alpha = 0.5) + 
    theme_bw() + theme(text = element_text(size = 20)) + 
    geom_hline(yintercept = sens_sj, size = 2, color = "#E8B88D") + 
    scale_fill_manual(values = c( "#E74C3C", "#6DBE45")) + 
    geom_text(x = 100, y = sens_sj + 0.05, label = paste0("Sensitivity: ", round(sens_sj, 2)), size= 8) + 
    labs(x = "Id of individual", y = "Posterior probability of infection", title = "`serojump`", 
        fill = "Inferred positive")  + theme(legend.position = "top")

p0 <- (p1 + p2 + p3) / (p4 + p5) + plot_layout(guides = "collect", heights = c(1,1)) + plot_annotation(tag_levels = "A") & theme(legend.position = "bottom") 
ggsave(here::here("outputs", "figs", "fig3.png"), width = 15, height = 15)
require(ggplot2)



p6 <- data.frame(
    `Spike four-fold` = sens_4_s,
    `NCP four-fold` = sens_4_n,
    `Spike seropos` = sens_p_s,
    `NCP seropos` = sens_p_n,
    `serojump` = sens_sj
) %>% pivot_longer(everything()) %>% 
    mutate(name = factor(name, levels = c("Spike.four.fold", "NCP.four.fold", "Spike.seropos", "NCP.seropos", "serojump"))) %>%
    ggplot() + 
        geom_col(aes(x = name, y = value), color = "black") + theme_bw() + theme(text = element_text(size = 20)) +
        guides(fill = "none")

((p1 + p2 + p3) / (p4 + p5) ) / p6 + plot_layout(guides = "collect", heights = c(1,1,2))  + 
    labs(x = "Heuristic used for serologically detecting infection", y = "Sensitivity", title = "Sensitivity of serological detection")
#ggsave(here::here("outputs", "figs", "fig3.png"), width = 20, height = 20)
