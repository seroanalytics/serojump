
generate_convergence_plot <- function(model_summary, title_i) {
    df_smi_df <- calcScaleModelIndicator(model_summary)


    fit <- model_summary$fit
    p1 <- (fit$post$mcmc  %>% bayesplot::mcmc_trace()) + theme_minimal() + theme(legend.position = "top") + ggtitle("A. Trace plots for fitted parameters")
    p2 <- fit$post$lpost %>% ggplot() + geom_line(aes(x = sample_no, y = lpost, color = chain_no))  + theme_minimal() + theme(legend.position = "top")  + ggtitle("B. Trace plots for log posterior") + 
        labs(x = "Sample number", y = "Log-posterior")

    p3 <- df_conver_stat <- summarise_draws(fit$post$mcmc ) %>% select(variable, rhat, ess_bulk, ess_tail) %>% 
        pivot_longer(!variable, names_to = "stat", values_to = "value") %>%
        ggplot() + 
            geom_col(aes(y = variable, x = value)) +
            facet_wrap(~stat, scales = "free") + theme_minimal() + ggtitle("C. Convergence diagnosis for fitted parameters") + 
            labs(x = "Value", y = "Parameter")
    pA <-  p1 / p2 / p3



    pdims_trace <- df_smi_df %>% 
        ggplot() + 
            geom_line(aes(x = .iteration, y = dims, color = as.character(.chain))) + 
            labs(x = "Iteration", y = "Model dimension", color = "Chain")   +
        ggtitle("D. Trace plots for transdimensional convergence: dimensions of model") 


    pdims_hist <- df_smi_df %>% 
        ggplot() + 
            geom_histogram(aes(x = dims, fill = as.character(.chain))) + 
            labs(x = "Model dimension", y = "Count", fill = "Chain") 

    p1 <- pdims_trace + pdims_hist + plot_layout(guides = "collect") & theme_minimal() & theme(legend.position = "top") 



    psmi_trace <- df_smi_df %>% 
        ggplot() + 
            geom_line(aes(x = .iteration, y = sMI, color = as.character(.chain))) + 
            labs(x = "Iteration", y = "Log2 of SMI", color = "Chain")  +
        ggtitle("E. Trace plots for transdimensional convergence: SMI of model") 


    psmi_hist <- df_smi_df %>% 
        ggplot() + 
            geom_histogram(aes(x = sMI, fill = as.character(.chain))) + 
            labs(x = "Log2 of SMI", y = "Count", fill = "Chain") 

    p2 <- psmi_trace + psmi_hist  + plot_layout(guides = "collect") & theme_minimal() & theme(legend.position = "top")

    p3 <- summarise_draws(df_smi_df) %>% select(variable, rhat, ess_bulk, ess_tail) %>% pivot_longer(!variable, names_to = "stat", values_to = "value") %>%
        filter(variable == "sMI") %>% 
        ggplot()+ 
            geom_col(aes(y = "", x = value)) + 
            facet_grid(cols = vars(stat), scales = "free") + theme_minimal() + ggtitle("F. Convergence diagnosis for transdimensional convergence") + 
            labs(x = "Value", y = "Model")

    pB <- p1 / p2 / p3 


    (pA | pB) + plot_annotation(title = paste0("CONVERGENCE DIAGNOSITICS FOR ", title_i)) &
    theme(title = element_text(size = 12))
}


plot_Rhat_time_alt <- function(model_summary, title_i) {

    outputfull <- model_summary$post

    model_outline <- model_summary$fit$model
    bio_all <- model_outline$infoModel$biomarkers

    fit_states_dt <- as.data.table(outputfull$fit_states)
    S <- fit_states_dt %>% filter(id == 1) %>% nrow

    ids <- fit_states_dt %>% group_by(id) %>% summarise(prob = sum(inf_ind) / S) %>% filter(prob > 0.5) %>% pull(id) %>% unique

    if (length(ids) == 0) {
        cat("No individuals have posterior prob of infection > 0.5")
    } else {
        # extract values here
        df_mcmc_time <- fit_states_dt %>% filter(id %in% ids) %>% filter(inf_ind == 1) %>% 
            select(id, chain_no, sample, inf_time, !!bio_all) %>% rename(chain = chain_no) 

        df_mcmc_time_wide <- df_mcmc_time %>% 
            select(id, sample, chain, inf_time) %>% unique %>%
            pivot_wider(!chain, names_from = "id", values_from = "inf_time") 

        cols <- ncol(df_mcmc_time_wide)

        df_summary_disc <- 
                map_df(2:cols,
            ~df_mcmc_time_wide %>% select(sample, .x) %>% drop_na %>% summarise_draws() %>% .[2, ]
        )

        p1 <- df_mcmc_time %>% 
            ggplot() +
                stat_pointinterval(aes(x = inf_time, y = as.character(id), color = as.character(chain)), 
                    position = position_dodge(0.4)) + theme_bw() + 
                    labs(x = "Time in study", y = "ID", color = "Chain number") + 
                    ggtitle("A. Trace plots for timing of infection for individuals \nwith posterior P(Z) > 0.5")

        p2 <- df_summary_disc %>% ggplot() + geom_col(aes(x = rhat, y = as.character(variable))) + theme_bw() + 
            geom_vline(xintercept = 1.1, color = "red", linetype = "dashed") + 
            labs(x = "Rhat", y = "ID") +
            ggtitle("B. Convergence diagnostics for timing of infection \nindividuals with posterior P(Z) > 0.5")

        p1 + p2 + plot_annotation(title = paste0("TIMING CONVERGENCE DIAGNOSITICS FOR ", title_i)) &
            theme(title = element_text(size = 12))

    }

}

# CASE STUDY 1: COP
model_summary <-  readRDS(here::here("outputs", "fits", "simulated_data_hpc", paste0("cop", "_", "0.1"), "model_summary.RDS"))
p1 <- generate_convergence_plot(model_summary, "SIMULATED DATA WITH COP AND 0.1 UNCERTAINTTY IN OBSERVATIONAL ERROR")
ggsave(here::here("outputs", "figs", "supp", "conv", "cop_0.1_full.png"))
p2 <- plot_Rhat_time_alt(model_summary, "SIMULATED DATA WITH COP AND 0.1 UNCERTAINTTY IN OBSERVATIONAL ERROR")
ggsave(here::here("outputs", "figs", "supp", "conv", "cop_0.1_time.png"))

# CASE STUDY 1: NO COP
model_summary <-  readRDS(here::here("outputs", "fits", "simulated_data_hpc", paste0("no_cop", "_", "0.1"), "model_summary.RDS"))
p1 <- generate_convergence_plot(model_summary, "SIMULATED DATA NO COP AND 0.1 UNCERTAINTTY IN OBSERVATIONAL ERROR")
ggsave(here::here("outputs", "figs", "supp", "conv", "no_cop_0.1_full.png"))
p2 <- plot_Rhat_time_alt(model_summary, "SIMULATED DATA NO COP AND 0.1 UNCERTAINTTY IN OBSERVATIONAL ERROR")
ggsave(here::here("outputs", "figs", "supp", "conv", "no_cop_0.1_time.png"))

# CASE STUDY 2: TRANSVIR, NO PCR
model_summary <-  readRDS(here::here("outputs", "fits", "transvir_data", "wave2_no_pcr", "model_summary.RDS"))
p1 <- generate_convergence_plot(model_summary, "EMPIRICAL DATA WITH NO PCR")
ggsave(here::here("outputs", "figs", "supp", "conv", "wave2_no_pcr_full.png"), height = 20)
p2 <- plot_Rhat_time_alt(model_summary,  "EMPIRICAL DATA WITH NO PCR")
ggsave(here::here("outputs", "figs", "supp", "conv", "wave2_no_pcr_time.png"))

# CASE STUDY 2: TRANSVIR, PCR
model_summary <-  readRDS(here::here("outputs", "fits", "transvir_data", "wave2_base", "model_summary.RDS"))
p1 <- generate_convergence_plot(model_summary, "EMPIRICAL DATA WITH PCR")
ggsave(here::here("outputs", "figs", "supp", "conv", "wave2_pcr_full.png"), height = 20)
p2 <- plot_Rhat_time_alt(model_summary,  "EMPIRICAL DATA WITH PCR")
ggsave(here::here("outputs", "figs", "supp", "conv", "wave2_pcr_time.png"))