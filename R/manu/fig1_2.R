plotPostFigs <- function(model_summary, save_info) {


    check_save_info(save_info)
    
    if (model_summary$fit$data$priorPredFlag) {
        dir.create(here::here("outputs", "fits", save_info$file_name, save_info$model_name,  "figs_pp", "post"), showWarnings = FALSE, recursive = TRUE)
        file_path <- here::here("outputs", "fits", save_info$file_name, save_info$model_name,  "figs_pp", "post")
    } else{
        dir.create(here::here("outputs", "fits", save_info$file_name, save_info$model_name,  "figs", "post"), showWarnings = FALSE, recursive = TRUE)
        file_path <- here::here("outputs", "fits", save_info$file_name, save_info$model_name,  "figs", "post")
    }

    # Plot CPRS for fitted trajecotries
    # Plot sensitivity and specificity of the recovered posteriors
    # Plot the distribution of timing of cases (K-S test?)
    # Plot the recovery of the COP
    plot_CPRS(model_summary, file_path)

}

devtools::load_all()
library(purrr)

# Plot all summary plots for each model 

values_i <- c(0.01, seq(0.05, 0.5, 0.05))

library(ggridges)
simulated_normal <- map_df(seq_along(values_i),
    ~data.frame(
        sd = values_i[.x],
        values = rnorm(1e4, 0, values_i[.x])
    )
)

require(ggdist)

p0 <- simulated_normal %>% 
   ggplot() + 
        stat_pointinterval(aes(x = values, y = as.character(sd)), point_size = 10, .width = c(0.5, 0.95))  + 
        theme_bw() + 
        theme(text = element_text(size = 15)) + labs(x = "Observational error", y = "Standard deviation of the error") 



vals <- c(values_i, values_i)
types <- c(rep("cop", 11), rep("no_cop", 11))

for (i in 1:22) {
    model_summary_i <-  readRDS(here::here("outputs", "fits", "simulated_data_hpc", paste0(types[i], "_", vals[i]), "model_summary.RDS"))

    ab_kin <-  readRDS(here::here("outputs", "fits", "simulated_data_hpc", paste0(types[i], "_", vals[i]), "figs", "post", "plt_data",  "ab_kinetics_trajectories_High.RDS"))
    exp_times <-  readRDS(here::here("outputs", "fits", "simulated_data_hpc", paste0(types[i], "_", vals[i]), "figs", "post", "plt_data",  "exposure_recov.RDS"))
    exp_times_A <-  readRDS(here::here("outputs", "fits", "simulated_data_hpc", paste0(types[i], "_", vals[i]), "figs", "post", "plt_data",  "exposure_time_recov.RDS"))
    exp_times_B <-  readRDS(here::here("outputs", "fits", "simulated_data_hpc", paste0(types[i], "_", vals[i]), "figs", "post", "plt_data",  "exposure_time_recov_diff.RDS"))

    # trajecotries for subset of individuals 
    p1 <- ab_kin[[1]] %>% 
        ggplot() +
            geom_line(aes(x = t, y = titre_traj, group = sample), color = "#960319", alpha = 0.05) + facet_wrap(vars(id)) + 
            geom_point(aes(x = times, y = titre), data = ab_kin[[2]] %>% filter(id %in% unique(ab_kin[[1]]$id)), shape = 21, size = 2, fill = "gray") + theme_bw() + 
            labs(x = "Time in study", y = "Titre value", color = "Exposure type") + ggtitle(paste0("Individual-level antibody kinetics for simulated titre (subset)"))  + 
            theme(text = element_text(size = 12))


    id_figD <- exp_times_B[[1]] %>% pull(id)
    p2 <- exp_times_B[[1]] %>% mutate(id = factor(id, levels = id_figD)) %>% 
        ggplot() + 
            geom_hline(yintercept = c(14, -14), size = 0.5, color = "gray") +
            geom_linerange(aes(x = id, ymin = .lower, ymax = .upper, color = inf_post), 
                    size = 3) + 
            geom_point(aes(x = id, y = inf_time), size = 3) + 
            labs(y = "Difference between model-predicted infection day \n time and simulated infection day", x = "ID", color = "Posterior probability of exposure") + theme_bw() + 
            theme(axis.text.x = element_blank(), legend.position = "bottom") + 
            ggtitle("Model error in recovering exposure times of infected individuals")  + 
            theme(text = element_text(size = 12)) 
        

    p3 <- exp_times_A[[2]] %>%
        ggplot() +
            geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") + 
            geom_linerange(aes(x = inf_time_sim, y = inf_time, ymin =  .lower, ymax = .upper), color = "black", alpha = 0.7 ) + 
            geom_point(aes(x = inf_time_sim, y = inf_time), color = "black", alpha = 0.7 ) + 
            labs(y = "Density", x = "Infection time") +
            ggtitle("Distribution of infection times between data and model") + 
            theme(text = element_text(size = 12))


    p0 <- p1 + (p2 / p3) + plot_annotation(tag_level = "A",
        title = paste0("Simulation recovery for ", types[i], " with observational error ", vals[i]), 
        theme = theme(plot.title = element_text(size = 20)))
    ggsave(here::here("outputs", "fits", "simulated_data_hpc", paste0(types[i], "_", vals[i]), "figs", "post", "manu_figure1.pdf"), height = 10, width = 15)
}


## Figure 1 

i <- 3
model_summary_i <-  readRDS(here::here("outputs", "fits", "simulated_data_hpc", paste0(types[i], "_", vals[i]), "model_summary.RDS"))

ab_kin <-  readRDS(here::here("outputs", "fits", "simulated_data_hpc", paste0(types[i], "_", vals[i]), "figs", "post", "plt_data",  "ab_kinetics_trajectories_High.RDS"))
exp_times <-  readRDS(here::here("outputs", "fits", "simulated_data_hpc", paste0(types[i], "_", vals[i]), "figs", "post", "plt_data",  "exposure_recov.RDS"))
exp_times_A <-  readRDS(here::here("outputs", "fits", "simulated_data_hpc", paste0(types[i], "_", vals[i]), "figs", "post", "plt_data",  "exposure_time_recov.RDS"))
exp_times_B <-  readRDS(here::here("outputs", "fits", "simulated_data_hpc", paste0(types[i], "_", vals[i]), "figs", "post", "plt_data",  "exposure_time_recov_diff.RDS"))

write.csv(ab_kin[[2]], here::here("outputs", "fits", "simulated_data_hpc", "export", "sim_sero.csv") )


# trajecotries for subset of individuals 
p1A <- ab_kin[[1]] %>% 
    ggplot() +
        geom_line(aes(x = t, y = titre_traj, group = sample), color = "#960319", alpha = 0.05) + facet_wrap(vars(id)) + 
        geom_point(aes(x = times, y = titre), data = ab_kin[[2]] %>% filter(id %in% unique(ab_kin[[1]]$id)), shape = 21, size = 2, fill = "gray") + theme_bw() + 
        labs(x = "Time in study", y = "Titre value (log)", color = "Exposure type") + ggtitle(paste0("Individual-level antibody kinetics for simulated titre\n (with COP)"))  + 
        theme(text = element_text(size = 12))


id_figD <- exp_times_B[[1]] %>% pull(id)
p2A <- exp_times_B[[1]] %>% mutate(id = factor(id, levels = id_figD)) %>% 
    ggplot() + 
        geom_hline(yintercept = c(14, -14), size = 0.5, color = "gray") +
        geom_linerange(aes(x = id, ymin = .lower, ymax = .upper, color = inf_post), 
                size = 3) + 
        geom_point(aes(x = id, y = inf_time), size = 3) + 
        labs(y = "Difference between model-predicted \ninfection day time and \nsimulated infection day", x = "ID", color = "Posterior probability of exposure") + theme_bw() + 
        theme(axis.text.x = element_blank(), legend.position = "bottom") + 
        ggtitle("Model error in recovering infection times (with COP)")  + 
        theme(text = element_text(size = 12)) + guides(color = "none")
    

p3A <- exp_times_A[[2]] %>% 
    ggplot() +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") + 
        geom_linerange(aes(x = inf_time_sim, y = inf_time, ymin =  .lower, ymax = .upper), color = "black", alpha = 0.7 ) + 
        geom_point(aes(x = inf_time_sim, y = inf_time), color = "black", alpha = 0.7 ) + 
        labs(y = "Recovered infection time", x = "Simulated infection time", fill = "") +
        theme_bw() + 
        scale_fill_manual(values = c("Data (with COP)" = "gray", "Model" = "#960319")) +
        ggtitle("Distribution of infection times between \nsimulated data (with COP) and model") + 
        theme(text = element_text(size = 12))


i <- 14

model_summary_i <-  readRDS(here::here("outputs", "fits", "simulated_data_hpc", paste0(types[i], "_", vals[i]), "model_summary.RDS"))

ab_kin <-  readRDS(here::here("outputs", "fits", "simulated_data_hpc", paste0(types[i], "_", vals[i]), "figs", "post", "plt_data",  "ab_kinetics_trajectories_High.RDS"))
exp_times <-  readRDS(here::here("outputs", "fits", "simulated_data_hpc", paste0(types[i], "_", vals[i]), "figs", "post", "plt_data",  "exposure_recov.RDS"))
exp_times_A <-  readRDS(here::here("outputs", "fits", "simulated_data_hpc", paste0(types[i], "_", vals[i]), "figs", "post", "plt_data",  "exposure_time_recov.RDS"))
exp_times_B <-  readRDS(here::here("outputs", "fits", "simulated_data_hpc", paste0(types[i], "_", vals[i]), "figs", "post", "plt_data",  "exposure_time_recov_diff.RDS"))

# trajecotries for subset of individuals 
p1B <- ab_kin[[1]] %>% 
    ggplot() +
        geom_line(aes(x = t, y = titre_traj, group = sample), color = "#960319", alpha = 0.05) + facet_wrap(vars(id)) + 
        geom_point(aes(x = times, y = titre), data = ab_kin[[2]] %>% filter(id %in% unique(ab_kin[[1]]$id)), shape = 21, size = 2, fill = "gray") + theme_bw() + 
        labs(x = "Time in study", y = "Titre value (log)", color = "Exposure type") + ggtitle(paste0("Individual-level antibody kinetics for simulated titre\n (no COP)"))  + 
        theme(text = element_text(size = 12))


id_figD <- exp_times_B[[1]] %>% pull(id)
p2B <- exp_times_B[[1]] %>% mutate(id = factor(id, levels = id_figD)) %>% 
    ggplot() + 
        geom_hline(yintercept = c(14, -14), size = 0.5, color = "gray") +
        geom_linerange(aes(x = id, ymin = .lower, ymax = .upper, color = inf_post), 
                size = 3) + 
        geom_point(aes(x = id, y = inf_time), size = 3) + 
        labs(y = "Difference between model-predicted \ninfection day time and \nsimulated infection day", x = "ID", color = "Posterior probability of exposure") + theme_bw() + 
        theme(axis.text.x = element_blank(), legend.position = "bottom") + 
        ggtitle("Model error in recovering infection times (no COP)")  + 
        theme(text = element_text(size = 12)) + guides(color = "none")
    

p3B <- exp_times_A[[2]] %>%
    ggplot() +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") + 
        geom_linerange(aes(x = inf_time_sim, y = inf_time, ymin =  .lower, ymax = .upper), color = "black", alpha = 0.7 ) + 
        geom_point(aes(x = inf_time_sim, y = inf_time), color = "black", alpha = 0.7 ) + 
        labs(y = "Recovered infection time", x = "Simulated infection time", fill = "") +
        theme_bw() + 
        scale_fill_manual(values = c("Data (with COP)" = "gray", "Model" = "#960319")) +
        ggtitle("Distribution of infection times between \nsimulated data (with COP) and model") + 
        theme(text = element_text(size = 12))


(p1A | p1B) / ( (p2A / p2B) | (p3A / p3B)) + plot_annotation(tag_level = "A") & theme(text = element_text(size = 15))
ggsave(here::here("outputs", "figs", "fig1.png"), width = 15, height = 17)
ggsave(here::here("outputs", "figs", "fig1.tiff"), width = 15, height = 17)




####################################################################################
############################# Obsesrvation model ############################
####################################################################################

library(loo)

crps_res <- list()


for (i in 1:22) {
    model_summary_i <-  readRDS(here::here("outputs", "fits", "simulated_data_hpc", paste0(types[i], "_", vals[i]), "model_summary.RDS"))

    pred <- model_summary_i$post$fit_obs %>% select(sVNT, row_id, sample_c, chain_no) %>% 
        unite("sample", c(sample_c, chain_no) ) %>%
        pivot_wider(names_from = row_id, values_from = sVNT) %>% select(!sample) %>% as.matrix
    pred_trim <- pred[, -seq(1, 1000, 5) ] # remove starting value

    obs <- model_summary_i$post$fit_obs %>% select(row_id, sVNT_data) %>% unique %>% pull(sVNT_data)
    obs_trim <- obs[-seq(1, 1000, 5)] # remove starting value
    crps_res[[i]] <- loo::crps(pred_trim, pred_trim, obs_trim)

}


df_full_ab <- map_df(1:22, 
    ~data.frame(
        uncert = vals[.x],
        type = types[.x],
        crps = 
            (1 - exp(crps_res[[.x]]$pointwise %>% mean))))


df_full_ab_mean <- data.frame(
    uncert = vals,
    type = types,
    mean_crps = sapply(crps_res, function(x)
        (1 - exp(x$estimates[1]))),
    sd_crps = sapply(crps_res, function(x)
       exp(x$estimates[2]) - 1))
    


p1 <- df_full_ab_mean %>% 
    ggplot() + 
        geom_line(aes(x = uncert, y = 1 - mean_crps, color = types), size = 4, alpha = 0.8, position = position_dodge(0.05)) + theme_bw() +
        geom_ribbon( aes(x = uncert, ymin = 1 - mean_crps - 2*sd_crps,  ymax = 1 - mean_crps + 2*sd_crps, fill = types), size = 2, alpha = 0.2, position = position_dodge(0.05)) + theme_bw() +
      #  geom_linerange( aes(x = uncert, ymin = mean_crps - 2*sd_crps,  ymax = mean_crps + 2*sd_crps, color = types), size = 2, alpha = 0.8, position = position_dodge(0.05)) + theme_bw() +
        labs(y = "Mean 1 - CRPS of observational across all bleeds", x = "Simulated uncertainty in the observational model", 
        color = "Model type") + 
        guides(fill = "none") + 
        scale_color_manual(values = c("#434c69", "#c09741")) +
        scale_fill_manual(values = c("#434c69", "#c09741")) +
        ggtitle("Accuracy in recovering antibody kinetics") + 
        scale_x_continuous(breaks = seq(0, 1, 0.1)) 



####################################################################################
############################# Sensitiviy and spec of infection status ############################
####################################################################################

sim_model_cop <- readRDS(file = here::here("outputs", "sim_data", "cesCOP_notd", "inputs.RDS"))
sim_model_no_cop <- readRDS(file = here::here("outputs", "sim_data", "cesNoCOP_notd", "inputs.RDS"))

metric_info <- list()
roc_data <- list()

vals <- c(values_i, values_i)
types <- c(rep("cop", 11), rep("no_cop", 11))


for(i in 1:22)  {
        model_summary_i <-  readRDS(here::here("outputs", "fits", "simulated_data_hpc", paste0(types[i], "_", vals[i]), "model_summary.RDS"))

        post_status <- model_summary_i$post$fit_states %>% 
            group_by(id) %>% summarise(post_ind = sum(inf_time > -1)) %>% as.data.frame %>%
            mutate(post_ind = round(post_ind / 1000))

        if (i <= 11) {
            true_linker <- data.frame(
                id = 1:200,
                true_inf = sim_model_cop$simpar$foe_pars %>% apply(1, sum) 
            )
        } else {
            true_linker <- data.frame(
                id = 1:200,
                true_inf = sim_model_no_cop$simpar$foe_pars %>% apply(1, sum) 
            )
        }

        post_status_linked <- post_status %>% left_join(true_linker)

        library("ROCit")

        measure <- measureit(score = post_status_linked$post_ind, class = post_status_linked$true_inf,
                            measure = c("ACC", "SENS", "SPEC", "FSCR"))

        roc_empirical <- rocit(score = post_status_linked$post_ind, class = post_status_linked$true_inf, negref = "0") 

        roc_data_tmp <- data.frame(
            TPR = roc_empirical$TPR,
            FPR = roc_empirical$FPR,
            uncert = vals[i],
            type = types[i]
        )
        roc_data <- bind_rows(roc_data, roc_data_tmp)

        vector_vals <- map(1:10, ~measure[[.x]][2])
        names(vector_vals) <- names(measure)
        metric_info_tmp <- vector_vals %>% as.data.frame %>% mutate(uncert = vals[i], type = types[i])
        metric_info <- bind_rows(metric_info, metric_info_tmp)

    }

p2A <- roc_data %>% filter(type == "cop") %>%
    ggplot() + 
        geom_line(aes(x = FPR, TPR, color = uncert, group = uncert), size = 1.8, alpha = 0.6) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray30") +
        theme_bw() + ggtitle("Simulated data with COP") + 
        labs(x = "False Positive Rate", y = "True Positive Rate", color = "Uncertainty") 


p2B <- roc_data %>% filter(type == "no_cop") %>%
    ggplot() + 
        geom_line(aes(x = FPR, TPR, color = uncert, group = uncert), size = 1.8, alpha = 0.6) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray30") +
        theme_bw() + ggtitle("Simulated data without COP") + 
        labs(x = "False Positive Rate", y = "True Positive Rate", color = "Uncertainty") 




p2C <- metric_info %>% 
    ggplot() + 
        geom_line(aes(x = uncert, y = FSCR, color = type), size = 4, alpha = 0.8) + theme_bw()  +
        labs(y = "F1-Score", x = "Simulated uncertainty in the observational model", 
        color = "Model type") + 
        scale_x_continuous(breaks = seq(0, 1, 0.1)) + ggtitle("Accuracy in recovering infection status") + 
            scale_color_manual(values = c("#434c69", "#c09741")) 



p2 <- (p2A / p2B +
    plot_layout(guides = "collect")) | p2C +
    plot_layout(guides = "collect") + plot_annotation(title = "Recovery of the correlate of risk with different levels of uncertainty") 

ggsave(here::here("outputs", "figs", "inf_status.png"), width = 20, height = 10)




####################################################################################
############################# Timing of infection recovery ############################
####################################################################################

sim_model_cop <- readRDS(file = here::here("outputs", "sim_data", "cesCOP_notd", "inputs.RDS"))
sim_model_no_cop <- readRDS(file = here::here("outputs", "sim_data", "cesNoCOP_notd", "inputs.RDS"))

crps_time <- list()

vals <- c(values_i, values_i)
types <- c(rep("cop", 11), rep("no_cop", 11))


for(i in 1:22) {
        cat(i, "\n")
        model_summary_i <-  readRDS(here::here("outputs", "fits", "simulated_data", paste0(types[i], "_", vals[i]), "model_summary.RDS"))

        post_status <- model_summary_i$post$fit_states %>% 
            filter( inf_time > -1) %>%
            select(id, inf_time, sample)
            
    
        if (i <= 11) {

            result <- apply( sim_model_cop$simpar$foe_pars[,,1], 1, 
            function(x) {
                y = which(x == 1)
                if(length(y) == 0) {
                    -1
                } else {
                    y
                }
            }) %>% unlist
            true_linker <- data.frame(
                id = 1:200,
                true_time = result
            )
        } else {    
            result <- apply( sim_model_no_cop$simpar$foe_pars[,,1], 1, 
            function(x) {
                y = which(x == 1)
                if(length(y) == 0) {
                    -1
                } else {
                    y
                }
            }) %>% unlist
            true_linker <- data.frame(
                id = 1:200,
                true_time = result
            )
        }

        post_status_linked <- post_status %>% left_join(true_linker)

        pred <- post_status_linked %>% select(!true_time) %>%
            pivot_wider(names_from = id, values_from = inf_time) %>% as.data.frame %>% select(!c(sample)) %>% as.matrix


      #  post_status_linked <- post_status %>% left_join(true_linker) 
        obs <- post_status_linked %>% select(id, true_time) %>% unique %>% pull(true_time)
        crps_time[[i]] <- loo::crps(pred, pred, obs)
    }

df_full_time <- map_df(1:22, 
    ~data.frame(
        uncert = vals[.x],
        type = types[.x],
        mean_crps = 
            (1 - exp(crps_time[[.x]]$pointwise)) %>% if_else(is.na(.), 1, .))
)

df_full_time_mean <- data.frame(
    uncert = vals,
    type = types,
    mean_crps = sapply(crps_time, function(x)
        (1 - exp(x$pointwise)) %>% if_else(is.na(.), 1, .) %>% mean ),
    sd_crps = sapply(crps_res, function(x)
       x$estimates[2])
)

p3 <- df_full_time_mean %>% 
    ggplot() + 
        geom_line(aes(x = uncert, y = 1 - mean_crps, color = types), size = 4, alpha = 0.8, position = position_dodge(0.05)) + theme_bw() +
        geom_ribbon( aes(x = uncert, ymin = 1 - mean_crps - 2*sd_crps,  ymax = 1 - mean_crps + 2*sd_crps, fill = types), size = 2, alpha = 0.2, position = position_dodge(0.05)) + theme_bw() +
        labs(y = "Mean 1 - CRPS of infection times", x = "Simulated uncertainty in the observational model", 
        color = "Model type") + 
         guides(fill = "none") + 
        ggtitle("Accuracy in recovering epidemic curve") + 
        scale_color_manual(values = c("#434c69", "#c09741")) +
        scale_fill_manual(values = c("#434c69", "#c09741")) +
        scale_x_continuous(breaks = seq(0, 1, 0.1)) 



p1 <- p1 + p3
ggsave(here::here("outputs", "figs", "crps_obs_time.png"), width = 20, height = 10)



####################################################################################
############################# Plot the recovery of the COP ############################
####################################################################################



library(rstan)


biomarker_protection <- function(biomarker_quantity, biomarker_prot_midpoint, biomarker_prot_width) {
    risk <- 1 - 1/(1 + exp(biomarker_prot_width * (biomarker_quantity - biomarker_prot_midpoint)))
    return(risk)
}

# with cop
scalar_cop <-  sim_model_cop$simpar$N_exp / 200
obs_cop <-  (biomarker_protection(seq(0, 5, length.out = 100), 2, 2))
result_inf_cop <- apply( sim_model_cop$simpar$foe_pars[,,1], 1, 
function(x) {
    y = which(x == 1)
    if(length(y) == 0) {
        0
    } else {
        1   
    }
}) %>% unlist


# withput cop
scalar_no_cop <-  sim_model_no_cop$simpar$N_exp / 200
obs_no_cop <- ( biomarker_protection(seq(0, 5, length.out = 100), 2, 0))
result_inf_no_cop <- apply( sim_model_no_cop$simpar$foe_pars[,,1], 1, 
function(x) {
    y = which(x == 1)
    if(length(y) == 0) {
        0
    } else {
        1   
    }
}) %>% unlist


crps_cop <- list()
df_post <- list()
for (i in 1:22) {
    model_summary_i <-  readRDS(here::here("outputs", "fits", "simulated_data_hpc", paste0(types[i], "_", vals[i]), "model_summary.RDS"))
    extract_info <- calculate_reference_titre_expectation(model_summary_i) 

    if (i <= 11) {
        scalar_i <- scalar_cop
    } else {
        scalar_i <- scalar_no_cop
    }

    data_list <- list(
        N = 200,
        L = scalar_i,
        x = extract_info$titre_val,
        y = extract_info$prop
    )

    stan_code <- "data {
            int<lower=1> N;          // Number of observations
            vector[N] x;             // Predictor variable (e.g., titre)
            vector[N] y;             // Continuous response variable
            real L;
        }
        transformed data {
            
            real x_min = min(x);
            real x_max = max(x);
            real step = (x_max - x_min) / 99;
            real midpoint = (x_min + x_max) / 2;
        }
        parameters {
            real<lower = 0> k;                  // Steepness of the curve
            real<lower = x_min, upper = x_max> x0;                 // Midpoint (inflection point)
            real<lower=0> sigma;     // Standard deviation of errors
        }
        model {
            vector[N] y_hat;
            
            // Logistic function for curve fitting
            for (n in 1:N) {
                y_hat[n] = L * (1 - 1 / (1 + exp(-k * (x[n] - x0))));
            }

            // Likelihood: Assume normal residuals
            y ~ normal(y_hat, sigma);

            // Priors
            L ~ uniform(0, 1);      // Prior for upper asymptote
            k ~ normal(0, 5);        // Prior for steepness
            x0 ~ normal(midpoint, midpoint / 2);       // Prior for midpoint
            sigma ~ normal(0, 1);    // Prior for standard deviation
        }
        generated quantities {
            vector[100] y_hat_rel;
            vector[100] y_hat_new;
            vector[100] y_prot;
            vector[100] x_new;
            
            // Posterior prediction
            for (i in 1:100) {
            x_new[i] = x_min + (i - 1) * step;
            y_prot[i] =  1 / (1 + exp(-k * (x_new[i] - x0)));
            y_hat_new[i] = L * (1 - 1 / (1 + exp(-k * (x_new[i] - x0))));
            y_hat_rel[i] = y_hat_new[i] / y_hat_new[1];
            }
        }"

   

    fit <- stan(model_code = stan_code, data = data_list, iter = 2000, chains = 4, cores = 4)

    mat_extract <- rstan::extract(fit, pars = "y_prot")$y_prot %>% as.matrix

    colnames(mat_extract) <- seq(min(data_list$x), max(data_list$x), length.out = 100)
    df_post_tmp <- mat_extract %>% as.data.frame %>% pivot_longer(everything()) %>% mutate(type = types[i], uncert = vals[i])
    df_post <- bind_rows(df_post, df_post_tmp)
    if (i <= 11) {
        obs_i <- obs_cop
    } else {
        obs_i <- obs_no_cop
    }

    crps_cop[[i]] <- crps(mat_extract, mat_extract, obs_i)

}

data_plot <- data.frame(
    name = seq(0, 5, length = 100),
    obs_cop = obs_cop,
    obs_no_cop = obs_no_cop
)

p4A <- df_post %>% group_by(name, type, uncert) %>% summarise(mean = mean(value)) %>% mutate(name = as.numeric(name)) %>%
    filter(type == "cop") %>%
    ggplot() + 
        geom_line(aes(x = name, y = mean, color = uncert, group = uncert), size = 1.8, alpha = 0.6) + theme_bw() +
        geom_line(data = data_plot, aes(x = name, y = obs_cop), color = "red", linetype = "dashed", size = 2, alpha = 0.7) +
        labs(y = "Protection probability", x = "Titre value at exposure (log)", 
        color = "Uncertainty") + 
        scale_x_continuous(breaks = seq(0, 5, 0.5)) + ylim(0, 1) + ggtitle("Simulated data with COP")


p4B <- df_post %>% group_by(name, type, uncert) %>% summarise(mean = mean(value)) %>% mutate(name = as.numeric(name)) %>%
    filter(type == "no_cop") %>%
    ggplot() + 
        geom_line(aes(x = name, y = mean, color = uncert, group = uncert), size = 1.8, alpha = 0.6) + theme_bw() +
        geom_line(data = data_plot, aes(x = name, y = obs_no_cop), color = "red", linetype = "dashed", size = 2, alpha = 0.7) +
        labs(y = "Protection probability", x = "Titre value at exposure (log)", 
        color = "Uncertainty") + 
        scale_x_continuous(breaks = seq(0, 5, 0.5)) + ylim(0, 1) + ggtitle("Simulated data without COP")

p4C <- data.frame(
    uncert = vals,
    type = types,
    cprs = crps_cop %>% map_dbl(~(1 - exp(.x$estimates[1]))) 
) %>% 
    ggplot() + 
        geom_line(aes(x = uncert, y = 1 - cprs, color = type), size = 4, alpha = 0.8) + theme_bw() + 
        labs(y = "Mean 1 - CRPS of correlate of risk", x = "Simulated uncertainty in the observational model", 
        color = "Model type")  + ggtitle("Accuracy in recovering correlate of protection") +
        scale_color_manual(values = c("#434c69", "#c09741")) 


p4 <- (p4A / p4B +
    plot_layout(guides = "collect")) | p4C +
    plot_layout(guides = "collect") + plot_annotation(title = "Recovery of the correlate of protection with different levels of uncertainty") 

p4
ggsave(here::here("outputs", "figs", "cor_rec.png"), width = 20, height = 10)


p1 / p2 / p4 + plot_annotation(tag_level = "A") & theme(text = element_text(size = 15))
ggsave(here::here("outputs", "figs", "crps_all.png"), width = 15, height = 18)














 stan_code <- "
        data {
            int<lower=0> N;
            real scal;
            vector[N] x;
            array[N] int y;
            int N_pred;
            vector[N_pred] x_pred; 
        }
        parameters {
            real beta_0;
            real beta_1;
        }
        model {
            y ~ bernoulli(scal * inv_logit(beta_0 + beta_1 * x));
            beta_0 ~ normal(0, 1.8);
            beta_1 ~ normal(0, 2);
        } 
        generated quantities {
            vector[N_pred] y_pred;
            for (n in 1:N_pred)
                y_pred[n] = scal * inv_logit(beta_0 + beta_1 * x_pred[n]);
        }
    "