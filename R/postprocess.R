#' @useDynLib serojump
#' @importFrom Rcpp sourceCpp
#' @import ggplot2
#' @import purrr
#' @import tidyr
#' @import dplyr
#' @import ggdist
#' @import patchwork
#' @import tidybayes
#' @import bench
#' @importFrom magrittr %>% %<>%
NULL

#' @title plotPostFigs
#' @param model_summary Fit from the serojump model
#' @param save_info Information about saving the figures
#' @export
plotPostFigs <- function(model_summary, save_info) {
    check_save_info(save_info)
    
    if (model_summary$fit$data$priorPredFlag) {
        dir.create(here::here("outputs", "fits", save_info$file_name, save_info$model_name,  "figs_pp", "post"), showWarnings = FALSE, recursive = TRUE)
        dir.create(here::here("outputs", "fits", save_info$file_name, save_info$model_name,  "figs_pp", "post", "plt_data"), showWarnings = FALSE, recursive = TRUE)
        file_path <- here::here("outputs", "fits", save_info$file_name, save_info$model_name,  "figs_pp", "post")
    } else{
        dir.create(here::here("outputs", "fits", save_info$file_name, save_info$model_name,  "figs", "post"), showWarnings = FALSE, recursive = TRUE)
        dir.create(here::here("outputs", "fits", save_info$file_name, save_info$model_name,  "figs", "post", "plt_data"), showWarnings = FALSE, recursive = TRUE)
        file_path <- here::here("outputs", "fits", save_info$file_name, save_info$model_name,  "figs", "post")
    }

  #  model_summary <- output_1
    plot_titre_obs(model_summary, file_path)
    plot_titre_exp(model_summary, file_path)
    plot_cop_rec(model_summary, file_path)
    plot_abkinetics_trajectories(model_summary, file_path)
    plot_inf_rec(model_summary, file_path)
    plot_exp_times_rec(model_summary, file_path)
    plot_abkinetics_trajectories_ind(model_summary, file_path)

}

postprocess_fit <- function(model_fit) {

    post <- model_fit$post
    data_t <- model_fit$data_t
    model <- model_fit$model
    settings <- model_fit$settings
    n_chains <- settings$numberChainRuns

    knowndate <- data_t$knownInfsDate
    initialTitreValue <- data_t$initialTitreValue

    knownid <- which(knowndate != -1)
    N <- data_t$N
    N_data <- data_t$N_data

    n_post <- post$mcmc[[1]] %>% nrow
    post_exp_combine <- post$jump 

    post_infexp <- 1:n_chains %>% map_df(
        function(i) {
            1:n_post %>% map_df( 
                function(j) {
                    data.frame(
                        id = 1:N,
                        sample_c = j,
                        chain_no = i,
                        inf_time = post_exp_combine[[i]][j, ],
                        inf_ind = as.numeric(post_exp_combine[[i]][j, ] > 0)
                    )
                }
            )
        }
    )

      # Get titre value at exposure

    post_titreexp <- 1:n_chains %>% map_df(
        function(i) {
            1:n_post %>% map_df( 
                function(j)
                {
                    post_titreexp_j <- post$titreexp[[i]][[j]] %>% as.data.frame %>% 
                        mutate(id = 1:N, sample_c = j, n_chains = i) 
                    post_titreexp_j
                }
            )
        }
    )
    names(post_titreexp) <- c(model$infoModel$biomarkers, "id", "sample_c", "chain_no")

    fit_states <- post_infexp %>% left_join(post_titreexp, by = c("id", "sample_c", "chain_no")) %>% 
        arrange(id, chain_no, sample_c) %>% group_by(id) %>% mutate(sample = row_number()) %>% ungroup

    # Observation model estimates
    post_obstitre <- 1:n_chains %>% map_df(
        function(i) {
            1:n_post %>% map_df( 
                function(j)
                {
                    post_obstitre_j <- post$obstitre[[i]][[j]] %>% as.data.frame %>% 
                        mutate(id = data_t$id_full, row = 1:N_data, sample_c = j, n_chains = i, time = data_t$times_full) 
                    post_obstitre_j <- bind_cols(post_obstitre_j, as.data.frame(data_t$titre_full))
                    post_obstitre_j
                }
            )
        }
    )
    names(post_obstitre) <- c(model$infoModel$biomarkers, "id", "row_id", "sample_c", "chain_no", "time", paste(model$infoModel$biomarkers, "_data", sep = "")) 
 

    postprocess <- list(
        fit_states = fit_states,
        fit_obs = post_obstitre,
        n_chains = n_chains,
        n_post = n_post)
    postprocess
    #if (!priorPred) {
    #    saveRDS(postprocess, file = here::here("outputs", "fits", filename, paste0("pp_", modelname, ".RDS")))
    #} else {
    #    saveRDS(postprocess, file = here::here("outputs", "fits", filename, paste0("pp_prior_", modelname, ".RDS")))
    #}
}


plot_exp_rec <- function(model_summary, file_path) {

    fitfull <- model_summary$fit    
    outputfull <- model_summary$post


    fit_states <- outputfull$fit_states
    filename <- outputfull$filename
    modelname <- outputfull$modelname

    data_t <- fitfull$data_t
    N <- data_t$N

    df_data <- data.frame(
        id = 1:N,
        start_titre = data_t$initialTitreValue
    )

    n_chains <- outputfull$n_chains
    n_post <- outputfull$n_post
    n_length <- n_chains * n_post
    # a) Infection over whole process
    fit_states_exp_prob <- fit_states %>% select(id, exp_ind) %>% summarise(exp_post = mean(exp_ind), .by = "id")

    figA <- fit_states_exp_prob %>% left_join(df_data) %>% 
        ggplot() + 
            geom_point(aes(x = start_titre, y = exp_post), shape = 21, size = 3, alpha = 0.6) + 
            theme_bw() + 
            ylim(0, 1) +
            labs(x = "Titre at start of study", y = "Posterior probability of exposure", fill = "Exposure and infection status") + 
            ggtitle("Posterior plots of probability of exposure compared to simulated data") 

    no_exp_fit_df <- fit_states %>% select(id, sample, exp_ind) %>% summarise(exp_ind_sum = sum(exp_ind), .by = "sample")

    fit_states %>% filter(sample == 2948) %>% as.data.frame

    figB <- no_exp_fit_df %>% summarise(n = n()  / n_length, .by = exp_ind_sum) %>%
        ggplot() + 
            geom_col(aes(x = exp_ind_sum, y = n), alpha = 0.8) + 
            theme_bw()  + 
            labs(y = "Posterior density", x = expression("Estimated number of exposures in epidemic, n"[Z])) + 
            ggtitle("Recovery of population-level exposure burden")

    figB + plot_annotation(tag_levels = "A") 
    ggsave(here::here("outputs", "fits", filename,  "figs", modelname, fig_folder, "exposure_recov.png"), height = 10, width = 10)

}

plot_exp_times_rec <- function(model_summary, file_path) {

    fitfull <- model_summary$fit    
    outputfull <- model_summary$post

    fit_states <- outputfull$fit_states
    filename <- outputfull$filename
    modelname <- outputfull$modelname

    data_t <- fitfull$data_t
    N <- data_t$N

    df_data <- data.frame(
        id = as.character(1:N),
        start_titre = data_t$initialTitreValue,
        known_inf = data_t$knownInfsVec,
        known_date =  data_t$knownInfsTimeVec
    ) %>% mutate(known_inf = ifelse(known_inf == 0, "Not known", "Known"))

    fit_states_exp_prob <- fit_states %>% select(id, inf_ind) %>% summarise(inf_post = mean(inf_ind), .by = "id")

    dataplt <- fit_states %>% filter(inf_ind == 1) %>% select(id, inf_time) %>% group_by(id) %>% mean_qi(inf_time) %>% left_join(fit_states_exp_prob)

    pid_order1 <- dataplt %>% arrange(inf_time) %>% pull(id)

    figB <- dataplt %>% mutate(id = factor(id, levels = pid_order1)) %>%
        left_join(df_data) %>% mutate(id = factor(id, levels = pid_order1)) %>%
        ggplot() + 
            geom_linerange(aes(y = id, xmin = .lower, xmax = .upper, color = inf_post)) + 
            geom_point(aes(y = id, x = inf_time, shape = known_inf, fill = known_inf), color = "black") + 
            scale_fill_manual(values = c("red", "black")) + 
            scale_shape_manual(values = c(21, 22)) + 
            guides(shape = guide_legend(title = "Known infection status"),
                fill = guide_legend(title = "Known infection status")) +
            labs(x = "Time in study", y = "ID", fill = "Posterior probability of infection", color = "Probability of infection") + theme_bw() + 
            theme(axis.text.y = element_blank()) + 
            ggtitle("Recovery of timings of infection over epidemic") 

    
    datasets <- fit_states %>% filter(inf_time > -1, inf_ind == 1)  %>% arrange(sample) %>% 
        split(cumsum(c(TRUE, diff(.$sample) != 0))) %>% map(~pull(.x, inf_time) )
    minmax <- datasets %>% unlist() %>% range
    break_vec <- ceiling(seq(0, minmax[2], length.out = 30))

    hist_data <- datasets %>%
        map(~ hist(., breaks = break_vec, plot = FALSE)) %>%
        map(~ {
            data.frame(
            count = .$counts,
            bin = .$mids
            )
        })
    mean_counts <- hist_data %>%
        map(~ .x %>% pull(count)) %>%  do.call(rbind, .) %>% colMeans
    std_counts <- do.call(rbind, hist_data %>% map(~ .$count)) %>% apply(2, sd)
    
    dataplt_inf <- fit_states %>% filter(inf_time > -1, inf_ind == 1) %>% select(id, inf_time) %>% group_by(id) %>% mean_qi(inf_time)

    df_data_t <- df_data %>% filter(known_inf == "Known") %>% select(id, type = known_inf, time = known_date ) %>% 
        mutate(id = as.numeric(id))
# Plotting
    figC <- df_data_t %>% ggplot() +
      geom_histogram(aes(x = time, y = after_stat(count), fill = "gray40"), size = 2, alpha = 0.5, color = "black") + 
    #geom_line(data = bind_rows(hist_data), aes(x = bin, y = count), alpha = 0.1, color = "blue") +
    geom_line(data = data.frame(bin = hist_data[[1]]$bin, count = mean_counts), aes(x = bin, y = count), color = "red", size = 1.5) +
    geom_ribbon(data = data.frame(bin = hist_data[[1]]$bin, 
                                    lower = mean_counts - std_counts, 
                                    upper = mean_counts + std_counts), 
                aes(x = bin, ymin = lower, ymax = upper, fill = "red"), 
                 alpha = 0.3) +
    scale_fill_manual(
         values =c('gray40'='gray40','red'='red'),
         labels = c('Known infections', 'Inferred total infections')) +
    theme_bw()  + 
            labs(x = "Time in study", y = "Number of people infected", fill = "") + 
            ggtitle("Recovery of infection timings") + 
    theme(legend.position = "bottom")
 
    saveRDS(list(df_data_t, hist_data, mean_counts, std_counts), here::here(file_path, "plt_data", "exposure_time.RDS"))

    figB / figC + plot_annotation(tag_levels = "A") 
    ggsave(here::here(file_path, "exposure_time_recov.png"), height = 10, width = 10)
}

plot_inf_rec <- function(model_summary, file_path) {

    fitfull <- model_summary$fit    
    outputfull <- model_summary$post

    fit_states <- outputfull$fit_states
    filename <- outputfull$filename
    modelname <- outputfull$modelname

    data_t <- fitfull$data_t
    N <- data_t$N

    df_data <- data_t$initialTitreValue %>% as.data.frame %>% mutate( id = 1:N,  known_inf = data_t$knownInfsVec) %>%
        mutate(known_inf = ifelse(known_inf == 0, "Not known", "Known"))

    # b) Infection given exposure 
    n_chains <- outputfull$n_chains
    n_post <- outputfull$n_post
    n_length <- n_chains * n_post

    no_inf_fit_df <- fit_states %>% filter(inf_ind == 1) %>% select(id, sample, inf_ind) %>% summarise(inf_ind_sum = sum(inf_ind), .by = "sample")

    figC <- no_inf_fit_df %>% summarise(n = n()  / n_length, .by = inf_ind_sum) %>%
        ggplot() + 
            geom_col(aes(x = inf_ind_sum, y = n), alpha = 0.8) + 
            theme_bw()  + 
            labs(y = "Posterior density", x = expression("Estimated number of infections in epidemic, n"[Z])) + 
            ggtitle("Recovery of population-level infection burden")

    saveRDS(list(no_inf_fit_df), here::here(file_path, "plt_data", "inf_recov.RDS"))

    figC + plot_annotation(tag_levels = "A")
    ggsave(here::here(file_path, "infection_recov.png"), height = 10, width = 10)
}

calculate_reference_titre_full <- function(model_summary) {
    
    fitfull <- model_summary$fit    
    outputfull <- model_summary$post

    fit_states <- outputfull$fit_states
    filename <- outputfull$filename
    modelname <- outputfull$modelname

    post <- fitfull$post
    data_t <- fitfull$data_t

    n_chains <- outputfull$n_chains
    n_post <- outputfull$n_post
    n_length <- n_chains * n_post
    chain_samples <- 1:n_chains %>% map(~c(rep(.x, n_post))) %>% unlist

    model_outline <- fitfull$model

    post_fit <- post$mcmc %>% lapply(as.data.frame) %>% do.call(rbind, .) %>% as.data.frame %>% mutate(chain = as.character(chain_samples ))

    n_post <- outputfull$n_post
    n_length <- n_chains * n_post
  

   df_full_rf <- fit_states %>% group_by(id) %>% mutate(sample = row_number())

}

#' @export
calculate_reference_titre_expectation <- function(model_summary) {
    
    fitfull <- model_summary$fit    
    outputfull <- model_summary$post

    fit_states <- outputfull$fit_states
    filename <- outputfull$filename
    modelname <- outputfull$modelname

    post <- fitfull$post
    data_t <- fitfull$data_t

    n_chains <- outputfull$n_chains
    n_post <- outputfull$n_post
    n_length <- n_chains * n_post
    chain_samples <- 1:n_chains %>% map(~c(rep(.x, n_post))) %>% unlist

    model_outline <- fitfull$model

    post_fit <- post$mcmc %>% lapply(as.data.frame) %>% do.call(rbind, .) %>% as.data.frame %>% mutate(chain = as.character(chain_samples ))

    n_post <- outputfull$n_post
    n_length <- n_chains * n_post
  

   df_full_info <- fit_states %>% group_by(id) %>% mutate(sample = row_number())


    df_exp_rf <- map_df(1:length(model_outline$observationalModel), 
        function(i) {
        biomarker <- model_outline$observationalModel[[i]]$biomarker

        df_full_info <- fit_states %>% mutate(titre_type = case_when(
            inf_time == -1 ~ "No inf",
            inf_time != -1 ~ "Inf")
        ) %>%  group_by(id, titre_type) %>% add_count(name = "n_rows") %>% 
        group_by(id, titre_type, n_rows) %>% mutate(prop = n_rows / (n_length)) %>%
         group_by(id, titre_type, prop) %>%
        mean_qi(!!sym(biomarker) ) %>% 
        complete(id = 1:data_t$N, titre_type, fill = list(prop = 0)) %>%
        arrange(!!biomarker) %>% as.data.frame 

        # Get response for each individual
        df_full_info_res <- df_full_info %>% filter(titre_type == "Inf") %>% select(id, prop)

        # Get covariate for each individual
        df_full_info_x <- df_full_info %>% group_by(id) %>% mutate(!!sym(biomarker) := weighted.mean(x = !!sym(biomarker), w = prop)) %>% 
         select(id, !!sym(biomarker))

        df_data <- data.frame(
            id = 1:data_t$N,
            start_titre = data_t$initialTitreValue %>% as.data.frame %>% pull(!!biomarker),
            known_inf = data_t$knownInfsVec
        ) %>% mutate(known_inf = ifelse(known_inf == 0, "Not known", "Known"))


        df_cor_info <- df_full_info_res %>% left_join(df_full_info_x) %>% left_join(df_data)  %>%
            rename(titre_val = !!(biomarker)) %>% mutate(biomarker = !!biomarker) 

        })
    df_exp_rf %>% unique
}


plot_cop_rec <- function(model_summary, file_path) {

    cop_exp_sum_plot_all <- calculate_reference_titre_expectation(model_summary)

    saveRDS(list(cop_exp_sum_plot_all), here::here(file_path, "plt_data", "cop_data.RDS"))

   figA <- cop_exp_sum_plot_all %>%
            ggplot() + geom_point(aes(x = titre_val, y = prop)) + theme_bw() + 
            geom_smooth(method = "lm", aes(x = titre_val, y = prop)) +
        ylim(0, 1) +
        facet_wrap(vars(biomarker)) + 
        labs(x = expression("Titre at exposure"), y = expression("Posterior probability of infection"))

    ggsave(here::here(file_path, "cor.png"), height = 10, width = 10)

}


plot_titre_exp <- function(model_summary, file_path) {

    fitfull <- model_summary$fit    
    outputfull <- model_summary$post

    fit_states <- outputfull$fit_states
    filename <- outputfull$filename
    modelname <- outputfull$modelname

    post <- fitfull$post
    data_t <- fitfull$data_t
    n_chains <- outputfull$n_chains
    n_post <- outputfull$n_post
    n_length <- n_chains * n_post
    model_outline <- fitfull$model

    cop_exp_sum_plot_sum <- map_df(1:length(model_outline$observationalModel), 
            function(i) {
            pars_extract <- model_outline$observationalModel[[i]]$pars
            biomarker <- model_outline$observationalModel[[i]]$biomarker

            titre_cop_sum <- fit_states %>% filter(!!sym(biomarker) != -1) %>%
                group_by(id) %>% mean_qi(!!sym(biomarker) ) %>% arrange(!!biomarker) %>%
                rename(titre_val = !!(biomarker))

            u_ids <- titre_cop_sum$id

            df_data_post <- fit_states %>% filter(inf_ind == 1) %>% select(id, inf_ind) %>%
                summarise(inf_post = mean(inf_ind), n = n(), .by = "id") %>%
                mutate(prop = n / n_length)

            cop_exp_sum_plot <- titre_cop_sum %>%
                left_join(df_data_post) %>%
                mutate(id = factor(id, levels = u_ids)) %>% filter(!is.na(n)) %>% 
                mutate(biomarker = biomarker)
            cop_exp_sum_plot
            }
    )

    p1 <- cop_exp_sum_plot_sum %>% ggplot() +
        geom_linerange(aes(y = id, xmin = .lower, xmax = .upper, color = inf_post, alpha = prop), size = 1) + 
        geom_point(aes(y = id, x = titre_val, alpha = prop), size = 1) + 
        labs(y = "ID", x = "Posterior of titre value at exposure", color = "Probability of infection") + 
        theme_bw() + theme(axis.text.y = element_blank()) + facet_wrap(vars(biomarker))
   # if (!is.null(scale_ab)) {
   #     p1 <- p1 + scale_x_continuous(breaks = scale_ab %>% as.numeric, labels = scale_ab %>% names)
   # }
    saveRDS(list(cop_exp_sum_plot_sum), here::here(file_path, "plt_data", "titre_exp_recovery.RDS"))

    p1 
    ggsave(here::here(file_path, "titre_exp_recovery.png"), height = 10, width = 10)

}   

plot_titre_obs <- function(model_summary, file_path) {
    fitfull <- model_summary$fit    
    outputfull <- model_summary$post


    filename <- outputfull$filename
    modelname <- outputfull$modelname

    model_outline <- fitfull$model

    biomarkers <- model_outline$infoModel$biomarkers
    data_plot <- outputfull$fit_obs %>% select(!c(biomarkers, sample_c, chain_no) ) %>% unique %>% 
        pivot_longer(paste0(biomarkers, "_data"), names_to = "biomarker", values_to = "titre") %>%
        mutate(biomarker = stringr::str_remove(biomarker, "_data"))

    
    model_plot <- outputfull$fit_obs %>% select(!paste0(biomarkers, "_data")) %>% 
        pivot_longer(biomarkers, names_to = "biomarker", values_to = "titre") %>%
        group_by(row_id, time, biomarker) %>% mean_qi(titre)

    p1 <- model_plot %>% 
        ggplot() + geom_point(aes(x = row_id, y = titre)) +
            geom_linerange(aes(x = row_id, ymin = .lower, ymax = .upper)) +
            geom_point(data = data_plot, aes(row_id, titre), color = "red") + facet_wrap(vars(biomarker)) 

    saveRDS(list(model_plot), here::here(file_path, "plt_data", "titre_obs.RDS"))

    ggsave(here::here(file_path, "titre_obs.png"), height = 10, width = 10)

}



plot_abkinetics_trajectories <- function(model_summary, file_path) {

    #model_summary <- output_1

    fitfull <- model_summary$fit    
    outputfull <- model_summary$post

    filename <- outputfull$filename
    modelname <- outputfull$modelname

    n_chains <- outputfull$n_chains
    n_post <- outputfull$n_post

    chain_samples <- 1:n_chains %>% map(~c(rep(.x, n_post))) %>% unlist

    data_t <- fitfull$data_t
    N <- data_t$N
    T_max <- (data_t$endTitreTime %>% max) + 10

    post <- fitfull$post
    par_tab <- fitfull$par_tab
    model_outline <- fitfull$model
    post_fit <- post$mcmc %>% lapply(as.data.frame) %>% do.call(rbind, .) %>% as.data.frame %>% mutate(chain = as.character(chain_samples ))


    posteriorsAllExposure <- map_df(1:length(model_outline$abkineticsModel),
        function(name1) {
            pars_extract <- model_outline$abkineticsModel[[name1]]$pars
            functionalForm <- model_outline$abkineticsModel[[name1]]$funcForm
            biomarker <- model_outline$abkineticsModel[[name1]]$biomarker
            exposureType <- model_outline$abkineticsModel[[name1]]$exposureType

           # compare <- bind_rows(
           #     post$mcmc %>% lapply(as.data.frame) %>% do.call(rbind, .) %>% as.data.frame  %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>%
           #         mutate(type = "Posterior distribution") %>% filter(param %in% pars_extract)
#
           # )

            ab_function <- function(T, pars) {
                1:T %>% map( 
                    ~functionalForm(0, .x, pars)
                ) %>% unlist
            }

            ## hierarchical effects
            hierFlag_value <- model_outline$abkineticsModel[[name1]]$hierFlag

            pars_extract_list <- calculate_param_list(pars_extract, hierFlag_value, model_outline, name1)

            df_trajectories <- map_df(1:length(pars_extract_list), 
                function(j) {
                    df_adjusted_par <- extract_parameters(model_outline, name1, post_fit, pars_extract_list, 1:nrow(post_fit), j)# %>% mutate(covar_value = j, exposure_type = exposureType, biomarker = biomarker)
                    T <- T_max
                    traj_post <- 1:(100) %>% purrr::map_df(
                        ~data.frame(
                            time = 1:T,
                            value = ab_function(T, as.numeric(df_adjusted_par[.x, ]))
                        )
                    )  %>% group_by(time) %>% mean_qi() %>% mutate(exposure_type = exposureType) %>% 
                                mutate(biomarker = biomarker, covar = j)
                }
            ) 
        }
    )
    

    saveRDS(list(posteriorsAllExposure), here::here(file_path, "plt_data", "ab_kinetics_recov.RDS"))

    if (max(posteriorsAllExposure$covar) == 1) {
        p1 <- posteriorsAllExposure %>%  
            ggplot() + 
            geom_ribbon(aes(ymin = .lower, ymax = .upper, x = time, fill = exposure_type), size = 3, alpha = 0.3) + 
            geom_line(aes(y = value, x = time, color = exposure_type), size = 2) + 
            theme_bw() + labs(x = "Time post-infection/since bleed, s", y = "Titre change", color = "type", fill = "type") + 
            ggtitle("Antibody kinetic trajectories") +
            facet_wrap(vars(biomarker))
    } else {
        p1 <- posteriorsAllExposure %>%  
            ggplot() + 
            geom_ribbon(aes(ymin = .lower, ymax = .upper, x = time, fill = paste(exposure_type, as.character(covar), sep = "_")), size = 3, alpha = 0.3) + 
            geom_line(aes(y = value, x = time, color =  paste(exposure_type, as.character(covar), sep = "_")), size = 2) + 
            theme_bw() + labs(x = "Time post-infection/since bleed, s", y = "Titre change", color = "type", fill = "type") + 
            ggtitle("Antibody kinetic trajectories") +
            facet_wrap(vars(biomarker))
    }

    p1
    ggsave(here::here(file_path, "ab_kinetics_trajectories.png"), height = 10, width = 10)

}

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## LARGE FUNCTION TO SIMULATE INDIVIDUAL-LEVEL TRAJECTORIES ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# Function to filter and sort data based on ID, sample, and biomarker
filter_sorted_data <- function(data, id_i, sample_i, biomarker_i) {
  data %>%
    filter(id == id_i, sample == sample_i, biomarker == biomarker_i) %>%
    arrange(time)
}

# Extract hierarchical parameter if applicable
extract_hierarchical_data_value <- function(model, id_key, index) {
  hier_flag <- model$abkineticsModel[[id_key]]$hierFlag
  
  if (isTRUE(hier_flag)) {
    return(model$abkineticsModel[[id_key]]$dataHier[index])
  }
  return(1)
}

# Adjust hierarchical parameters
apply_hierarchical_adjustment <- function(post_fit, base_params, hier_params, model_outline, group_idx ) {
    #post_fit = mock_post_fit_1
    #param_list = result_hier[[group_idx]]
    #base_params = mock_model_1$abkineticsModel[[id_key]]$parsBase
    #hier_params =  mock_model_1$abkineticsModel[[id_key]]$parsHier
    #model_summary = mock_model_summary
  
  for (param in hier_params) {
   #param <- hier_params[1]
    boundaries <- model_outline$infoModel$logitBoundaries %>%
      filter(par_name == param)
    upper <- boundaries$ub
    lower <- boundaries$lb

    param_names <- c(param, paste0("z_", param, "_", group_idx), paste0("sigma_", param) )
    
    post_fit <- post_fit %>% mutate(
      !!sym(param_names[1]) := logit_inverse(
        !!sym(param_names[1]) +
        !!sym(param_names[2]) *
        !!sym(param_names[3])
      ) * (upper - lower) + lower
    )
  }
  return(post_fit %>% select(all_of(base_params)))
}

# Calculate extracted parameters list based on hierarchy flag
calculate_param_list <- function(params, hier_flag, model, id_key) {
  if (isTRUE(hier_flag)) {
  #  hier_data <- model$abkineticsModel[[id_key]]$dataHier
    hier_params <- model$abkineticsModel[[id_key]]$parsHier
    base_params <- model$abkineticsModel[[id_key]]$parsBase
    num_groups <- model$abkineticsModel[[id_key]]$dataHierN
    
    param_list <- map(seq_len(num_groups), function(i) {
      map(base_params, ~ if (. %in% hier_params) c(., paste0("z_", ., "_", i), paste0("sigma_", .)) else .) %>% unlist()
    })
  } else {
    param_list <- list(params)
  }
  return(param_list)
}

# Extract parameters based on hierarchy
extract_parameters <- function(model_outline, id_key, samples, param_list, sample_idx, group_idx) {
    hier_flag <- model_outline$abkineticsModel[[id_key]]$hierFlag
    param_list_group <- param_list[[group_idx]]

  if (isTRUE(hier_flag)) {
    hier_params <- model_outline$abkineticsModel[[id_key]]$parsHier
    base_params <- model_outline$abkineticsModel[[id_key]]$parsBase
    return(apply_hierarchical_adjustment(samples, base_params, hier_params, model_outline, group_idx)[sample_idx, ])
  }
  return(samples[sample_idx, param_list[[group_idx]], drop = FALSE])
}

# Compute trajectory for a given time vector
compute_trajectory <- function(anchor_titre, time_vector, model_function, params) {
  map_dbl(seq_len(time_vector), ~ model_function(anchor_titre, .x, params))
}

# Prepare posterior samples from MCMC output
prepare_posterior_samples <- function(post, chain_samples) {
  post$mcmc %>% 
    lapply(as.data.frame) %>% 
    do.call(rbind, .) %>% 
    as.data.frame() %>%
    mutate(chain = as.character(chain_samples))
}

# Compute individual trajectories
compute_individual_trajectories <- function(id, bio_markers, samples, exposure_order, model_outline, param_samples, T_max, data_fit_list) {
    require(data.table)

         #   id = selected_ids
         #   bio_markers = bio_markers
         #   samples = sample_ids
         #   exposure_order = df_inferred_exp
         #   model = params$model_outline
        #    param_samples = param_samples
         #   T_max = params$T_max
        #    data_fit_list = data_fit_list
#filter(id == 162, sample == 29)
  map_df(id, function(i) {
  #  cat("i: ", i, "\n")
   # i <- id[1]
    map_df(bio_markers, function(bio) {
     #   bio <- bio_markers[1]
     map_df(samples, function(s) {
      #  cat("s: ", s, "\n")
     #   s <- samples[1]
        exp_data <- filter_sorted_data(exposure_order, i, s, bio) %>% mutate(row_id = row_number())
             #   cat("exp_data: ", exp_data, "\n")

        times <- c(exp_data$time, T_max)
        time_diff <- diff(times)
        titre_traj <- c()
        titre_traj_type <- c()
        titre_anchor <- NULL
        titre_base <- data_fit_list[[i]] %>% filter(row_num == 1, bio == bio) %>% pull(titre)
        
        for (row_idx in seq_len(nrow(exp_data))) {
            exp_type <- exp_data[row_idx, ]$exp_type
            id_key <- which(model_outline$abkineticsModel %>% map(~.x$exposureType ) %>% unlist() == exp_data[row_idx, ]$exp_type)
            model_func <- model_outline$abkineticsModel[[id_key]]$funcForm
            base_params <- model_outline$abkineticsModel[[id_key]]$pars
            
            group_idx <- extract_hierarchical_data_value(model_outline, id_key, i)

            param_list <- calculate_param_list(base_params, model_outline$abkineticsModel[[id_key]]$hierFlag, model_outline, id_key)
            model_params <- extract_parameters(model_outline, id_key, param_samples, param_list, s, group_idx) %>% as.numeric
            
            if ((exp_data[row_idx, ] %>% pull(row_id)) == 1) {
                titre_traj_type <-  c(rep(NA, max(times[1] - 1, 0)), rep(exp_type, time_diff[row_idx]))
                titre_traj <- rep(NA, times[1])
                new_titre <- compute_trajectory(titre_base, time_diff[row_idx], model_func, model_params)
                titre_traj <- c(titre_traj, new_titre)
                titre_anchor <- ifelse(length(new_titre) == 0, titre_base, new_titre[length(new_titre)])
            } else {
                titre_traj_type <- c(titre_traj_type, rep(exp_type, time_diff[row_idx]))
                new_titre <- compute_trajectory(titre_anchor, time_diff[row_idx], model_func, model_params)
                titre_traj <- c(titre_traj, new_titre)
                titre_anchor <- ifelse(length(new_titre) == 0, titre_anchor, new_titre[length(new_titre)])
            }
        }
        titre_traj_type[length(new_titre)] <- exp_type
    #    cat("titre_traj: ", titre_traj, "\n")
        
        data.table(
          id = i,
          sample = s,
          t = 1:T_max,
          biomarker = bio,
          exp_type = titre_traj_type,
          titre_trajectory = titre_traj
        ) %>% filter(!is.na(titre_trajectory))
      })
    })
  })
}

# model_outline = model_summary$fit$model

# Initialize model parameters
initialize_parameters <- function(model_summary, S) {
  fit <- model_summary$fit
  post <- model_summary$post
  
  list(
    filename = post$filename,
    modelname = post$modelname,
    n_chains = post$n_chains,
    n_post = post$n_post,
    chain_samples = rep(1:post$n_chains, each = post$n_post),
    M = post$n_chains * post$n_post,
    data_t = fit$data_t,
    N = fit$data_t$N,
    T_max = max(fit$data_t$endTitreTime) + 10,
    post = fit$post,
    post_alt = post,
    par_tab = fit$par_tab,
    model_outline = fit$model,
    S = S
  )
}

# Extract exposure order from model outline
get_exposure_order_known <- function(model_outline, N) {
  exposures <- model_outline$infoModel$exposureInfo %>% map(~ .x$exposureType) %>% unlist
  map_df(seq_along(exposures), function(e) {
    data.frame(
      id = 1:N,
      time = model_outline$infoModel$exposureInfo[[e]]$known_inf,
      exp_type = exposures[e]
    )
  }) %>% filter(time > -1)
}

sample_fitted_states <- function(outputfull, sample_s) {
    require(data.table)
    fit_states_dt <- data.table::as.data.table(outputfull$fit_states)
    fit_states_dt %>% filter(sample %in% sample_s)
}


compute_exposure_order <- function(params, bio_all, sample_s) {
   # params
   # bio_all <- bio_markers
    # sample_s <- sample_ids
    N <- params$N
    S <- params$S
    exposures_fit <- params$model_outline$infoModel$exposureFitted
    post <- params$post_alt

    fit_states_dt <- sample_fitted_states(post, sample_s)
    df_know_exp <- get_exposure_order_known(params$model_outline, N)

    future_map(1:N, function(i) {
        df_know_exp_i <- df_know_exp %>% filter(id == i, time > -1, !exp_type %in% exposures_fit)
        fit_i <- fit_states_dt %>% filter(id == i)

        map(1:S, function(s) {
            fit_i_s <- fit_i %>% filter(sample == sample_s[s])

            if (fit_i_s$inf_ind == 1) {
                fit_i_s_long <- fit_i_s %>%
                    pivot_longer(c(!!bio_all), names_to = "biomarker", values_to = "titre") %>%
                    mutate(exp_type = !!exposures_fit) %>%
                    select(id, sample, time = inf_time, exp_type, biomarker, titre)
            } else {
                fit_i_s_long <- data.table()
            }

            if (nrow(df_know_exp_i) > 0) {
                new_rows_list <- lapply(1:nrow(df_know_exp_i), function(j) {
                    data.table(
                        id = i,
                        sample = sample_s[s],
                        time = df_know_exp_i$time[j],
                        exp_type = df_know_exp_i$exp_type[j],
                        biomarker = bio_all,
                        titre = NA
                    )
                })
                fit_i_s_long <- rbindlist(c(list(fit_i_s_long), new_rows_list), fill = TRUE) %>%
                    arrange(time)
            }
            fit_i_s_long
        }) %>% rbindlist
    }) %>% rbindlist
}

extract_serological_data <- function(data_t, bio_all, N) {
    data_fit <- map_df(1:N, function(i) {
        map_df(1:length(bio_all), function(b) {
            titre_b <- c(data_t$titre_list[[i]][[b]])
            times <- c(data_t$times_list[[i]])
            data.frame(
                id = i,
                type = "sero",
                titre = titre_b,
                times = times,
                bio = bio_all[b],
                row_num = 1:length(titre_b)
            )
        })
    })
    data_fit_list <- split(data_fit, data_fit$id)
}

summarize_exposure_intensity <- function(df_exposure_order, params, bio_all) {
    S <- params$S
    exposures_fit <- params$model_outline$infoModel$exposureFitted
    df_exposure_order %>%
        filter(exp_type %in% exposures_fit) %>%
        summarise(prob = n() / (S * length(bio_all)), .by = id) %>%
        mutate(type = cut(prob, c(0, 0.25, 0.75, 1), labels = c("Low", "Medium", "High")))
}

# Main function to plot antibody kinetics trajectories
plot_abkinetics_trajectories_ind <- function(model_summary, file_path, parallel = FALSE) {
    params <- initialize_parameters(model_summary, 100)
    
   # if (parallel) {
   #     plan(multisession, workers = 8)
   # }
    
    # Extract the biomarkers, the sampled numbers and the resulting posterior values 
    bio_markers <- params$model_outline$infoModel$biomarkers
    sample_ids <- sample(1:params$M, params$S)
    param_samples <- prepare_posterior_samples(params$post, params$chain_samples)
    
    # Get the exposure order for each indivdual and each chain in the posterior
    df_inferred_exp <- compute_exposure_order(params, bio_markers, sample_ids)
    df_inferred_exp_sum <- summarize_exposure_intensity(df_inferred_exp, params, bio_markers)

    # Get the original serologicla data
    data_fit_list <- extract_serological_data(params$data_t, bio_markers, params$N)

    # Get the trajectory types
    # need to add low, medium and high into this 
    types_traj <- unique(df_inferred_exp_sum$type)
  
    for (traj_type in types_traj) {
       # traj_type <- types_traj[2]
        cat("Calculating trajectories for", traj_type, "\n")
        
        selected_ids <- df_inferred_exp_sum %>% filter(type == traj_type) %>% pull(id) %>% unique()
        selected_ids <- sample(selected_ids, min(20, length(selected_ids)))

     #   df_inferred_exp %>% filter(id %in% selected_ids) %>% filter(id == 162, sample == 772)
        
        traj_model <- compute_individual_trajectories(
            id = selected_ids,
            bio_markers = bio_markers,
            samples = sample_ids,
            exposure_order = df_inferred_exp,
            model_outline = params$model_outline,
            param_samples = param_samples,
            T_max = params$T_max,
            data_fit_list = data_fit_list
        )

        saveRDS(list(traj_model), here::here(file_path, "plt_data", paste0("ab_kinetics_recov_indiv_ ", traj_type, ".RDS")) )

        data_sero_filter <- data_fit_list %>% bind_rows %>% filter(id %in% selected_ids )
        
        plot <- ggplot(traj_model) +
            geom_point(data = data_sero_filter, aes(x = times, y = titre)) +
            geom_line(  aes(x = t, y = titre_trajectory, color = exp_type, group = sample), alpha = 0.2) +
            facet_wrap(vars(id)) +
            theme_bw() +
            labs(title = paste("Antibody kinetics for", traj_type), x = "Time", y = "Titre value")
            
        ggsave(file.path(file_path, paste0("ab_kinetics_trajectories_", traj_type, ".png")), plot, width = 10, height = 8)
    }
}
