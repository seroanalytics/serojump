#' @useDynLib serojump
#' @importFrom Rcpp sourceCpp
#' @import ggplot2
#' @import purrr
#' @import tidyr
#' @import dplyr
#' @import ggdist
#' @import patchwork
#' @import tidybayes
#' @import furrr
#' @importFrom magrittr %>% %<>%
NULL




#' Postprocess the posteriors of the simulation study
#' @param modelname The name of the model
#' @param modelname_sim The name of the simulation model
#' @param obs_er The observation error
#' @param n_chains The number of chains
#' @export 
postprocessFigs_simInf <- function(modelname, modelname_sim, obs_er, n_chains) {

    postprocess_run_simInf( modelname, modelname_sim, obs_er, n_chains, TRUE)
    postprocess_run_simInf( modelname, modelname_sim, obs_er, n_chains)


    outputfull <- readRDS(file = here::here("outputs", "fits", "local", "sim_study_inf", modelname_sim, modelname, paste0("pp_obs_", obs_er, ".RDS")))
    modelname_sim_local <- outputfull$modelname_sim_local
    fitfull <- readRDS(here::here("outputs", "fits", modelname_sim_local, modelname, paste0("fit_", "obs_", obs_er, ".RDS")))
    fig_folder <- "figs"

    plot_titre_obs_simInf(outputfull, fitfull, fig_folder)
    plot_abkinetic_simInf(outputfull, fitfull, fig_folder)
   # plot_inf_rec_sim(outputfull)
    plot_inf_rec_simInf(outputfull, fitfull, fig_folder)
    plot_exp_times_rec_simInf(outputfull, fitfull, fig_folder)
    plot_cop_rec_simInf(outputfull, fitfull, fig_folder)

    outputfull <- readRDS(file = here::here("outputs", "fits", "local", "sim_study_inf", modelname_sim, modelname, paste0("pp_prior_obs_", obs_er, ".RDS")))
    modelname_sim_local <- outputfull$modelname_sim_local
    fitfull <- readRDS(here::here("outputs", "fits", modelname_sim_local, modelname, paste0("fit_prior_", "obs_", obs_er, ".RDS")))
    fig_folder <- "figs_pp"

    plot_titre_obs_simInf(outputfull, fitfull, fig_folder)
    plot_abkinetic_simInf(outputfull, fitfull, fig_folder)
   # plot_inf_rec_sim(outputfull)
    plot_inf_rec_simInf(outputfull, fitfull, fig_folder)
    plot_exp_times_rec_simInf(outputfull, fitfull, fig_folder)
    plot_cop_rec_simInf(outputfull, fitfull, fig_folder)
}

postconvergeFigs_simInf <- function(modelname, modelname_sim, obs_er) {

    outputfull <- readRDS(file = here::here("outputs", "fits", "local", "sim_study_inf", modelname_sim, modelname, paste0("pp_obs_", obs_er, ".RDS")))
    modelname_sim <- outputfull$modelname_sim
    modelname_sim_local <- outputfull$modelname_sim_local
    fitfull <- readRDS(here::here("outputs", "fits", modelname_sim_local, modelname, paste0("fit_", "obs_", obs_er, ".RDS")))
    filepath <- here::here("outputs", "fits", modelname_sim_local, modelname, "figs", paste0("obs_", obs_er) )

    seroJumpPostDaig(fitfull, filepath)

    outputfull <- readRDS(file = here::here("outputs", "fits", "local", "sim_study_inf", modelname_sim, modelname, paste0("pp_prior_obs_", obs_er, ".RDS")))
    modelname_sim <- outputfull$modelname_sim
    modelname_sim_local <- outputfull$modelname_sim_local
    fitfull <- readRDS(here::here("outputs", "fits", modelname_sim_local, modelname, paste0("fit_prior_", "obs_", obs_er, ".RDS")))
    filepath <- here::here("outputs", "fits", modelname_sim_local, modelname, "figs_pp", paste0("obs_", obs_er) )

    seroJumpPostDaig(fitfull, filepath)

}



postprocess_run_simInf <- function(modelname,  modelname_sim, obs_er, n_chains, priorPred = FALSE) {

    dir.create(here::here("outputs", "sim_data", modelname), showWarnings = FALSE)
    dir.create(here::here("outputs", "sim_data", modelname, "figs"), showWarnings = FALSE)
    dir.create(here::here("outputs", "sim_data", modelname, "figs", paste0("obs_", obs_er)), showWarnings = FALSE)
    modelA <- readRDS(here::here("outputs", "sim_data", modelname_sim, "inputs.RDS"))
    res <- readRDS(file = here::here("outputs", "sim_data", modelname_sim, paste0("sim_data_", obs_er, ".rds")))

    if (!priorPred) {
        fitfull <- readRDS(here::here("outputs", "fits", "local", "sim_study_inf", modelname_sim, modelname, paste0("fit_", "obs_", obs_er, ".RDS")))
    } else {
        fitfull <- readRDS(here::here("outputs", "fits", "local", "sim_study_inf", modelname_sim, modelname, paste0("fit_prior_", "obs_", obs_er, ".RDS")))
    }

    post <- fitfull$post
    data_t <- fitfull$data_t
    knowndate <- data_t$knownInfsDate
    model <- fitfull$model

    knownid <- which(knowndate != -1)
    N <- data_t$N
    N_data <- data_t$N_data
    initialTitreValue <- data_t$initialTitreValue

    n_post <- post$mcmc[[1]] %>% nrow
    true_inf_time <- res$immune_histories_long %>% filter(value == 1) %>% 
        complete(i = 1:N, fill = list(value = 0, t = -1)) %>% rename(id = i, inf_time = t, inf_ind_sim = value) %>% 
        select(id, inf_ind_sim)

    exp_true <- modelA$simpar$exp %>% as.data.frame %>% 
        mutate(i = 1:N) %>% pivot_longer(!i, names_to = "t", values_to = "exp") %>% mutate(t = as.numeric(substr(t, 2, 4))) %>% 
        filter(exp == 1) %>% complete(i = 1:N, fill = list(t = -1, exp = 0)) %>% 
        rename(id = i, exp_ind_sim = exp, exp_time_sim = t)

    sim_states <- left_join(exp_true, true_inf_time)
    sim_states$init_titre <- data_t$initialTitreValue
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

    post_titreexp %>% filter(id == 4)
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
        modelname = modelname,
        modelname_sim_local = paste0("local/sim_study_inf/", modelname_sim),
        modelname_sim = modelname_sim,
        obs_er = obs_er,
        fit_states = fit_states,
        fit_obs = post_obstitre,
        sim_states = sim_states,
        n_chains = n_chains,
        n_post = n_post)

    if (!priorPred) {
        saveRDS(postprocess, file = here::here("outputs", "fits", "local", "sim_study_inf", modelname_sim, modelname, paste0("pp_obs_", obs_er, ".RDS")))
    } else {
        saveRDS(postprocess, file = here::here("outputs", "fits", "local", "sim_study_inf", modelname_sim, modelname, paste0("pp_prior_obs_", obs_er, ".RDS")))
    }
}

plot_titre_obs_simInf <- function(outputfull, fitfull, fig_folder) {

    filename <- outputfull$filename
    modelname <- outputfull$modelname
    modelname_sim <- outputfull$modelname_sim
    modelname_sim_local <- outputfull$modelname_sim_local
    obs_er <- outputfull$obs_er

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
  #  if (!is.null(scale_ab)) {
   #     p1 <- p1 + scale_y_continuous(breaks = scale_ab %>% as.numeric, labels = scale_ab %>% names)
  #  }
    dir.create(here::here("outputs", "fits", modelname_sim_local, modelname, fig_folder, paste0("obs_", obs_er)), recursive = TRUE)

    ggsave(here::here("outputs", "fits", modelname_sim_local, modelname, fig_folder, paste0("obs_", obs_er),"titre_obs.png"), height = 10, width = 10)

}


plot_abkinetic_simInf <- function(outputfull, fitfull, fig_folder) {
    require(posterior)
    require(bayesplot)
    require(ggdist)


    modelname <- outputfull$modelname
    modelname_sim <- outputfull$modelname_sim
    obs_er <- outputfull$obs_er
    file2 <- paste0("obs_", obs_er)
    modelname_sim_local <- outputfull$modelname_sim_local

    n_chains <- outputfull$n_chains
    n_post <- outputfull$n_post

    chain_samples <- 1:n_chains %>% map(~c(rep(.x, n_post))) %>% unlist

    modelA <- readRDS(here::here("outputs", "sim_data", modelname_sim, "inputs.RDS"))
    res <- readRDS(file = here::here("outputs", "sim_data", modelname_sim, paste0("sim_data_", obs_er, ".rds")))

    data_t <- fitfull$data_t
    post <- fitfull$post
    par_tab <- fitfull$par_tab
    model_outline <- fitfull$model
    T <- modelA$sim$T


    true_ab_par <- data.frame( 
        param = c("a", "b", "c"),
        value = c(modelA$simpar$a, modelA$simpar$b, modelA$simpar$c)
    )

    res$kinetics_parameters %>% filter(name %in% c("a", "b", "c")) %>% select(name, value) %>% rename(param = name)

    compare <- bind_rows(
        post$mcmc %>% lapply(as.data.frame) %>% do.call(rbind, .) %>% as.data.frame %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>%
            mutate(type = "Posterior distribution") %>% filter(param %in% c("a", "b", "c")) ,
        res$kinetics_parameters %>% filter(name %in% c("a", "b", "c")) %>% select(name, value)  %>% rename(param = name) %>%  mutate(type = "Simulated distribution") ,
        purrr::map_df(1:n_post,
            ~model_outline$samplePriorDistributions(par_tab)
        ) %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>%  mutate(type = "Prior distribution")  %>%
        filter(param %in% c("a", "b", "c")) %>% filter(type != "Prior distribution")
    )


    ab_function <- function(a, b, c, T) {
        1:T %>% map( 
            function(t) {
                if (t < 14) {
                    titre_init <- 0 +  (log(exp(a) +  exp(c)) * (t) / 14);
                } else {
                    titre_init <- 0 + (log(exp(a) * exp(-b/10 * ((t) - 14)) + exp(c)));
                }
            }
        ) %>% unlist
    }


    post_fit <- post$mcmc %>% lapply(as.data.frame) %>% do.call(rbind, .) %>% as.data.frame %>% mutate(chain = as.character(chain_samples ))
    a_post <- post_fit[["a"]]
    b_post <- post_fit[["b"]]
    c_post <- post_fit[["c"]]

   # c_post <- c_slope_post #(c_trans_post - exp(- c_slope_post * b_post))

   # p1 <- post_fit %>% ggplot() + geom_point(aes(a, b, color = chain))
   # p2 <- post_fit %>% ggplot() + geom_point(aes(a, c_slope, color = chain))
   # p3 <- post_fit %>% ggplot() + geom_point(aes(b, c_slope, color = chain))

   # (p1 + p2) / (p3 )


    traj_post <- 1:(n_post * n_chains) %>% purrr::map_df(
        ~data.frame(
            time = 1:T,
            value = ab_function(a_post[.x], b_post[.x], c_post[.x], T)
        )
    )
    traj_post_summ <- traj_post %>% group_by(time) %>% mean_qi()

    a_sim <- res$kinetics_parameters %>% filter(name == "a") %>% pull(value)
    b_sim <- res$kinetics_parameters %>% filter(name == "b") %>% pull(value)
    c_sim <- res$kinetics_parameters %>% filter(name == "c") %>% pull(value)

    traj_true <- 1:(length(a_sim)) %>% purrr::map_df(
        ~data.frame(
            time = 1:T,
            value = ab_function(a_sim[.x], b_sim[.x], c_sim[.x], T),
            sim = .x
        )
    )

    true_ab_par <- data.frame( 
        param = c("a", "b", "c"),
        value = c(mean(a_sim), mean(b_sim), mean(c_sim))
    )

    p1 <- compare %>% filter(type == "Posterior distribution") %>% 
        ggplot() + 
            geom_density(aes(x = value, fill = type), alpha = 0.5) +
            geom_vline(data = true_ab_par, aes(xintercept = value), size = 2, alpha = 0.7) +
            scale_fill_manual(values = c("red", "gray")) +
            facet_wrap(vars(param), scale = "free") + theme_bw() + 
            labs(x = "Value", y = "Density", fill = "Type") + theme_bw() + 
            ggtitle("Antibody kinetics parameters")

    p2 <- traj_post_summ %>%  
        ggplot() + 
        geom_line(data = traj_true, aes(y = value, x = time, group = sim, color = "black"),  size = 0.4, alpha = 0.4) + 
        geom_ribbon(aes(ymin = .lower, ymax = .upper, x = time), fill = "red", size = 3, alpha = 0.5) + 
        geom_line(aes(y = value, x = time, color = "red"), size = 2) + 
        theme_bw() + labs(x = expression("Time post-infection, s", y = "Titre boost, f"[ab])) + 
        ggtitle("Antibody trajectories post-infection,") + 
        scale_colour_manual(name = "Line type", 
         values =c('black'='black','red'='red'), labels = c('Simualted trajectory','Posterior trajectory'))

    require(patchwork)

    p1 / p2 + plot_annotation(title = "Simualtion recovery of the antibody kinetics")
    ggsave(here::here("outputs", "fits", modelname_sim_local, modelname, fig_folder, paste0("obs_", obs_er), "ab_kinetics_recov.png"), height = 10, width = 10)

}

plot_inf_rec_simInf <- function(outputfull, fitfull, fig_folder) {

    sim_states <- outputfull$sim_states %>% mutate(
            inf = case_when(
                inf_ind_sim == 1 ~ "Infected",
                inf_ind_sim == 0 ~ "Not infected"
            )
        )
    fit_states <- outputfull$fit_states
    modelname <- outputfull$modelname
    modelname_sim <- outputfull$modelname_sim
    modelname_sim_local <- outputfull$modelname_sim_local

    obs_er <- outputfull$obs_er

    data_t <- fitfull$data_t

    df_data <- data.frame(
        id = 1:200,
        start_titre = data_t$initialTitreValue
    )

    # a) Infection over whole process
    fit_states_inf_prob <- fit_states %>% select(id, inf_ind) %>% summarise(inf_post = mean(inf_ind), .by = "id")

    no_inf_sim <- sim_states %>% select(id, inf_ind_sim) %>% pull(inf_ind_sim) %>% sum
    no_inf_fit_df <- fit_states %>% select(id, sample, inf_ind) %>% summarise(inf_ind_sum = sum(inf_ind), .by = "sample")

    figB <- no_inf_fit_df %>% summarise(n = n()  / 4000, .by = inf_ind_sum) %>%
        ggplot() + 
            geom_col(aes(x = inf_ind_sum, y = n), alpha = 0.8) + 
            theme_bw()  + 
            geom_vline(xintercept = no_inf_sim, size = 3, color = "white") +
            geom_vline(xintercept = no_inf_sim, size = 2, color = "red") +
            labs(y = "Posterior density", x = expression("Estimated number of exposures in epidemic, n"[Z])) + 
            ggtitle("Recovery of population-level infection burden")

    figB + plot_annotation(tag_levels = "A") 
    ggsave(here::here("outputs", "fits", modelname_sim_local, modelname, fig_folder, paste0("obs_", obs_er), "exposure_recov.png"), height = 10, width = 10)

}



plot_exp_times_rec_simInf <- function(outputfull, fitfull, fig_folder) {

    sim_states <- outputfull$sim_states %>% mutate(
            inf = case_when(
                inf_ind_sim == 1 ~ "Infected",
                inf_ind_sim == 0 ~ "Not infected"
            )
        )
    fit_states <- outputfull$fit_states
    modelname <- outputfull$modelname
    modelname_sim_local <- outputfull$modelname_sim_local
    obs_er <- outputfull$obs_er

    data_t <- fitfull$data_t

    df_data <- data.frame(
        id = 1:200,
        start_titre = data_t$initialTitreValue
    )
    fit_states_inf_prob <- fit_states %>% select(id, inf_ind) %>% summarise(inf_post = mean(inf_ind), .by = "id")

    dataplt <- left_join(
        sim_states %>% filter(inf_ind_sim == 1) %>% select(id, exp_time_sim, inf),
        fit_states %>% filter(inf_ind == 1) %>% select(id, inf_time) %>% group_by(id) %>% mean_qi(inf_time)
    ) %>% left_join(fit_states_inf_prob)

    pid_order1 <- dataplt %>% arrange(exp_time_sim) %>% pull(id)

    figB <- dataplt %>% mutate(id = factor(id, levels = pid_order1)) %>% filter(inf == "Infected") %>% 
        ggplot() + 
            geom_linerange(aes(y = id, xmin = .lower, xmax = .upper, color = inf_post)) + 
            geom_point(aes(y = id, x = inf_time)) + 
            geom_point(aes(y = id, x = exp_time_sim), shape = 17, color = "red", size = 1.5, alpha = 0.8) + 
            labs(x = "Time in study", y = "ID", fill = "Posterior probability of exposure") + theme_bw() + 
            facet_grid(cols = vars(inf)) + 
            theme(axis.text.y = element_blank()) + 
            ggtitle("Recovery of timings of infection over epidemic") 
 
    figB + plot_annotation(tag_levels = "A") 
    ggsave(here::here("outputs", "fits", modelname_sim_local, modelname, fig_folder, paste0("obs_", obs_er), "exposure_time_recov.png"), height = 10, width = 10)


    dataplt_error <- dataplt %>% mutate(id = factor(id, levels = pid_order1)) %>% filter(inf == "Infected") %>% 
        mutate(across(inf_time:.upper, ~.x - exp_time_sim)) %>% arrange(inf_time)
    id_figD <- dataplt_error %>% pull(id)
    figD <- dataplt_error %>% mutate(id = factor(id, levels = id_figD)) %>% mutate() %>%
        ggplot() + 
            geom_hline(yintercept = c(14, -14), size = 0.5, color = "gray") +
            geom_linerange(aes(x = id, ymin = .lower, ymax = .upper, color = inf_post), 
                 size = 3) + 
            geom_point(aes(x = id, y = inf_time), size = 3) + 
            labs(y = "Difference between model-predicted infection day \n time and simulated infection day", x = "ID", color = "Posterior probability of exposure") + theme_bw() + 
            theme(axis.text.x = element_blank(), legend.position = "top") + 
            ggtitle("Model error in recovering exposure itime of infected individuals") 
    figD
    ggsave(here::here("outputs", "fits", modelname_sim_local, modelname, fig_folder, paste0("obs_", obs_er), "exposure_time_recov_error.png"), height = 10, width = 10)

    vector_need <- c((dataplt_error %>% filter(abs(inf_time) <= 14) %>% nrow) / (dataplt_error %>% nrow),
        dataplt_error %>% filter(abs(inf_time) <= 14) %>% nrow,
        dataplt_error %>% nrow)
    write.csv(vector_need, file = here::here("outputs", "fits", modelname_sim_local, modelname, fig_folder, paste0("obs_", obs_er), "exp_info.csv"))
}



plot_cop_rec_simInf <- function(outputfull, fitfull, fig_folder) {

    outputfull
    fit_states <- outputfull$fit_states
    filename <- outputfull$filename
    modelname <- outputfull$modelname
    obs_er <- outputfull$obs_er
    modelname_sim_local <- outputfull$modelname_sim_local

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


    biomarker <- "IgG"
    df_data <- data.frame(
        id = 1:data_t$N,
        start_titre = data_t$initialTitreValue %>% as.data.frame %>% pull(!!biomarker),
        known_inf = data_t$knownInfsVec
    ) %>% mutate(known_inf = ifelse(known_inf == 0, "Not known", "Known"))

    titre_cop_sum <- fit_states %>% filter(!!sym(biomarker) != -1) %>%
            group_by(id) %>% mean_qi(!!sym(biomarker) ) %>% arrange(!!biomarker) %>%
            rename(titre_val = !!(biomarker))
    u_ids <- titre_cop_sum$id

    df_data_post <- fit_states  %>% select(id, inf_ind) %>%
                summarise(inf_post = mean(inf_ind), n = n(), .by = "id") %>%
            left_join(df_data) %>% mutate(prop = n / n_length)

    cop_exp_sum_plot <- titre_cop_sum %>% left_join(df_data_post) %>% mutate(id = factor(id, levels = u_ids)) %>%
        filter(!is.na(n)) %>% mutate(biomarker = !!biomarker)

    cop_exp_sum_plot %>% ggplot() + 
        geom_smooth(aes(x = titre_val, y = inf_post)) +
        labs(x = "Initial titre", y = "Posterior probability of infection", color = "Known infection status") + 
         geom_linerange(data = cop_exp_sum_plot,
            aes(y = inf_post, xmin = .lower, xmax = .upper), size = 1) + 
        geom_point(data = cop_exp_sum_plot, aes(x = titre_val, y = inf_post, color = known_inf))
    ggsave(here::here("outputs", "fits", modelname_sim_local, modelname, fig_folder, paste0("obs_", obs_er), "cop_recov.png"), height = 10, width = 10)

  #  if (!is.null(scale_ab)) {
   #     figA <- figA + scale_x_continuous(breaks = scale_ab %>% as.numeric, labels = scale_ab %>% names)
   # }


}