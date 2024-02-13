postprocess_run <- function(modelname, modelname_sim, obs_er, n_chains) {

    library(ggdist)
    dir.create(here::here("outputs", "sim_data", modelname), showWarnings = FALSE)
    dir.create(here::here("outputs", "sim_data", modelname, "figs"), showWarnings = FALSE)
    dir.create(here::here("outputs", "sim_data", modelname, "figs", paste0("obs_", obs_er)), showWarnings = FALSE)
    modelA <- readRDS(here::here("outputs", "sim_data", modelname_sim, "inputs.RDS"))
    res <- readRDS(file = here::here("outputs", "sim_data", modelname_sim, paste0("sim_data_", obs_er, ".rds")))
    fitfull <- readRDS(here::here("outputs", "fits", modelname_sim, modelname, paste0("fit_", "obs_", obs_er, ".RDS")))

    post <- fitfull$post
    data_t <- fitfull$data_t

    knowndate <- data_t$knownInfsDate

    knownid <- which(knowndate != -1)
    N <- data_t$N
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
    post_exp_combine <- post$jump %>% combine
    post_inf_combine <- post$inf %>% combine

    fit_states <- map_df(1:N, 
        function(x) {
            data.frame(
                id = x,
                sample = 1:(n_post * n_chains),
                chain_no = c(rep(1, n_post), rep(2, n_post), rep(3, n_post), rep(4, n_post)),
                exp_ind = as.numeric(post_exp_combine[, x] > 0),
                exp_time = post_exp_combine[, x],
                inf_ind = post_inf_combine[, x]
            )
        }
    )
    postprocess <- list(
        modelname = modelname,
        modelname_sim = modelname_sim,
        obs_er = obs_er,
        fit_states = fit_states,
        sim_states = sim_states,
        n_chains = n_chains,
        n_post = n_post)

    saveRDS(postprocess, file = here::here("outputs", "fits", modelname_sim, modelname, paste0("pp_obs_", obs_er, ".RDS")))

}


plot_trace <- function(outputfull) {
    require(posterior)
    require(bayesplot)
    require(ggdist)

    modelname_sim <- outputfull$modelname_sim
    modelname <- outputfull$modelname
    obs_er <- outputfull$obs_er

    fitfull <- readRDS(here::here("outputs", "fits", modelname_sim, modelname, paste0("fit_", "obs_", obs_er, ".RDS")))
    post <- fitfull$post

    post$mcmc %>% combine %>% summary
    p1 <- post$mcmc %>% mcmc_trace
    p2 <- post$lpost %>% ggplot() + geom_line(aes(x = sample_no, y = lpost, color = chain_no))
    p1 / p2

    dir.create(here::here("outputs", "fits", modelname_sim, modelname, "figs"))
    dir.create(here::here("outputs", "fits", modelname_sim, modelname, "figs", paste0("obs_", obs_er)))

    ggsave(here::here("outputs", "fits", modelname_sim, modelname, "figs", paste0("obs_", obs_er), "trace_plots.png"))
}


plot_abkinetics <- function(outputfull) {
    require(posterior)
    require(bayesplot)
    require(ggdist)


    modelname <- outputfull$modelname
    modelname_sim <- outputfull$modelname_sim
    obs_er <- outputfull$obs_er

    n_chains <- outputfull$n_chains
    n_post <- outputfull$n_post

    chain_samples <- 1:n_chains %>% map(~c(rep(.x, n_post))) %>% unlist

    modelA <- readRDS(here::here("outputs", "sim_data", modelname_sim, "inputs.RDS"))
    res <- readRDS(file = here::here("outputs", "sim_data", modelname_sim, paste0("sim_data_", obs_er, ".rds")))
    fitfull <- readRDS(here::here("outputs", "fits", modelname_sim, modelname, paste0("fit_", "obs_", obs_er, ".RDS")))
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
        post$mcmc %>% combine %>% as.data.frame %>% mutate(c = c_slope) %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>%
            mutate(type = "Posterior distribution") %>% filter(param %in% c("a", "b", "c")) ,
        res$kinetics_parameters %>% filter(name %in% c("a", "b", "c")) %>% select(name, value)  %>% rename(param = name) %>%  mutate(type = "Simulated distribution") ,
        map_df(1:n_post,
            ~model_outline$samplePriorDistributions(par_tab)
        ) %>% rename(c = c_slope) %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>%  mutate(type = "Prior distribution")  %>%
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


    post_fit <- post$mcmc %>% combine %>% as.data.frame %>% mutate(chain = as.character(chain_samples ))
    a_post <- post_fit[["a"]]
    b_post <- post_fit[["b"]]
    c_slope_post <- post_fit[["c_slope"]]

    c_post <- c_slope_post #(c_trans_post - exp(- c_slope_post * b_post))

    p1 <- post_fit %>% ggplot() + geom_point(aes(a, b, color = chain))
    p2 <- post_fit %>% ggplot() + geom_point(aes(a, c_slope, color = chain))
    p3 <- post_fit %>% ggplot() + geom_point(aes(b, c_slope, color = chain))

    (p1 + p2) / (p3 )


    traj_post <- 1:(n_post * n_chains) %>% map_df(
        ~data.frame(
            time = 1:T,
            value = ab_function(a_post[.x], b_post[.x], c_post[.x], T)
        )
    )
    traj_post_summ <- traj_post %>% group_by(time) %>% mean_qi()

    a_sim <- res$kinetics_parameters %>% filter(name == "a") %>% pull(value)
    b_sim <- res$kinetics_parameters %>% filter(name == "b") %>% pull(value)
    c_sim <- res$kinetics_parameters %>% filter(name == "c") %>% pull(value)

    traj_true <- 1:(length(a_sim)) %>% map_df(
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
    ggsave(here::here("outputs", "fits",  modelname_sim, modelname, "figs", paste0("obs_", obs_er), "ab_kinetics_recov.png"), height = 10, width = 10)

}

plot_exp_rec <- function(outputfull) {

    sim_states <- outputfull$sim_states %>% mutate(
            expinf = case_when(
                exp_ind_sim == 1 & inf_ind_sim == 1 ~ "Exposed and infected",
                exp_ind_sim == 1 & inf_ind_sim == 0 ~ "Exposed and not infected",
                exp_ind_sim == 0 ~ "Not exposed"
            )
        )
    fit_states <- outputfull$fit_states
    modelname <- outputfull$modelname
    modelname_sim <- outputfull$modelname_sim
    obs_er <- outputfull$obs_er

    fitfull <- readRDS(here::here("outputs", "fits", modelname_sim, modelname, paste0("fit_", "obs_", obs_er, ".RDS")))
    data_t <- fitfull$data_t

    df_data <- data.frame(
        id = 1:200,
        start_titre = data_t$initialTitreValue
    )

    # a) Infection over whole process
    fit_states_exp_prob <- fit_states %>% select(id, exp_ind) %>% summarise(exp_post = mean(exp_ind), .by = "id")

    figA <- left_join(
        sim_states %>% select(id, expinf),
        fit_states_exp_prob
    ) %>% left_join(df_data) %>% 
        ggplot() + 
            geom_point(aes(x = start_titre, y = exp_post, fill = expinf), shape = 21, size = 3, alpha = 0.6) + 
            scale_fill_manual(values = c("red", "blue", "green")) + 
            theme_bw() + 
            ylim(0, 1) +
            labs(x = "Titre at start of study", y = "Posterior probability of exposure", fill = "Exposure and infection status") + 
            ggtitle("Posterior plots of probability of exposure compared to simulated data") 

    no_exp_sim <- sim_states %>% select(id, exp_ind_sim) %>% pull(exp_ind_sim) %>% sum
    no_exp_fit_df <- fit_states %>% select(id, sample, exp_ind) %>% summarise(exp_ind_sum = sum(exp_ind), .by = "sample")

    figB <- no_exp_fit_df %>% summarise(n = n()  / 4000, .by = exp_ind_sum) %>%
        ggplot() + 
            geom_col(aes(x = exp_ind_sum, y = n), alpha = 0.8) + 
            theme_bw()  + 
            geom_vline(xintercept = no_exp_sim, size = 3, color = "white") +
            geom_vline(xintercept = no_exp_sim, size = 2, color = "red") +
            labs(y = "Posterior density", x = expression("Estimated number of exposures in epidemic, n"[Z])) + 
            ggtitle("Recovery of population-level exposure burden")

    figB + plot_annotation(tag_levels = "A") 
    ggsave(here::here("outputs", "fits", modelname_sim, modelname, "figs", paste0("obs_", obs_er), "exposure_recov.png"), height = 10, width = 10)

}

plot_exp_times_rec <- function(outputfull) {

    sim_states <- outputfull$sim_states %>% mutate(
            expinf = case_when(
                exp_ind_sim == 1 & inf_ind_sim == 1 ~ "Exposed and infected",
                exp_ind_sim == 1 & inf_ind_sim == 0 ~ "Exposed and not infected",
                exp_ind_sim == 0 ~ "Not exposed"                
            )
        )
    fit_states <- outputfull$fit_states
    modelname <- outputfull$modelname
    modelname_sim <- outputfull$modelname_sim
    obs_er <- outputfull$obs_er

    fitfull <- readRDS(here::here("outputs", "fits", modelname_sim, modelname, paste0("fit_", "obs_", obs_er, ".RDS")))
    data_t <- fitfull$data_t

    df_data <- data.frame(
        id = 1:200,
        start_titre = data_t$initialTitreValue
    )
    fit_states_exp_prob <- fit_states %>% select(id, exp_ind) %>% summarise(exp_post = mean(exp_ind), .by = "id")

    dataplt <- left_join(
        sim_states %>% filter(exp_ind_sim == 1) %>% select(id, exp_time_sim, expinf),
        fit_states %>% filter(exp_ind == 1) %>% select(id, exp_time) %>% group_by(id) %>% mean_qi(exp_time)
    ) %>% left_join(fit_states_exp_prob)

    pid_order1 <- dataplt %>% arrange(exp_time_sim) %>% pull(id)

    figB <- dataplt %>% mutate(id = factor(id, levels = pid_order1)) %>% filter(expinf == "Exposed and infected") %>% 
        ggplot() + 
            geom_linerange(aes(y = id, xmin = .lower, xmax = .upper, color = exp_post)) + 
            geom_point(aes(y = id, x = exp_time)) + 
            geom_point(aes(y = id, x = exp_time_sim), shape = 17, color = "red", size = 1.5, alpha = 0.8) + 
            labs(x = "Time in study", y = "ID", fill = "Posterior probability of exposure") + theme_bw() + 
            facet_grid(cols = vars(expinf)) + 
            theme(axis.text.y = element_blank()) + 
            ggtitle("Recovery of timings of exposure over epidemic") 


    dataplt_inf <- left_join(
        sim_states %>% filter(inf_ind_sim == 1) %>% select(id, exp_time_sim),
        fit_states %>% filter(inf_ind == 1) %>% select(id, exp_time) %>% group_by(id) %>% mean_qi(exp_time)
    )

    figC <- dataplt_inf %>% select(id, exp_time_sim, exp_time) %>% 
        rename(Data = exp_time_sim, `Model recovery` = exp_time) %>%
        pivot_longer(!id, names_to = "type", values_to = "time") %>%
        ggplot() + 
            geom_density(aes(x = time, y = ..count.., fill = type), size = 2, alpha = 0.5, shape = 21) + 
            theme_bw()  + 
            xlim(0, 120) + 
            labs(x = "Time in study", y = "Number of people infected daily", fill = "") + 
            ggtitle("Recovery of infection timings")

 
    figB / figC + plot_annotation(tag_levels = "A") 
    ggsave(here::here("outputs", "fits", modelname_sim, modelname, "figs", paste0("obs_", obs_er), "exposure_time_recov.png"), height = 10, width = 10)


    dataplt_error <- dataplt %>% mutate(id = factor(id, levels = pid_order1)) %>% filter(expinf == "Exposed and infected") %>% 
        mutate(across(exp_time:.upper, ~.x - exp_time_sim)) %>% arrange(exp_time)
    id_figD <- dataplt_error %>% pull(id)
    figD <- dataplt_error %>% mutate(id = factor(id, levels = id_figD)) %>% mutate() %>%
        ggplot() + 
            geom_hline(yintercept = c(14, -14), size = 0.5, color = "gray") +
            geom_linerange(aes(x = id, ymin = .lower, ymax = .upper, color = exp_post), 
                 size = 3) + 
            geom_point(aes(x = id, y = exp_time), size = 3) + 
            labs(y = "Difference between model-predicted infection day \n time and simulated infection day", x = "ID", color = "Posterior probability of exposure") + theme_bw() + 
            theme(axis.text.x = element_blank(), legend.position = "top") + 
            ggtitle("Model error in recovering exposure itime of infected individuals") 
    figD
    ggsave(here::here("outputs", "fits", modelname_sim, modelname, "figs", paste0("obs_", obs_er), "exposure_time_recov_error.png"), height = 10, width = 10)

    vector_need <- c((dataplt_error %>% filter(abs(exp_time) <= 14) %>% nrow) / (dataplt_error %>% nrow),
        dataplt_error %>% filter(abs(exp_time) <= 14) %>% nrow,
        dataplt_error %>% nrow)
    write.csv(vector_need, file = here::here("outputs", "fits", modelname_sim, modelname, "figs", paste0("obs_", obs_er), "exp_info.csv"))
}


plot_inf_rec <- function(outputfull) {
    sim_states <- outputfull$sim_states
    fit_states <- outputfull$fit_states
    modelname <- outputfull$modelname
    modelname_sim <- outputfull$modelname_sim
    obs_er <- outputfull$obs_er

    fitfull <- readRDS(here::here("outputs", "fits", modelname_sim, modelname, paste0("fit_", "obs_", obs_er, ".RDS")))
    data_t <- fitfull$data_t

    df_data <- data.frame(
        id = 1:200,
        start_titre = data_t$initialTitreValue
    )



    # b) Infection given exposure 

    figB <- left_join(
        sim_states %>% filter(exp_ind_sim == 1) %>% select(id, inf_ind_sim)  %>% mutate(inf_status_data = recode(inf_ind_sim, "1" = "Infected", "0" = "Not infected")),
        fit_states %>% filter(exp_ind == 1) %>% select(id, inf_ind) %>% summarise(inf_post = mean(inf_ind), .by = "id")
    ) %>% left_join(df_data) %>% 
        ggplot() + 
            geom_point(aes(x = start_titre, y = inf_post, fill = inf_status_data), size = 2, alpha = 0.7, shape = 21) + 
            theme_bw() +
            labs(x = expression("Titre at start of study Y"[j]^0), y = expression("Expectation of posterior distribution of infection E[Z"[j]*"]"), fill = "Simulated infection status") + 
            ggtitle("Recovery of individual-level infection status")


    no_inf_sim <- sim_states %>% filter(exp_ind_sim == 1) %>% select(id, inf_ind_sim) %>% pull(inf_ind_sim) %>% sum
    no_inf_fit_df <- fit_states %>% filter(exp_ind == 1) %>% select(id, sample, inf_ind) %>% summarise(inf_ind_sum = sum(inf_ind), .by = "sample")

    figC <- no_inf_fit_df %>% summarise(n = n()  / 4000, .by = inf_ind_sum) %>%
        ggplot() + 
            geom_col(aes(x = inf_ind_sum, y = n), alpha = 0.8) + 
            theme_bw()  + 
            geom_vline(xintercept = no_inf_sim, size = 3, color = "white") +
            geom_vline(xintercept = no_inf_sim, size = 2, color = "red") +
            labs(y = "Posterior density", x = expression("Estimated number of infections in epidemic, n"[Z])) + 
            ggtitle("Recovery of population-level infection burden")

    dataplt_inf <- left_join(
        sim_states %>% filter(inf_ind_sim == 1) %>% select(id, exp_time_sim),
        fit_states %>% filter(inf_ind == 1) %>% select(id, exp_time) %>% group_by(id) %>% mean_qi(exp_time)
    )

    #figD <- dataplt_inf %>% select(id, exp_time_sim, exp_time) %>% 
    #    rename(Data = exp_time_sim, `Model recovery` = exp_time) %>%
    #    pivot_longer(!id, names_to = "type", values_to = "time") %>%
   # ggplot() + 
    #        geom_density(aes(x = time, fill = type), size = 2, alpha = 0.5, shape = 21) + 
     #       theme_bw()  + 
      #      labs(x = "Time in study", y = "Density", fill = "") + 
       #     ggtitle("Recovery of infection timings")

    (figB) / figC + plot_annotation(tag_levels = "A")
    ggsave(here::here("outputs", "fits", modelname_sim, modelname, "figs", paste0("obs_", obs_er), "infection_recov.png"), height = 10, width = 10)
}


plot_cop_rec <- function(outputfull) {

    modelname <- outputfull$modelname
    modelname_sim <- outputfull$modelname_sim
    obs_er <- outputfull$obs_er
    fit_states <- outputfull$fit_states

    sim_states <- outputfull$sim_states
    sim_df <- sim_states %>% filter(exp_ind_sim == 1) %>% select(id, inf_ind_sim)  %>% mutate(inf_status_data = recode(inf_ind_sim, "1" = "Infected", "0" = "Not infected"))

    fitfull <- readRDS(here::here("outputs", "fits", modelname_sim, modelname, paste0("fit_", "obs_", obs_er, ".RDS")))
    post <- fitfull$post

    reps <-  post$mcmc %>% combine %>% as.data.frame %>% mutate(c = c_slope) %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>%
            mutate(type = "Posterior distribution") %>% filter(param %in% c("beta0", "beta1"))

    b0_rep <- reps %>% filter(param == "beta0") %>% pull(value)
    b1_rep <- reps %>% filter(param == "beta1") %>% pull(value)

    biomarker_protection <- function(biomarker_quantity, biomarker_prot_midpoint, biomarker_prot_width) {
        risk <- 1 - 1/(1 + exp(biomarker_prot_width * (biomarker_quantity - biomarker_prot_midpoint)))
        return(risk)
    }
    modelA <- readRDS(here::here("outputs", "sim_data", modelname_sim, "inputs.RDS"))


    require(tidybayes)
    df_here <- data.frame(
        titre = seq(0, 4, 0.1),
        true_data = (1 - biomarker_protection(seq(0, 4, 0.1), 2, modelA$protection_curve))
    )

    lol <- 1:1000 %>% map_df(
        ~data.frame(
            t = seq(0, 4, 0.1),
            vals = 1.0 / (1.0 + exp(-(b0_rep[.x] + b1_rep[.x] * seq(0, 4, 0.1))))
        )
    )

    lol %>% ggplot() + 
        geom_smooth(data = sim_states %>% filter(exp_ind_sim == 1), aes(x = init_titre, y = inf_ind_sim), color = "white", se = FALSE, method = "glm",  formula = y ~ x, method.args = list(family = "binomial"), size = 5, alpha = 0.3) +
        geom_smooth(data = sim_states %>% filter(exp_ind_sim == 1), aes(x = init_titre, y = inf_ind_sim, color = "Simulated COP curve"), se = FALSE, method = "glm",  formula = y ~ x, method.args = list(family = "binomial"), size = 5, alpha = 0.3) +
        stat_lineribbon(aes(x = t, y = vals, color = "Recovered COP curve"), fill = "blue", alpha = 0.5, .width = c(0.5, 0.95)) + 
        theme_bw() + 
        ylim(0, 1) +
        scale_color_manual(values = c("red", "black")) +
        labs(x = expression("Titre at start of study Y"[j]^0), y = expression("Posterior probability of infection for correlate of protection, f"[cop]*"(Y"[j]^0*", "*theta[cop]*")"), color = "Curve type")
    ggsave(here::here("outputs", "fits", modelname_sim, modelname, "figs", paste0("obs_", obs_er), "cop_recov.png"), height = 10, width = 10)

}


postprocessFigs <- function(modelname, modelname_sim, obs_er, n_chains) {
   # modelname <- "inferExp"
   # modelname_sim <-  "cesNoCOP_notd"
    #obs_er <- "0.3"
   # n_chains <- 4

    postprocess_run( modelname, modelname_sim, obs_er, n_chains)
    outputfull <- readRDS(file = here::here("outputs", "fits", modelname_sim, modelname, paste0("pp_obs_", obs_er, ".RDS")))

    plot_trace(outputfull)
    plot_abkinetics(outputfull)
    plot_inf_rec(outputfull)
    plot_exp_rec(outputfull)
    plot_exp_times_rec(outputfull)
    plot_cop_rec(outputfull)
}

