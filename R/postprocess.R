#' @useDynLib rjmc
#' @importFrom Rcpp sourceCpp
#' @import ggplot2
#' @import purrr
#' @import tidyr
#' @import dplyr
#' @import ggdist
#' @import patchwork
#' @import tidybayes
#' @importFrom magrittr %>% %<>%
NULL
postprocess_run <- function(filename, modelname, n_chains) {

    library(ggdist)
    fitfull <- readRDS(here::here("outputs", "fits", filename, paste0("fit_", modelname, ".RDS")))

    post <- fitfull$post
    data_t <- fitfull$data_t

    knowndate <- data_t$knownInfsDate

    knownid <- which(knowndate != -1)
    N <- data_t$N
    initialTitreValue <- data_t$initialTitreValue

    n_post <- post$mcmc[[1]] %>% nrow
    post_exp_combine <- post$jump %>% combine
    post_inf_combine <- post$inf %>% combine

    fit_states <- purrr::map_df(1:N, 
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
        filename = filename,
        modelname = modelname,
        fit_states = fit_states,
        n_chains = n_chains,
        n_post = n_post)

    saveRDS(postprocess, file = here::here("outputs", "fits", filename, paste0("pp_", modelname, ".RDS")))

}



plot_trace <- function(outputfull) {
    require(posterior)
    require(bayesplot)
    require(ggdist)

    filename <- outputfull$filename
    modelname <- outputfull$modelname

    fitfull <- readRDS(here::here("outputs", "fits", filename, paste0("fit_", modelname, ".RDS")))
    post <- fitfull$post

    post$mcmc %>% combine %>% summary
    p1 <- post$mcmc %>% mcmc_trace
    p2 <- post$lpost %>% ggplot() + geom_line(aes(x = sample_no, y = lpost, color = chain_no))
    p1 / p2

    ggsave(here::here("outputs", "fits", filename, "figs", modelname, "trace_plots.png"))
}


plot_abkinetics <- function(outputfull) {
    require(posterior)
    require(bayesplot)
    require(ggdist)

    filename <- outputfull$filename
    modelname <- outputfull$modelname
    obs_er <- outputfull$obs_er

    n_chains <- outputfull$n_chains
    n_post <- outputfull$n_post

    chain_samples <- 1:n_chains %>% map(~c(rep(.x, n_post))) %>% unlist

    fitfull <- readRDS(here::here("outputs", "fits", filename, paste0("fit_", modelname, ".RDS")))
    data_t <- fitfull$data_t
    N <- data_t$N

    post <- fitfull$post
    par_tab <- fitfull$par_tab
    model_outline <- fitfull$model


    compare <- bind_rows(
        post$mcmc %>% combine %>% as.data.frame  %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>%
             mutate(type = "Posterior distribution") %>% filter(param %in% c("a", "b", "c")),
        purrr::map_df(1:n_post,
            ~model_outline$samplePriorDistributions(par_tab)
        )  %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>%  mutate(type = "Prior distribution")  %>%
        filter(param %in% c("a", "b", "c"))
    )


    ab_function <- function(a, b, c, T) {
        1:T %>% map( 
            function(t) {
                if (t < 14) {
                    titre_init <- 0 + (log(exp(a) +  exp(c)) * (t) / 14);
                } else {
                    titre_init <- 0 + (log(exp(a) * exp(-b/10 * ((t) - 14)) + exp(c)));
                }
            }
        ) %>% unlist
    }


    post_fit <- post$mcmc %>% combine %>% as.data.frame %>% mutate(chain = as.character(chain_samples ))
    a_post <- post_fit[["a"]]
    b_post <- post_fit[["b"]]
    c_post <- post_fit[["c"]]


    T <- 300
    traj_post <- 1:(n_post * n_chains) %>% purrr::map_df(
        ~data.frame(
            time = 1:T,
            value = ab_function(a_post[.x], b_post[.x], c_post[.x], T)
        )
    )
    traj_post_summ <- traj_post %>% group_by(time) %>% mean_qi()

    p1 <- compare %>% filter(type == "Posterior distribution") %>% 
        ggplot() + 
            geom_density(aes(x = value, fill = type), alpha = 0.5) +
            scale_fill_manual(values = c("red", "gray")) +
            facet_wrap(vars(param), scale = "free") + theme_bw() + 
            labs(x = "Value", y = "Density", fill = "Type") + theme_bw() + 
            ggtitle("Antibody kinetics parameters")
    p2 <- traj_post_summ %>%  
        ggplot() + 
        geom_ribbon(aes(ymin = .lower, ymax = .upper, x = time), fill = "red", size = 3, alpha = 0.5) + 
        geom_line(aes(y = value, x = time, color = "red"), size = 2) + 
        theme_bw() + labs(x = expression("Time post-infection, s", y = "Titre boost, f"[ab])) + 
        ggtitle("Antibody trajectories post-infection,") + 
        scale_colour_manual(name = "Line type", 
         values =c('black'='black','red'='red'), labels = c('Simualted trajectory','Posterior trajectory'))

    require(patchwork)

    p1 / p2 + plot_annotation(title = "Simualtion recovery of the antibody kinetics")
    ggsave(here::here("outputs", "fits", filename, "figs", modelname, "ab_kinetics_recov.png"), height = 10, width = 10)

}

plot_exp_rec <- function(outputfull) {


    fit_states <- outputfull$fit_states
    filename <- outputfull$filename
    modelname <- outputfull$modelname

    fitfull <- readRDS(here::here("outputs", "fits", filename, paste0("fit_", modelname, ".RDS")))
    data_t <- fitfull$data_t
    N <- data_t$N

    df_data <- data.frame(
        id = 1:N,
        start_titre = data_t$initialTitreValue
    )

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

    figB <- no_exp_fit_df %>% summarise(n = n()  / 4000, .by = exp_ind_sum) %>%
        ggplot() + 
            geom_col(aes(x = exp_ind_sum, y = n), alpha = 0.8) + 
            theme_bw()  + 
            labs(y = "Posterior density", x = expression("Estimated number of exposures in epidemic, n"[Z])) + 
            ggtitle("Recovery of population-level exposure burden")

    figB + plot_annotation(tag_levels = "A") 
    ggsave(here::here("outputs", "fits", filename, "figs", modelname, "exposure_recov.png"), height = 10, width = 10)

}


plot_exp_times_rec <- function(outputfull) {

    fit_states <- outputfull$fit_states
    filename <- outputfull$filename
    modelname <- outputfull$modelname

    fitfull <- readRDS(here::here("outputs", "fits", filename, paste0("fit_", modelname, ".RDS")))
    data_t <- fitfull$data_t
    N <- data_t$N

    df_data <- data.frame(
        id = as.character(1:N),
        start_titre = data_t$initialTitreValue,
        known_inf = data_t$knownInfsVec,
        known_date =  data_t$knownInfsDate
    ) %>% mutate(known_inf = ifelse(known_inf == 0, "Not known", "Known"))

    fit_states_exp_prob <- fit_states %>% select(id, exp_ind) %>% summarise(exp_post = mean(exp_ind), .by = "id")

    dataplt <- fit_states %>% filter(exp_ind == 1) %>% select(id, exp_time) %>% group_by(id) %>% mean_qi(exp_time) %>% left_join(fit_states_exp_prob)

    pid_order1 <- dataplt %>% arrange(exp_time) %>% pull(id)

    figB <- dataplt %>% mutate(id = factor(id, levels = pid_order1)) %>%
        left_join(df_data) %>% mutate(id = factor(id, levels = pid_order1)) %>%
        ggplot() + 
            geom_linerange(aes(y = id, xmin = .lower, xmax = .upper, color = exp_post)) + 
            geom_point(aes(y = id, x = exp_time, shape = known_inf, fill = known_inf), color = "black") + 
            scale_fill_manual(values = c("red", "black")) + 
            scale_shape_manual(values = c(21, 22)) + 
            guides(shape = guide_legend(title = "Known infection status"),
                fill = guide_legend(title = "Known infection status")) +
            labs(x = "Time in study", y = "ID", fill = "Posterior probability of exposure", color = "Probability of exposure") + theme_bw() + 
            theme(axis.text.y = element_blank()) + 
            ggtitle("Recovery of timings of exposure over epidemic") 

    
    datasets <- fit_states %>% filter(exp_time > -1, inf_ind == 1)  %>% arrange(sample) %>% 
        split(cumsum(c(TRUE, diff(.$sample) != 0))) %>% map(~pull(.x, exp_time) )
    minmax <- datasets %>% unlist() %>% range
    break_vec <- ceiling(seq(0, minmax[2], length.out = 30))

    hist_data <- datasets %>%
        map(~ hist(., breaks = break_vec, plot = FALSE)) %>%
        map(~ {
            data_frame(
            count = .$counts,
            bin = .$mids
            )
        })
    mean_counts <- hist_data %>%
        map(~ .x %>% pull(count)) %>%  do.call(rbind, .) %>% colMeans
    std_counts <- do.call(rbind, hist_data %>% map(~ .$count)) %>% apply(2, sd)
    
    dataplt_inf <- fit_states %>% filter(exp_time > -1, inf_ind == 1) %>% select(id, exp_time) %>% group_by(id) %>% mean_qi(exp_time)

    df_data_t <- df_data %>% filter(known_inf == "Known") %>% select(id, type = known_inf, time = known_date ) %>% 
        mutate(id = as.numeric(id))
# Plotting
    figC <- df_data_t %>% ggplot() +
      geom_histogram(aes(x = time, y = ..count.., fill = "gray40"), size = 2, alpha = 0.5, color = "black") + 
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
 
    figB / figC + plot_annotation(tag_levels = "A") 
    ggsave(here::here("outputs", "fits", filename, "figs", modelname, "exposure_time_recov.png"), height = 10, width = 10)
}


plot_inf_rec <- function(outputfull) {
    fit_states <- outputfull$fit_states
    filename <- outputfull$filename
    modelname <- outputfull$modelname


    fitfull <- readRDS(here::here("outputs", "fits", filename, paste0("fit_", modelname, ".RDS")))
    data_t <- fitfull$data_t
    N <- data_t$N

    df_data <- data.frame(
        id = 1:N,
        start_titre = data_t$initialTitreValue,
        known_inf = data_t$knownInfsVec
    ) %>% mutate(known_inf = ifelse(known_inf == 0, "Not known", "Known"))

    # b) Infection given exposure 

    figB <- fit_states %>% filter(exp_ind == 1) %>% select(id, inf_ind) %>% summarise(inf_post = mean(inf_ind), .by = "id") %>%
            left_join(df_data) %>% 
        ggplot() + 
            geom_count(aes(x = start_titre, y = inf_post, shape = known_inf, color = known_inf), size = 2, alpha = 0.7) + 
            theme_bw() +
            labs(x = expression("Titre at start of study Y"[j]^0),
                y = expression("Expectation of posterior distribution of infection E[Z"[j]*"]"), shape = "", color = "") + 
            ggtitle("Recovery of individual-level infection status")

    no_inf_fit_df <- fit_states %>% filter(exp_ind == 1) %>% select(id, sample, inf_ind) %>% summarise(inf_ind_sum = sum(inf_ind), .by = "sample")

    figC <- no_inf_fit_df %>% summarise(n = n()  / 4000, .by = inf_ind_sum) %>%
        ggplot() + 
            geom_col(aes(x = inf_ind_sum, y = n), alpha = 0.8) + 
            theme_bw()  + 
            labs(y = "Posterior density", x = expression("Estimated number of infections in epidemic, n"[Z])) + 
            ggtitle("Recovery of population-level infection burden")

    dataplt_inf <- 
        fit_states %>% filter(inf_ind == 1) %>% select(id, exp_time) %>% group_by(id) %>% mean_qi(exp_time)

    #figD <- dataplt_inf %>% select(id, exp_time_sim, exp_time) %>% 
    #    rename(Data = exp_time_sim, `Model recovery` = exp_time) %>%
    #    pivot_longer(!id, names_to = "type", values_to = "time") %>%
   # ggplot() + 
    #        geom_density(aes(x = time, fill = type), size = 2, alpha = 0.5, shape = 21) + 
     #       theme_bw()  + 
      #      labs(x = "Time in study", y = "Density", fill = "") + 
       #     ggtitle("Recovery of infection timings")

    (figB) / figC + plot_annotation(tag_levels = "A")
    ggsave(here::here("outputs", "fits", filename, "figs", modelname, "infection_recov.png"), height = 10, width = 10)
}


plot_cop_rec <- function(outputfull) {

    fit_states <- outputfull$fit_states
    filename <- outputfull$filename
    modelname <- outputfull$modelname

    fitfull <- readRDS(here::here("outputs", "fits", filename, paste0("fit_", modelname, ".RDS")))
    post <- fitfull$post

    reps <-  post$mcmc %>% combine %>% as.data.frame  %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>%
            mutate(type = "Posterior distribution") %>% filter(param %in% c("beta0", "beta1"))

    b0_rep <- reps %>% filter(param == "beta0") %>% pull(value)
    b1_rep <- reps %>% filter(param == "beta1") %>% pull(value)

    biomarker_protection <- function(biomarker_quantity, biomarker_prot_midpoint, biomarker_prot_width) {
        risk <- 1 - 1/(1 + exp(biomarker_prot_width * (biomarker_quantity - biomarker_prot_midpoint)))
        return(risk)
    }

    lol <- 1:1000 %>% purrr::map_df(
        ~data.frame(
            t = seq(0, 4, 0.1),
            vals = 1.0 / (1.0 + exp(-(b0_rep[.x] + b1_rep[.x] * seq(0, 4, 0.1))))
        )
    )

    lol %>% ggplot() + 
        stat_lineribbon(aes(x = t, y = vals, color = "Recovered COP curve"), fill = "blue", alpha = 0.5, .width = c(0.5, 0.95)) + 
        theme_bw() + 
        ylim(0, 1) +
        labs(x = expression("Titre at start of study Y"[j]^0), y = expression("Posterior probability of infection for correlate of protection, f"[cop]*"(Y"[j]^0*", "*theta[cop]*")"), color = "Curve type")
    ggsave(here::here("outputs", "fits", filename, "figs", modelname, "cop_recov.png"), height = 10, width = 10)

}


plot_abkinetics_delta <- function(outputfull) {
    require(posterior)
    require(bayesplot)
    require(ggdist)

    filename <- outputfull$filename
    modelname <- outputfull$modelname
    obs_er <- outputfull$obs_er

    n_chains <- outputfull$n_chains
    n_post <- outputfull$n_post

    chain_samples <- 1:n_chains %>% map(~c(rep(.x, n_post))) %>% unlist

    fitfull <- readRDS(here::here("outputs", "fits", filename, paste0("fit_", modelname, ".RDS")))
    data_t <- fitfull$data_t
    N <- data_t$N

    post <- fitfull$post
    par_tab <- fitfull$par_tab
    model_outline <- fitfull$model


    compare <- bind_rows(
        post$mcmc %>% combine %>% as.data.frame  %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>%
             mutate(type = "Posterior distribution") %>% filter(param %in% c("a_d", "b_d", "c_d")),
        purrr::map_df(1:n_post,
            ~model_outline$samplePriorDistributions(par_tab)
        )  %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>%  mutate(type = "Prior distribution")  %>%
        filter(param %in% c("a_d", "b_d", "c_d"))
    )


    ab_function <- function(a, b, c, T) {
        1:T %>% map( 
            function(t) {
                if (t < 14) {
                    titre_init <- 0 + (log(exp(a) +  exp(c)) * (t) / 14);
                } else {
                    titre_init <- 0 + (log(exp(a) * exp(-b/10 * ((t) - 14)) + exp(c)));
                }
            }
        ) %>% unlist
    }


    post_fit <- post$mcmc %>% combine %>% as.data.frame %>% mutate(chain = as.character(chain_samples ))
    a_post <- post_fit[["a_d"]]
    b_post <- post_fit[["b_d"]]
    c_post <- post_fit[["c_d"]]


    T <- 300
    traj_post <- 1:(n_post * n_chains) %>% purrr::map_df(
        ~data.frame(
            time = 1:T,
            value = ab_function(a_post[.x], b_post[.x], c_post[.x], T)
        )
    )
    traj_post_summ <- traj_post %>% group_by(time) %>% mean_qi()

    p1 <- compare %>% filter(type %in% c("Posterior distribution", "Prior distribution")) %>% 
        ggplot() + 
            geom_density(aes(x = value, fill = type), alpha = 0.5) +
            scale_fill_manual(values = c("red", "gray")) +
            facet_wrap(vars(param), scale = "free") + theme_bw() + 
            labs(x = "Value", y = "Density", fill = "Type") + theme_bw() + 
            ggtitle("Antibody kinetics parameters")
    p2 <- traj_post_summ %>%  
        ggplot() + 
            geom_hline(yintercept = 0, size = 2, color = "red") +
        geom_ribbon(aes(ymin = .lower, ymax = .upper, x = time), fill = "red", size = 3, alpha = 0.5) + 
        geom_line(aes(y = value, x = time, color = "red"), size = 2) + 
        theme_bw() + labs(x = expression("Time post-infection, s", y = "Titre boost, f"[ab])) + 
        ggtitle("Antibody trajectories post-infection,") + 
        scale_colour_manual(name = "Line type", 
         values =c('black'='black','red'='red'), labels = c('Simualted trajectory','Posterior trajectory'))

    require(patchwork)

    p1 / p2 + plot_annotation(title = "Simualtion recovery of the antibody kinetics for delta")
    ggsave(here::here("outputs", "fits", filename, "figs", modelname, "ab_kinetics_recov_delta.png"), height = 10, width = 10)

}

plot_abkinetics_vax <- function(outputfull) {
    require(posterior)
    require(bayesplot)
    require(ggdist)

    filename <- outputfull$filename
    modelname <- outputfull$modelname
    obs_er <- outputfull$obs_er

    n_chains <- outputfull$n_chains
    n_post <- outputfull$n_post

    chain_samples <- 1:n_chains %>% map(~c(rep(.x, n_post))) %>% unlist

    fitfull <- readRDS(here::here("outputs", "fits", filename, paste0("fit_", modelname, ".RDS")))
    data_t <- fitfull$data_t
    N <- data_t$N

    post <- fitfull$post
    par_tab <- fitfull$par_tab
    model_outline <- fitfull$model


    compare <- bind_rows(
        post$mcmc %>% combine %>% as.data.frame  %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>%
             mutate(type = "Posterior distribution") %>% filter(param %in% c("a_vax", "b_vax", "c_vax")),
        purrr::map_df(1:n_post,
            ~model_outline$samplePriorDistributions(par_tab)
        )  %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>%  mutate(type = "Prior distribution")  %>%
        filter(param %in% c("a_vax", "b_vax", "c_vax"))
    )


    ab_function <- function(a, b, c, T) {
        1:T %>% map( 
            function(t) {
                if (t < 14) {
                    titre_init <- 0 + (log(exp(a) +  exp(c)) * (t) / 14);
                } else {
                    titre_init <- 0 + (log(exp(a) * exp(-b/10 * ((t) - 14)) + exp(c)));
                }
            }
        ) %>% unlist
    }


    post_fit <- post$mcmc %>% combine %>% as.data.frame %>% mutate(chain = as.character(chain_samples ))
    a_post <- post_fit[["a_vax"]]
    b_post <- post_fit[["b_vax"]]
    c_post <- post_fit[["c_vax"]]


    T <- 300
    traj_post <- 1:(n_post * n_chains) %>% purrr::map_df(
        ~data.frame(
            time = 1:T,
            value = ab_function(a_post[.x], b_post[.x], c_post[.x], T)
        )
    )
    traj_post_summ <- traj_post %>% group_by(time) %>% mean_qi()

    p1 <- compare %>% filter(type %in% c("Posterior distribution", "Prior distribution")) %>% 
        ggplot() + 
            geom_density(aes(x = value, fill = type), alpha = 0.5) +
            scale_fill_manual(values = c("red", "gray")) +
            facet_wrap(vars(param), scale = "free") + theme_bw() + 
            labs(x = "Value", y = "Density", fill = "Type") + theme_bw() + 
            ggtitle("Antibody kinetics parameters")
    p2 <- traj_post_summ %>%  
        ggplot() + 
        geom_ribbon(aes(ymin = .lower, ymax = .upper, x = time), fill = "red", size = 3, alpha = 0.5) + 
        geom_line(aes(y = value, x = time, color = "red"), size = 2) + 
        theme_bw() + labs(x = expression("Time post-infection, s", y = "Titre boost, f"[ab])) + 
        ggtitle("Antibody trajectories post-infection,") + 
        scale_colour_manual(name = "Line type", 
         values =c('black'='black','red'='red'), labels = c('Simualted trajectory','Posterior trajectory'))

    require(patchwork)

    p1 / p2 + plot_annotation(title = "Simualtion recovery of the antibody kinetics for vaccination")
    ggsave(here::here("outputs", "fits", filename, "figs", modelname, "ab_kinetics_recov_vax.png"), height = 10, width = 10)

}


plot_abkinetics_pdelta <- function(outputfull) {
    require(posterior)
    require(bayesplot)
    require(ggdist)

    filename <- outputfull$filename
    modelname <- outputfull$modelname
    obs_er <- outputfull$obs_er

    n_chains <- outputfull$n_chains
    n_post <- outputfull$n_post

    chain_samples <- 1:n_chains %>% map(~c(rep(.x, n_post))) %>% unlist

    fitfull <- readRDS(here::here("outputs", "fits", filename, paste0("fit_", modelname, ".RDS")))
    data_t <- fitfull$data_t
    N <- data_t$N

    post <- fitfull$post
    par_tab <- fitfull$par_tab
    model_outline <- fitfull$model


    compare <- bind_rows(
        post$mcmc %>% combine %>% as.data.frame  %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>%
             mutate(type = "Posterior distribution") %>% filter(param %in% c("a_pd", "b_pd", "c_pd")),
        purrr::map_df(1:n_post,
            ~model_outline$samplePriorDistributions(par_tab)
        )  %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>%  mutate(type = "Prior distribution")  %>%
        filter(param %in% c("a_pd", "b_pd", "c_pd"))
    )


    ab_function <- function(a, b, c, T) {
        1:T %>% map( 
            function(t) {
                if (t < 14) {
                    titre_init <- 0 + (log(exp(a) +  exp(c)) * (t) / 14);
                } else {
                    titre_init <- 0 + (log(exp(a) * exp(-b/10 * ((t) - 14)) + exp(c)));
                }
            }
        ) %>% unlist
    }


    post_fit <- post$mcmc %>% combine %>% as.data.frame %>% mutate(chain = as.character(chain_samples ))
    a_post <- post_fit[["a_pd"]]
    b_post <- post_fit[["b_pd"]]
    c_post <- post_fit[["c_pd"]]


    T <- 300
    traj_post <- 1:(n_post * n_chains) %>% purrr::map_df(
        ~data.frame(
            time = 1:T,
            value = ab_function(a_post[.x], b_post[.x], c_post[.x], T)
        )
    )
    traj_post_summ <- traj_post %>% group_by(time) %>% mean_qi()

    p1 <- compare %>% filter(type %in% c("Posterior distribution", "Prior distribution")) %>% 
        ggplot() + 
            geom_density(aes(x = value, fill = type), alpha = 0.5) +
            scale_fill_manual(values = c("red", "gray")) +
            facet_wrap(vars(param), scale = "free") + theme_bw() + 
            labs(x = "Value", y = "Density", fill = "Type") + theme_bw() + 
            ggtitle("Antibody kinetics parameters")
    p2 <- traj_post_summ %>%  
        ggplot() + 
        geom_ribbon(aes(ymin = .lower, ymax = .upper, x = time), fill = "red", size = 3, alpha = 0.5) + 
        geom_line(aes(y = value, x = time, color = "red"), size = 2) + 
        theme_bw() + labs(x = expression("Time post-infection, s", y = "Titre boost, f"[ab])) + 
        ggtitle("Antibody trajectories post-infection,") + 
        scale_colour_manual(name = "Line type", 
         values =c('black'='black','red'='red'), labels = c('Simualted trajectory','Posterior trajectory'))

    require(patchwork)

    p1 / p2 + plot_annotation(title = "Simualtion recovery of the antibody kinetics for pre-delta infection")
    ggsave(here::here("outputs", "fits", filename, "figs", modelname, "ab_kinetics_recov_pdelta.png"), height = 10, width = 10)

}



postprocessFigs <- function(filename, modelname, n_chains) {
    postprocess_run( filename, modelname, n_chains)
    outputfull <- readRDS(file = here::here("outputs", "fits", filename, paste0("pp_", modelname, ".RDS")))

    plot_trace(outputfull)
    plot_abkinetics(outputfull)
    plot_inf_rec(outputfull)
    plot_exp_rec(outputfull)
    plot_exp_times_rec(outputfull)
    plot_cop_rec(outputfull)
}


postprocessFigs_wave2 <- function(filename, modelname, n_chains) {
    postprocess_run( filename, modelname, n_chains)
    outputfull <- readRDS(file = here::here("outputs", "fits", filename, paste0("pp_", modelname, ".RDS")))

    plot_trace(outputfull)
    plot_abkinetics_delta(outputfull)
    plot_abkinetics_vax(outputfull)
    plot_abkinetics_pdelta(outputfull)
    plot_inf_rec(outputfull)
    plot_exp_rec(outputfull)
    plot_exp_times_rec(outputfull)
    plot_cop_rec(outputfull)
}