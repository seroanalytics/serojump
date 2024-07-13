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
postprocess_runInf <- function(fitfull, filename, modelname, n_chains, priorPred = FALSE) {

    library(ggdist)

    #filename <- "hpc/nih_2024_inf/p3"
    #modelname <- "h3"
    #n_chain <- 4
    #fitfull <- fitfull_h3 

   # fitfull$post$jump[[1]][,781 ]

    post <- fitfull$post
    data_t <- fitfull$data_t
    model <- fitfull$model

    knowndate <- data_t$knownInfsDate
    initialTitreValue <- data_t$initialTitreValue

    knownid <- which(knowndate != -1)
    N <- data_t$N
    N_data <- data_t$N_data

    n_post <- post$mcmc[[1]] %>% nrow
    post_exp_combine <- post$jump 
    post_inf_combine <- post$inf 

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
    outputfull$fit_states %>% filter(id == 111)


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
        modelname = modelname,
        filename = filename,
        fit_states = fit_states,
        fit_obs = post_obstitre,
        n_chains = n_chains,
        n_post = n_post)

    if (!priorPred) {
        saveRDS(postprocess, file = here::here("outputs", "fits", filename, paste0("pp_", modelname, ".RDS")))
    } else {
        saveRDS(postprocess, file = here::here("outputs", "fits", filename, paste0("pp_prior_", modelname, ".RDS")))
    }

}

plot_abkineticsInf <- function(outputfull, fitfull, fig_folder) {
    require(posterior)
    require(bayesplot)
    require(ggdist)

    filename <- outputfull$filename
    modelname <- outputfull$modelname
    obs_er <- outputfull$obs_er

    n_chains <- outputfull$n_chains
    n_post <- outputfull$n_post

    chain_samples <- 1:n_chains %>% map(~c(rep(.x, n_post))) %>% unlist

    data_t <- fitfull$data_t
    N <- data_t$N

    post <- fitfull$post
    par_tab <- fitfull$par_tab
    model_outline <- fitfull$model


    compare <- bind_rows(
        post$mcmc %>% lapply(as.data.frame) %>% do.call(rbind, .) %>% as.data.frame  %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>%
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


    post_fit <- post$mcmc %>% lapply(as.data.frame) %>% do.call(rbind, .) %>% as.data.frame %>% mutate(chain = as.character(chain_samples ))
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
    ggsave(here::here("outputs", "fits", filename,  "figs", modelname, fig_folder, "ab_kinetics_recov.png"), height = 10, width = 10)

}

plot_exp_recInf <- function(outputfull, fitfull, fig_folder) {


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

plot_exp_times_recInf <- function(outputfull, fitfull, fig_folder) {

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
 
    figB / figC + plot_annotation(tag_levels = "A") 
    ggsave(here::here("outputs", "fits", filename, "figs", modelname, fig_folder,  "exposure_time_recov.png"), height = 10, width = 10)
}

plot_inf_recInf <- function(outputfull, fitfull, fig_folder, scale_ab = NULL) {
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


    #figD <- dataplt_inf %>% select(id, exp_time_sim, exp_time) %>% 
    #    rename(Data = exp_time_sim, `Model recovery` = exp_time) %>%
    #    pivot_longer(!id, names_to = "type", values_to = "time") %>%
   # ggplot() + 
    #        geom_density(aes(x = time, fill = type), size = 2, alpha = 0.5, shape = 21) + 
     #       theme_bw()  + 
      #      labs(x = "Time in study", y = "Density", fill = "") + 
       #     ggtitle("Recovery of infection timings")

    figC + plot_annotation(tag_levels = "A")
    ggsave(here::here("outputs", "fits", filename,  "figs", modelname, fig_folder, "infection_recov.png"), height = 10, width = 10)
}


plot_cop_recInf <- function(outputfull, fitfull, fig_folder, scale_ab = NULL) {

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
  
    cop_exp_sum_plot_all <- map_df(1:length(model_outline$observationalModel), 
        function(i) {
        biomarker <- model_outline$observationalModel[[i]]$biomarker

        titre_cop_sum <- fit_states %>% filter(!!sym(biomarker) != -1) %>%
            group_by(id) %>% mean_qi(!!sym(biomarker) ) %>% arrange(!!biomarker) %>%
            rename(titre_val = !!(biomarker))
        u_ids <- titre_cop_sum$id

        df_data <- data.frame(
            id = 1:data_t$N,
            start_titre = data_t$initialTitreValue %>% as.data.frame %>% pull(!!biomarker),
            known_inf = data_t$knownInfsVec
        ) %>% mutate(known_inf = ifelse(known_inf == 0, "Not known", "Known"))

        df_data_post <- fit_states %>% filter(inf_ind == 1) %>% select(id, inf_ind) %>%
                 summarise(inf_post = mean(inf_ind), n = n(), .by = "id") %>%
                left_join(df_data) %>% mutate(prop = n / n_length)

        cop_exp_sum_plot <- titre_cop_sum %>% left_join(df_data_post) %>% mutate(id = factor(id, levels = u_ids)) %>%
            filter(!is.na(n)) %>% mutate(biomarker = !!biomarker)

        })

    figA <- cop_exp_sum_plot_all %>% 
        ggplot() + 
        # geom_smooth(aes(x = titre_val, y = inf_post)) +
        geom_linerange(data = cop_exp_sum_plot_all, 
            aes(y = prop, xmin = .lower, xmax = .upper, alpha = prop), size = 1) + 
        geom_point(data = cop_exp_sum_plot_all, aes(x = titre_val, y = prop, alpha = prop)) + 
        theme_bw() + 
        ylim(0, 1) +
        facet_wrap(vars(biomarker)) + 
        labs(x = expression("Titre at exposure"), y = expression("Posterior probability of infection for correlate of protection, f"[cop]*"(Y"[j]^0*", "*theta[cop]*")"), color = "Curve type")

    if (!is.null(scale_ab)) {
        figA <- figA + scale_x_continuous(breaks = scale_ab %>% as.numeric, labels = scale_ab %>% names)
    }


   # nih_inf_raw <- read_excel(path = here::here("data", "nih_2024_XX", "2022_2023_Flu_Swabs.xlsx") )
   # pid_cols <- data_titre %>% select(pid, id)
   # swab_info <- nih_inf_raw %>% filter(year == 2023) %>% left_join(pid_cols) %>% select(pid, id, swab_virus) %>% unique
   # swab_info %>% left_join(cop_exp_sum_plot_all %>% mutate(id = as.numeric(id)), by = "id" ) %>% 
   #       ggplot() + 
   #         geom_linerange( 
   #         aes(y = prop, xmin = .lower, xmax = .upper, color = swab_virus), size = 1, alpha = 0.7) + 
   #         geom_point(aes(x = titre_val, y = prop, color = swab_virus), alpha = 0.7, size = 3) + 
   #         theme_bw() 

    ggsave(here::here("outputs", "fits", filename,  "figs", modelname, fig_folder, "cop_recov.png"), height = 10, width = 10)

}


plot_titre_expInf <- function(outputfull, fitfull, fig_folder, scale_ab = NULL) {
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
    if (!is.null(scale_ab)) {
        p1 <- p1 + scale_x_continuous(breaks = scale_ab %>% as.numeric, labels = scale_ab %>% names)
    }

    p1 
    ggsave(here::here("outputs", "fits", filename,  "figs", modelname, fig_folder, "titre_exp_recovery.png"), height = 10, width = 10)

}   

plot_titre_obsInf <- function(outputfull, fitfull, fig_folder, scale_ab = NULL) {
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
    if (!is.null(scale_ab)) {
        p1 <- p1 + scale_y_continuous(breaks = scale_ab %>% as.numeric, labels = scale_ab %>% names)
    }
    ggsave(here::here("outputs", "fits", filename,  "figs", modelname, fig_folder, "titre_obs.png"), height = 10, width = 10)

}


postprocessFigsInf <- function(filename, modelname, n_chains, scale_ab = NULL) {
   # filename <- "test/nih_2024"
   # modelname <- "h3_profile"
   # n_chains <- 4
   #filename <- "hpc/nih_2024"
  # modelname <- "h3"
   #filename <-  "local/nih_2024_inf/test"
  # modelname <- "h3"
  # n_chains <- 4
    filename <- "hpc/nih_2024_inf/p3"
    modelname <- "h3"
    n_chains <- 4

    dir.create(here::here("outputs", "fits", filename,  "figs", modelname), recursive = TRUE, showWarnings = FALSE)
    fitfull_pp <- readRDS(here::here("outputs", "fits", filename, paste0("fit_prior_", modelname, ".RDS")))
    fitfull <- readRDS(here::here("outputs", "fits", filename, paste0("fit_", modelname, ".RDS")))

    postprocess_runInf(fitfull_pp, filename, modelname, n_chains, TRUE)
    postprocess_runInf(fitfull, filename, modelname, n_chains)

    outputfull_pp <- readRDS(file = here::here("outputs", "fits", filename, paste0("pp_prior_", modelname, ".RDS")))
    outputfull <- readRDS(file = here::here("outputs", "fits", filename, paste0("pp_", modelname, ".RDS")))

 #   scale_ab <- NULL
    postconvergeFigs_Inf(fitfull_pp, filename, modelname, TRUE)
    postconvergeFigs_Inf(fitfull, filename, modelname)

    fig_folder <- "pp"
    plot_titre_obsInf(outputfull_pp, fitfull_pp, fig_folder, scale_ab)
    plot_titre_expInf(outputfull_pp, fitfull_pp, fig_folder, scale_ab)
    plot_abkinetics_trajectoriesInf(outputfull_pp, fitfull_pp, fig_folder)
    plot_inf_recInf(outputfull_pp, fitfull_pp, fig_folder,  scale_ab)
    plot_exp_times_recInf(outputfull_pp, fitfull_pp, fig_folder)
    plot_cop_recInf(outputfull_pp, fitfull_pp, fig_folder, scale_ab)




    fig_folder <- "post"
    plot_titre_obsInf(outputfull, fitfull, fig_folder, scale_ab)
    plot_titre_expInf(outputfull, fitfull, fig_folder, scale_ab)
    plot_abkinetics_trajectoriesInf(outputfull, fitfull, fig_folder)
    plot_inf_recInf(outputfull,fitfull, fig_folder,  scale_ab)
    plot_exp_times_recInf(outputfull, fitfull, fig_folder)
    plot_cop_recInf(outputfull, fitfull, fig_folder, scale_ab)

}

postconvergeFigs_Inf <- function(fitfull, filename, modelname, priorPred = FALSE) {

    if (!priorPred) {
        filepath <- here::here("outputs", "fits", filename, "figs",  modelname, "post")
        dir.create(filepath, showWarnings = FALSE, recursive = TRUE)
    } else {
        filepath <- here::here("outputs", "fits", filename, "figs",  modelname, "pp")
        dir.create(filepath, showWarnings = FALSE, recursive = TRUE)
    }
    seroJumpPostDaig(fitfull, filepath)

}


plot_abkinetics_trajectoriesInf <- function(outputfull, fitfull, fig_folder) {

    require(posterior)
    require(bayesplot)
    require(ggdist)

    filename <- outputfull$filename
    modelname <- outputfull$modelname

    n_chains <- outputfull$n_chains
    n_post <- outputfull$n_post

    chain_samples <- 1:n_chains %>% map(~c(rep(.x, n_post))) %>% unlist

    data_t <- fitfull$data_t
    N <- data_t$N

    post <- fitfull$post
    par_tab <- fitfull$par_tab
    model_outline <- fitfull$model
    post_fit <- post$mcmc %>% lapply(as.data.frame) %>% do.call(rbind, .) %>% as.data.frame %>% mutate(chain = as.character(chain_samples ))

    post$jump
    ###
    posteriorsAllExposure <- map_df(1:length(model_outline$abkineticsModel),
        function(name1) {
            pars_extract <- model_outline$abkineticsModel[[name1]]$pars
            functionalForm <- model_outline$abkineticsModel[[name1]]$funcForm
            biomarker <- model_outline$abkineticsModel[[name1]]$biomarker
            exposureType <- model_outline$abkineticsModel[[name1]]$exposureType
            cat(pars_extract)
            compare <- bind_rows(
                post$mcmc %>% lapply(as.data.frame) %>% do.call(rbind, .) %>% as.data.frame  %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>%
                    mutate(type = "Posterior distribution") %>% filter(param %in% pars_extract)
            #  purrr::map_df(1:n_post,
            #      ~model_outline$samplePriorDistributions(par_tab)
            #  )  %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>%  mutate(type = "Prior distribution")  %>%
            #  filter(param %in% pars_extract)
            )

            ab_function <- function(T, pars) {
                1:T %>% map( 
                    ~functionalForm(0, .x, pars)
                ) %>% unlist
            }

            post_par <- data.frame(post_fit[[pars_extract[1]]])
            if (length(pars_extract) > 1) {
                for (i in 2:length(pars_extract)) {
                    post_par <- cbind(post_par, post_fit[[pars_extract[i]]])
                }   
            }

            T <- 300
            traj_post <- 1:(100) %>% purrr::map_df(
                ~data.frame(
                    time = 1:T,
                    value = ab_function(T, as.numeric(post_par[.x, ]))
                )
            )
            traj_post_summ <- traj_post %>% group_by(time) %>% mean_qi() %>% mutate(exposure_type = exposureType) %>% 
                mutate(biomarker = biomarker)
        }
    )

    p1 <- posteriorsAllExposure %>%  
        ggplot() + 
        geom_ribbon(aes(ymin = .lower, ymax = .upper, x = time, fill = exposure_type), size = 3, alpha = 0.3) + 
        geom_line(aes(y = value, x = time, color = exposure_type), size = 2) + 
        theme_bw() + labs(x = expression("Time post-infection/since bleed, s", y = "Titre change, f"[ab])) + 
        ggtitle("Antibody kinetic trajectories") +
        facet_wrap(vars(biomarker))
    ggsave(here::here("outputs", "fits", filename, "figs",  modelname, fig_folder, "ab_kinetics_trajectories.png"), height = 10, width = 10)

}



plot_abkinetics_trajectories2Inf <- function(outputfull, fitfull, fig_folder) {

    require(posterior)
    require(bayesplot) 
    require(ggdist)

    filename <- outputfull$filename
    modelname <- outputfull$modelname

    n_chains <- outputfull$n_chains
    n_post <- outputfull$n_post

    chain_samples <- 1:n_chains %>% map(~c(rep(.x, n_post))) %>% unlist
    M <- length(chain_samples)
    data_t <- fitfull$data_t
    N <- data_t$N


    post <- fitfull$post
    par_tab <- fitfull$par_tab
    model_outline <- fitfull$model
    post_fit <- post$mcmc %>% lapply(as.data.frame) %>% do.call(rbind, .) %>% as.data.frame %>% mutate(chain = as.character(chain_samples ))

    ## FOR individual i
    require(data.table)
    bio_all <- model_outline$infoModel$biomarkers

    # Get serological data to plot!
    data_fit <- map_df(1:N, 
        function(i) {
            map_df(1:length(bio_all), 
                function(b) {
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
                }
            )
        }
    )

    # Extrat the order of the exposure for each individual


    exposures <- model_outline$infoModel$exposureInfo %>% map(~.x$exposureType) %>% unlist
    fit_states_dt <- as.data.table(outputfull$fit_states)
    df_know_exp <- map_df(1:length(exposures),
        function(e) { 
            exp <- exposures[e]
            known_inf_e  <- model_outline$infoModel$exposureInfo[[e]]$known_inf
            data.frame(
                id = 1:N,
                time = known_inf_e,
                type = exp
            )
        }
    )
    known_types <- df_know_exp$type


    S <- 100
    sample_s <- sample(1:M, S)
    fit_states_dt_trim <- fit_states_dt[sample %in% sample_s]

    exposures_fit <- model_outline$infoModel$exposureFitted
    exposures <-  model_outline$exposureTypes

    require(future)

    plan(multisession, workers = 8)

 
    cat("Get order of all entries")

    df_exposure_order <- future_map(1:N, 
        function(i) {
            df_know_exp_i <- df_know_exp %>% filter(id == i) %>% filter(time > -1)
            fit_i <- fit_states_dt[id == i]

            map (1:S, 
            function(s) {
                # Get the inferred inection time first
                fit_i_s <- fit_i[sample == sample_s[s]]
                if (fit_i_s$inf_ind == 1){
                    fit_i_s_long <- melt(fit_i_s, measure.vars = bio_all, variable.name = "biomarker", value.name = "titre")
                    fit_i_s_long[, exp_type := exposures_fit]
                    fit_i_s_long <- fit_i_s_long[, .(id, sample, time = inf_time, exp_type, biomarker, titre)]
                } else {
                    fit_i_s_long <- data.table()
                }

               if (nrow(df_know_exp_i) > 0) {
                    new_rows_list <- lapply(1:nrow(df_know_exp_i), function(j) {
                        data.table(
                            id = i,
                            sample = sample_s[s],
                            time = df_know_exp_i$time[j],
                            exp_type = df_know_exp_i$type[j],
                            biomarker = bio_all,
                            titre = NA
                        )
                    })
                    fit_i_s_long <- rbindlist(c(list(fit_i_s_long), new_rows_list), fill = TRUE)
                }

                setorder(fit_i_s_long, time)
                fit_i_s_long
            }
        ) %>% rbindlist
            
        }
    ) %>% rbindlist
      # Perform final mutations and groupings
    df_exposure_order[, time := round(time, 0)]
    df_exposure_order[, row_id := seq_len(.N), by = .(id, biomarker, sample)]
    df_exposure_order <- df_exposure_order[order(id, biomarker, time)]


    df_exposure_order[id == i & sample == s & biomarker == bio_i]   


    bio_map_ab <- model_outline$abkineticsModel %>% map(~.x$biomarker) %>% unlist
    exp_map_ab <- model_outline$abkineticsModel %>% map(~.x$exposureType) %>% unlist
    df_map_ab <- data.frame(
        k = 1:length(bio_map_ab),
        bio = bio_map_ab,
        exp = exp_map_ab
    )

    par_mean <- bind_rows(
        post$mcmc %>% lapply(as.data.frame) %>% do.call(rbind, .) %>% as.data.frame  %>% summarise(across(everything(), mean))  %>% 
            mutate(type = "Posterior distribution") 
    )

    par_sample <- bind_rows(
        post$mcmc %>% lapply(as.data.frame) %>% do.call(rbind, .) %>% as.data.frame  
    )

    # Get the simulated trajectories for each individual and exposure type and biomarker
    

    i <- 1
    plan(multisession, workers = 12)


    setDT(data_fit)
    setDT(df_map_ab)

    # Pre-compute data_fit and df_map_ab selections
    data_fit_list <- split(data_fit, data_fit$id)
    df_map_ab_list <- split(df_map_ab, df_map_ab$bio)

    # Get known exposures
    exposures_know <- setdiff(exposures, exposures_fit)
    df_ids_plot_known <- map_df(1:length(exposures_know),
        function(i) {
            df_exposure_order_k <- df_exposure_order %>% filter(exp_type == exposures_know[i]) 
            ids_take <- df_exposure_order_k %>% pull(id) %>% unique
            if (length(ids_take) > 20) {
                ids_take <- sample(ids_take, 20)
            }
            df_ids_plot <- data.frame(
                id = ids_take,
                type = exposures_know[i]
            )
        }
    )

    id_skip <- df_ids_plot_known$id

    # Get inferred trajectories for missed infections for each individual 
    df_exposure_order_intense <- df_exposure_order %>% filter(!id %in% id_skip) %>% select(!biomarker) %>% filter(exp_type %in% exposures_fit) %>%
        summarise(prob = n() / (S * length(bio_all)), .by = id) %>% mutate(
            type = cut(prob, c(0, 0.25, 0.75, 1), labels = c("Low", "Medium", "High"))
        )
    type_infer <- c("Low", "Medium", "High")
    df_ids_plot_infer <- map_df(1:3,
        function(i) {
            df_exposure_order_k <- df_exposure_order_intense %>% filter(type == type_infer[i]) 
            ids_take <- df_exposure_order_k %>% pull(id) %>% unique
            if (length(ids_take) > 20) {
                ids_take <- sample(ids_take, 20)
            }
            df_ids_plot <- data.frame(
                id = ids_take,
                type = type_infer[i]
            )
        }
    )

    df_ids_plot <- bind_rows(df_ids_plot_known, df_ids_plot_infer)

    ########

    ####
    df_ids_plot_i <- df_ids_plot %>% filter(type == "High")

    id_i <- 217

    df_exposure_order[id == id_i]
    data_fit %>% filter(id == id_i)
    data_t$initialTitreValue[id_i, ]
    data_t$initialTitreTime[id_i]
    data_t$titre_list[[id_i]]
    data_t$times_list[[id_i]]
    data_t$id_full[[id_i]]
    data_t$T

    df_traj_post_ind <- map(
        df_ids_plot_i$id,
        function(i) {
                map(bio_all, 
                    function(bio_i) {
                        map(sample_s,
                            function(s) {
                                df_exposure_order_i <- df_exposure_order[id == i & sample == s & biomarker == bio_i] %>% arrange(time, .by_group = TRUE)
                        
                                times <- c(df_exposure_order_i[["time"]], 365)
                                timesince_vec <- times %>% diff

                                titre_traj <- NULL
                                titre_anchor <- NULL

                                for (j in seq_len(nrow(df_exposure_order_i))) {
                                    exp_type_i <- df_exposure_order_i[j, exp_type]
                                    id_key <- df_map_ab_list[[bio_i]][exp == exp_type_i]$k
                                    ab_func <- model_outline$abkineticsModel[[id_key]]$funcForm
                                    pars_extract <- model_outline$abkineticsModel[[id_key]]$pars
                                    par_in <- as.numeric(par_sample[s, pars_extract])
                                 #   cat(par_in)
                                    titre_start <- data_fit_list[[i]][row_num == 1, titre]

                                    if (df_exposure_order_i[j, row_id] == 1) {

                                        time_anchor <- 1
                                        titre_start <- data_fit_list[[i]][row_num == 1, titre]
                                        titre_traj <- rep(NA, times[1])
                                        traj_type <- rep(NA, times[1])
                                        vector_titre <- unlist(lapply(seq_len(timesince_vec[j]), function(t) ab_func(titre_start, t, par_in)))

                                        titre_traj <- c(titre_traj, vector_titre)
                                        traj_type <- c(traj_type, rep(exp_type_i, length(vector_titre)))
                                        titre_anchor <- ifelse(length(vector_titre) == 0, titre_start, vector_titre[length(vector_titre)])
                                        time_anchor <- length(titre_traj)
                                    } else {

                                        vector_titre <- unlist(lapply(seq_len(timesince_vec[j]), function(t) ab_func(titre_anchor, t, par_in)))
                                        titre_traj <- c(titre_traj, vector_titre)
                                        traj_type <- c(traj_type, rep(exp_type_i, length(vector_titre)))
                                        titre_anchor <- ifelse(length(vector_titre) == 0, titre_start, vector_titre[length(vector_titre)])
                                    }
                                }

                                data.table(
                                    id = i,
                                    sample = s,
                                    t = 1:365,
                                    biomarker = bio_i,
                                    type = traj_type,
                                    titre_traj = titre_traj
                                ) %>% filter(!is.na(titre_traj))
                            }
                        ) %>% rbindlist
                    }
                )  %>% rbindlist
        }
    ) %>% rbindlist



    p1 <- df_traj_post_ind  %>% filter( biomarker == bio_all[1]) %>%
        ggplot() +
            geom_line(aes(x = t, y = titre_traj, color = type, group = sample), alpha = 0.05) + facet_wrap(vars(id)) + 
            geom_point(aes(x = times, y = titre), data = data_fit %>% filter(id %in% df_ids_plot_i$id), shape = 21, size = 2, fill = "gray") + theme_bw()
            fig_folder <- "post"
    ggsave(here::here("outputs", "fits", filename, "figs",  modelname, fig_folder, "ab_kinetics_trajectories_2.png"), height = 10, width = 10)



}
