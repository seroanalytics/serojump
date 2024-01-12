extract_long_post <- function(mcmc_out, ...) {
    mcmc_out %>% combine %>% as.data.frame %>%
        dplyr::select(...) %>% 
        pivot_longer(everything(), names_to = "var", values_to = "value")
}
extract_post <- function(mcmc_out, ...) {
    mcmc_out %>% combine %>% as.data.frame %>%
        dplyr::select(...) 
}

add_par_df <- function(...) {
    # Convert the ellipsis arguments to a list
    args_list <- list(...)
    names(args_list) <- c("par_name", "lb", "ub", "dist", "dist_par1", "dist_par2")
    df <- as.data.frame(args_list)
    df$dist_par1 <- as.character(df$dist_par1)
    df$dist_par2 <- as.character(df$dist_par2)
    df$part_type <- "prior"
  return(df)
}

cal_lprior_non_centered <- function(par_tab, params) {
    p <- 0
    P <- nrow(par_tab)
    names(params) <- par_tab$par_name
    for (i in 1:P) {
        dist_name <- paste0("d", par_tab[i, 4])
        my_dist_name <- get(dist_name)
        if (par_tab[i, 4] == "exp") {
            p <- p + do.call(my_dist_name, list(as.numeric(params[i]), as.numeric(par_tab[i, 5]), log = TRUE) )
        } else {
            p <- p + do.call(my_dist_name, list(as.numeric(params[i]), as.numeric(par_tab[i, 5]), as.numeric(par_tab[i, 6]), log = TRUE) )
        }
    }
    p
}


get_sample_non_centered <- function(par_tab) {
    P <- nrow(par_tab)
    s <- vector(mode = "numeric", length = P)
    names(s) <- par_tab$par_name
    for (i in 1:P) {
        dist_name <- paste0("r", par_tab[i, 4])
        my_dist_name <- get(dist_name)
        if (par_tab[i, 4] == "exp") {
            s[i] <- do.call(my_dist_name, list(1, as.numeric(par_tab[i, 5])) )
            while (check_boundaries(s[i], par_tab[i, 2],  par_tab[i, 3])) {
                s[i] <- do.call(my_dist_name, list(1, as.numeric(par_tab[i, 5]) ))
            }
        } else {
            s[i] <- do.call(my_dist_name, list(1, as.numeric(par_tab[i, 5]),  as.numeric(par_tab[i, 6])) )
            while (check_boundaries(s[i], par_tab[i, 2],  par_tab[i, 3])) {
                s[i] <- do.call(my_dist_name, list(1,  as.numeric(par_tab[i, 5]),  as.numeric(par_tab[i, 6])) )
            }
        }
    }
    s
}

check_boundaries <- function(x, lb, ub) {
    (x < lb) | (x > ub);
}


clean_simulated_rjmcmc <- function(modelname_sim, obs_er, prob_known, known_exp = FALSE) {

    modeli <- readRDS(here::here("outputs", "sim_data", modelname_sim, "inputs.RDS"))
    res <- readRDS(file = here::here("outputs", "sim_data", modelname_sim, paste0("sim_data_", obs_er, ".rds")))

    name <- modeli$name
    protection_curve <- modeli$protection_curve
    pre_season_titre <- modeli$simpar$start_titre
    N <- modeli$simpar$N    
    T <- modeli$simpar$T    
    #endregionplot_subset_individuals_history(res$biomarker_states, res$immune_histories_long, subset=30, demography)
    #ggsave(here::here("outputs", "sim", modelname, "plot_subset_individuals_history.pdf"))

    # Determine known infections
    no_inf <- sum(res$immune_histories %>% apply(1, sum) > 0)
    no_known <- round(no_inf * prob_known)
    known_post <- rdunif(no_known, 1, no_inf)
    known_times <- which(t(res$immune_histories[which(res$immune_histories %>% apply(1, sum) == 1), ][known_post, ]) == 1) %% T
    known <- rep(0, N)
    dayinf <- rep(-1, N)
    known[which(res$immune_histories %>% apply(1, sum) == 1)][known_post] <- 1
    dayinf[which(res$immune_histories %>% apply(1, sum) == 1)][known_post] <- known_times

    # Add in titre daya and time of bleed data
    initialTitreValue <- res$observed_biomarker_states %>% group_by(i) %>% filter(t == min(t)) %>% .[["value"]]
    initialTitreTime <- res$observed_biomarker_states %>% group_by(i) %>% filter(t == min(t)) %>% .[["t"]]
    endTitreTime <- res$observed_biomarker_states %>% group_by(i) %>% filter(t == max(t)) %>% .[["t"]]

    if (length(initialTitreValue) > N) {
        stop("ERROR: initialTitreValue wrong size ")
    }
    if (length(initialTitreTime) > N) {
        stop("ERROR: initialTitreTime wrong size ")
    }
    if (length(endTitreTime) > N) {
        stop("ERROR: endTitreTime wrong size ")
    }

    df_immune_hist <- res$immune_histories %>% as.data.frame %>% mutate(i = 1:N) %>%
        pivot_longer(!i, names_to = "t", values_to = "inf") %>%
        mutate(t = as.numeric(substr(t, 2, 4)))

    res$observed_biomarker_states %>% 
        ggplot() + geom_point(aes(x = t, y = i, color = value)) + 
        geom_point(data = df_immune_hist %>% filter(inf == 1), aes(t, i), color = "red") + theme_bw()


    observed_biomarker_statesStudy <- res$observed_biomarker_states %>% filter(t != 1)

    titre_obs <- observed_biomarker_statesStudy$observed
    titre_true <- observed_biomarker_statesStudy$value
    times_full <- observed_biomarker_statesStudy$t
    id_full <- observed_biomarker_statesStudy$i
    N_data <- length(id_full)

    if (known_exp) {
        knownExpVec <- modeli$simpar$exp %>% as.data.frame %>% 
        mutate(i = 1:N) %>% pivot_longer(!i, names_to = "t", values_to = "exp") %>% mutate(t = as.numeric(substr(t, 2, 4))) %>% 
        filter(exp == 1) %>% complete(i = 1:N, fill = list(t = -1, exp = -1)) %>% pull(t)
    } else {
        knownExpVec <- NA
    }

    # Things needed, data for likelihood and sampler
    data_t <- list(
        N = N,
        knownExpVec = knownExpVec,
        knownInfsVec = known,        
        knownInfsN = sum(known),
        knownInfsDate = dayinf,

        N_data = N_data,
        initialTitreValue = initialTitreValue,
        initialTitreTime = initialTitreTime,
        endTitreTime = endTitreTime,

        titre_full = titre_true,
        times_full = times_full,
        id_full = id_full
    )
}


createModelRJCMC <- function(ab_ll, par_tab, exp_prior) {
    model_type <- list()
    
    model_type$evaluateLogLikelihood <- ab_ll 
    model_type$lowerParSupport_fitted <- par_tab$lb
    model_type$upperParSupport_fitted <- par_tab$ub

    model_type$namesOfParameters <- par_tab$par_name

    model_type$samplePriorDistributions = function(datalist) {
        get_sample_non_centered(par_tab)
    }

    model_type$evaluateLogPrior <- function(params, jump, datalist) {
        cal_lprior_non_centered(par_tab, params)
    }

    model_type$initialiseJump <- function(datalist) {
        c(datalist$knownInfsDate)
    }

    model_type$exposureFunctionSample <- function() {
        dist_name <- paste0("r", exp_prior[1, 3])
        my_dist_name <- get(dist_name)
        s <- do.call(my_dist_name, list(1, as.numeric(exp_prior[1, 4]),  as.numeric(exp_prior[1, 5])) )
        s
    }

    model_type$exposureFunctionDensity <- function(jump_i) {
        dist_name <- paste0("d", exp_prior[1, 3])
        my_dist_name <- get(dist_name)
        d <- do.call(my_dist_name, list(as.numeric(jump_i), as.numeric(exp_prior[1, 4]), as.numeric(exp_prior[1, 5])) )
        d
    }

    model_type
}


createModelRJCMCFull <- function(ab_ll, par_tab, exp_prior, cop_func) {
    model_type <- list()
    
    model_type$evaluateLogLikelihood <- ab_ll 
    model_type$lowerParSupport_fitted <- par_tab$lb
    model_type$upperParSupport_fitted <- par_tab$ub

    model_type$namesOfParameters <- par_tab$par_name

    model_type$samplePriorDistributions = function(datalist) {
        get_sample_non_centered(par_tab)
    }

    model_type$evaluateLogPrior <- function(params, jump, datalist) {
        cal_lprior_non_centered(par_tab, params)
    }

    model_type$copFunction <- cop_func

    model_type$initialiseJump <- function(datalist) {
        c(datalist$knownInfsDate)
    }

    model_type$exposureFunctionSample <- function() {
        dist_name <- paste0("r", exp_prior[1, 3])
        my_dist_name <- get(dist_name)
        s <- do.call(my_dist_name, list(1, as.numeric(exp_prior[1, 4]),  as.numeric(exp_prior[1, 5])) )
        s
    }

    model_type$exposureFunctionDensity <- function(jump_i) {
        dist_name <- paste0("d", exp_prior[1, 3])
        my_dist_name <- get(dist_name)
        d <- do.call(my_dist_name, list(as.numeric(jump_i), as.numeric(exp_prior[1, 4]), as.numeric(exp_prior[1, 5])) )
        d
    }

    model_type
}



safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
scales::show_col(safe_colorblind_palette)