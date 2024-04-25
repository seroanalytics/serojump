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

add_par_pool_df <- function(...) {
    # Convert the ellipsis arguments to a list
    args_list <- list(...)
    names(args_list) <- c("par_name", "length", "mean", "dist_name", "dist_sd", "dist_sd_par1", "dist_sd_par2")
    
    list_par_pool <- list(mode = "list", length = args_list$length)
    for (i in 1:args_list$length) {
        args_list_temp <- list(par_name = paste0(args_list$par_name, "_", i), lb = -5, ub = 5, dist = "norm",
            dist_par1 = as.character(args_list$mean), dist_par2 = args_list$dist_name)
        list_par_pool[[i]] <- as.data.frame(args_list_temp)
        list_par_pool[[i]]$part_type <- "hprior"
    }
    list_par_reg <- list(args_list$dist_name, 0, 5, args_list$dist_sd, args_list$dist_sd_par1, args_list$dist_sd_par2)
    names(list_par_reg) <-  c("par_name", "lb", "ub", "dist", "dist_par1", "dist_par2")
    list_par_reg$dist_par1 <- as.character(list_par_reg$dist_par1)
    list_par_reg$dist_par2 <- as.character(list_par_reg$dist_par2)
    list_par_reg$part_type <- "prior"

    df <- bind_rows(as.data.frame(list_par_reg)) %>% bind_rows(bind_rows(list_par_pool) )
    return(df)
}

# Can add heirarchical priors with the functions below 
# add_par_pool_df_non_centered("boost_naive_pvnt", length = 5, 0, "boost_naive_sigma_pvnt", "exp", 5, NA)

add_par_pool_df_non_centered <- function(...) {
    # Convert the ellipsis arguments to a list
    args_list <- list(...)
    names(args_list) <- c("par_name", "length", "mean", "dist_name", "dist_sd", "dist_sd_par1", "dist_sd_par2")
    
    list_par_pool <- list(mode = "list", length = args_list$length)
    for (i in 1:args_list$length) {
        args_list_temp <- list(par_name = paste0(args_list$par_name, "_", i), lb = -5, ub = 5, dist = "norm",
            dist_par1 = as.character(args_list$mean), dist_par2 = "1")
        list_par_pool[[i]] <- as.data.frame(args_list_temp)
        list_par_pool[[i]]$part_type <- "hprior"
    }
    list_par_reg <- list(args_list$dist_name, 0, 5, args_list$dist_sd, args_list$dist_sd_par1, args_list$dist_sd_par2)
    names(list_par_reg) <-  c("par_name", "lb", "ub", "dist", "dist_par1", "dist_par2")
    list_par_reg$dist_par1 <- as.character(list_par_reg$dist_par1)
    list_par_reg$dist_par2 <- as.character(list_par_reg$dist_par2)
    list_par_reg$part_type <- "prior"

    df <- bind_rows(as.data.frame(list_par_reg)) %>% bind_rows(bind_rows(list_par_pool) )
    return(df)
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


generate_data_alt <- function(data_titre_model, biomarkers, known_exp = NULL) {
    #data_titre_model <- data_sero
    N <- data_titre_model$id %>% unique %>% length  
    N_data <- nrow(data_titre_model)

    titre_true <-  data_titre_model %>% select(all_of(biomarkers)) %>% as.matrix
    times_full <- data_titre_model$time
    id_full <- data_titre_model$id

    initialTitreTime <- data_titre_model %>% group_by(id) %>% filter(time == min(time)) %>% unique %>% .[["time"]]
    endTitreTime <- data_titre_model %>% group_by(id) %>% filter(time == max(time)) %>% unique %>% .[["time"]]

    initialTitreValue <- data_titre_model %>% group_by(id) %>% filter(time == min(time)) %>% unique %>% ungroup %>% select(all_of(biomarkers)) %>% as.matrix
    endTitreValue <- data_titre_model %>% group_by(id) %>% filter(time == max(time)) %>% unique %>% ungroup %>% select(all_of(biomarkers)) %>% as.matrix

    T <- max(endTitreTime)

    # Make titre and times into lists of vectors for each individual
    titre_list <- list()
    titre_list_b <- list()
    times_list <- list()
    for (i in 1:N) {
        data_titre_model_i <- data_titre_model %>% filter(id == i)
        for (b in 1:length(biomarkers)) {
            titre_list_b[[b]] <- data_titre_model_i %>% pull(!!biomarkers[b])
        }
        titre_list[[i]] <- titre_list_b
        times_list[[i]] <- data_titre_model_i %>% pull(time)
    }

    # this is just for the simulated data
    if (!is.null(known_exp)) {
        knownExpVec <- known_exp
    } else {
        knownExpVec <- NA
    }
    cat(endTitreValue)
    data_t <- list(
        N = N,
        T = T,
        N_data = N_data,
        initialTitreValue = initialTitreValue,
        initialTitreTime = initialTitreTime,
        endTitreTime = endTitreTime,
        endTitreValue = endTitreValue,
        titre_full = titre_true,
        times_full = times_full,
        titre_list = titre_list,
        times_list = times_list,
        id_full = id_full,
        knownExpVec = knownExpVec
    )
    data_t
}

#' @brief This function creates data_t, which is a list of data to use in the rjmc model
#' 
#' This function requires a `data_titre_model` data.frame from the user, which is the observed titre data, and a `data_inf_model` data.frame from the user, which is the observed infection data. The function will then create a list of data to use in the rjmc model.
#' @param data_titre_model is a data.frame with column headings ['id', 'titre', 'time'] where 'id' is the individual, 'titre' is the titre value and 'time' is the time of the titre. If not defined correctly it will throw an error.
#' @param data_inf_model is a data.frame with column headings ['id', 'inf_time'] where 'id' is the individual and 'inf_time' is the time of infection. If not defined it will assume no known infection
#' @param known_exp add a known exposure vector to the model. If not defined it will assume no known exposure
#' @return returns a list of data
#' 
generate_data_t <- function(data_titre_model, data_inf_model = NULL, known_exp = NULL) {

    N <- data_titre_model$id %>% unique %>% length  
    N_data <- nrow(data_titre_model)
    titre_true <- data_titre_model$titre
    times_full <- data_titre_model$time
    id_full <- data_titre_model$id
    initialTitreValue <- data_titre_model %>% group_by(id) %>% filter(time == min(time)) %>% .[["titre"]]
    initialTitreTime <- data_titre_model %>% group_by(id) %>% filter(time == min(time)) %>% .[["time"]]
    endTitreTime <- data_titre_model %>% group_by(id) %>% filter(time == max(time)) %>% .[["time"]]
    endTitreValue <- data_titre_model %>% group_by(id) %>% filter(time == max(time)) %>% .[["titre"]]

    # Determine known infections
    data_inf_model_temp <- data.frame(id = 1:N, known = 0, inf_time = -1)

    if (is.null(data_inf_model)){
        ids_known <- vector()
        inftime_known <- vector()
        knownInfsN <- 0
    } else {
        ids_known <- data_inf_model$id
        inftime_known <- data_inf_model$inf_time 
        knownInfsN <- nrow(data_inf_model)
    }

    for (i in seq_len(length(ids_known))) {
        data_inf_model_temp[data_inf_model_temp$id == ids_known[i], "known"] <- 1
        data_inf_model_temp[data_inf_model_temp$id == ids_known[i], "inf_time"] <- inftime_known[i]
    }

    knownInfsVec <- data_inf_model_temp$known
    knownInfsDate <- data_inf_model_temp$inf_time


    if (length(initialTitreValue) > N) {
        stop("ERROR: initialTitreValue wrong size ")
    }
    if (length(initialTitreTime) > N) {
        stop("ERROR: initialTitreTime wrong size ")
    }
    if (length(endTitreTime) > N) {
        stop("ERROR: endTitreTime wrong size ")
    }
    if (!is.null(known_exp)) {
        knownExpVec <- known_exp %>% as.data.frame %>% 
        mutate(i = 1:N) %>% pivot_longer(!i, names_to = "t", values_to = "exp") %>% mutate(t = as.numeric(substr(t, 2, 4))) %>% 
        filter(exp == 1) %>% complete(i = 1:N, fill = list(t = -1, exp = -1)) %>% pull(t)
    } else {
        knownExpVec <- NA
    }

    # Things needed, data for likelihood and sampler
    data_t <- list(
        N = N,
        T = max(endTitreTime),
        N_data = N_data,
        initialTitreValue = initialTitreValue,
        endTitreValue = endTitreValue,
        initialTitreTime = initialTitreTime,
        endTitreTime = endTitreTime,
        titre_full = titre_true,
        times_full = times_full,
        id_full = id_full,

        knownExpVec = knownExpVec,
        knownInfsVec = knownInfsVec,
        knownInfsN = knownInfsN,
        knownInfsDate = knownInfsDate
    )
    data_t
}


#' @brief This function create the rjmc model.
#' 
#' By taking the defined prior distributions, this creates the rjmc model including the support, names of parameters, prior pdf, prior smapling function, the intial conditions for the exposure timings and adds the correlate of protection.
#' 
#' @param ab_ll the likelihood functions define by a function called `evaluateLogLikelihood`
#' @param par_tab a data.frame describing the prior distributions
#' @param cop_func a function for the correlate of protection
#' @return returns a rjmc model
#' 
#' @see createModelRJCMCFull() and generate_data_t()
createModelRJCMCFull <- function(ab_ll, par_tab) {
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
        init_exposure <- c(datalist$knownInfsDate)
        init_exposure
    }

    model_type
}

#' @brief This function adds the exposure prior to the rjmc model.
#' 
#' This function takes two types of priors. Either a functional prior or an empirical prior.
#' 
#' @param model_type rjmc_model created from createModelRJCMCFull()
#' @param data_t data included in the model created from generate_data_t()
#' @param exp_prior a data.frame describing either a functional or empirical prior
#' @param type define a type of exposure prior, `func` for functional and 'empirical' for empirical. If blank it will assume a uniform distribution between 1 and T.
#' @return returns a rjmc model
#' 
#' @see createModelRJCMCFull() and generate_data_t()
addExposurePrior <- function(model_type, data_t, exp_prior, type = NULL) {
    if (is.null(type)) {
        cat("Exposure rate is not defined over the time period. Defaulting to uniform distribution between 1 and ", data_t$T, ". \n")
        model_type$exposureFunctionSample <- function() {
            s <- runif(1, 1, data_t$T)
        }
        model_type$exposureFunctionDensity <- function(jump_i) {
            d <- log(1/data_t$T)
            d
        }
    } else if (type == "func") {
        # Code to check form of exp_prior
        addExposurePrior_checkfunction(exp_prior)

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
    } else if (type == "empirical") {
            
        addExposurePrior_checkempirical(exp_prior, data_t)

        # Code to check form of exp_prior
        model_type$exposureFunctionSample <- function() {
            sample(exp_prior$day, 1, prob = exp_prior$prob)
        }

        model_type$exposureFunctionDensity <- function(jump_i) {
            exp_prior$prob[max(round(jump_i, 0), 1)] %>% log
        }
    } else {
        cat("'type' argument must be either NULL, 'func' or 'empirical'. \n")
    }
    model_type
}


addFunctionTitreExp <- function(model_type, calculate_titre_exp_func = NULL) {

    if (is.null(calculate_titre_exp_func)) {
        cat("Function to calculate titre exposure not defined. Defaulting to titre value at first bleed for each individual \n")
        model_type$calculateTitreExp <- function(params, jump, jump_inf, covariance, datalist) {
            N <- datalist$N
            titre_exp <- vector(mode = "numeric", length = N)
            for (i in 1:N) {
                if (jump[i] == -1) {
                    titre_exp[i] <- -1
                } else {
                    titre_exp[i] <- datalist$initialTitreValue[i]
                }
            }
            titre_exp
        }
    } else{
         model_type$calculateTitreExp <- calculate_titre_exp_func
    }
    model_type
}

rdunif <- function(n, lower_bound, upper_bound) {
  if (!is.numeric(n) || !is.numeric(lower_bound) || !is.numeric(upper_bound)) {
    stop("Arguments must be numeric")
  }
  if (lower_bound >= upper_bound) {
    stop("Lower bound must be less than upper bound")
  }
  
  return(sample(lower_bound:upper_bound, n, replace = TRUE))
}



#safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
#                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
#scales::show_col(safe_colorblind_palette)

clean_simulated_rjmcmc <- function(modelname_sim, obs_er, prob_known, known_exp) {
    modeli <- readRDS(here::here("outputs", "sim_data", modelname_sim, "inputs.RDS"))
    res <- readRDS(file = here::here("outputs", "sim_data", modelname_sim, paste0("sim_data_", obs_er, ".rds")))
    data_titre_model <- res$observed_biomarker_states %>% select(i, t, value) %>% rename(id = i, time = t, titre = value) 
    data_titre_model
}

check_priors <- function(exp_prior_w2, gambia_exp_w2) {
    p1 <- exp_prior_w2 %>% 
        ggplot() + 
            geom_col(aes(x = day, y = prob), color = "red") + 
            theme_bw() + 
            labs(x = "Day of study", y = "Prior exposure probability")
    p2 <- gambia_exp_w2 %>% filter(inf_d_time > -1) %>% 
        ggplot() + 
            geom_histogram(aes(x = inf_d_time)) + 
            theme_bw() + 
            labs(x = "Day of study", y = "Known infection times ")
    p1 / p2 + plot_annotation(tag_levels = "A")
    ggsave(here::here("outputs", "fits", "test", "transvir", "figs", "wave2", "prior_comp.png"))
}


