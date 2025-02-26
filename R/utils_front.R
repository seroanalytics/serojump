extract_long_post <- function(mcmc_out, ...) {
    mcmc_out %>% combine %>% as.data.frame %>%
        dplyr::select(...) %>% 
        pivot_longer(everything(), names_to = "var", values_to = "value")
}
extract_post <- function(mcmc_out, ...) {
    mcmc_out %>% combine %>% as.data.frame %>%
        dplyr::select(...) 
}

#' @title addPrior
#' @description This function adds information about the prior distributions of the parameter in the serojump model. It takes a list of the form: 
#' - par_name: The name of the parameter as as string.
#' - lb: The numerical lower bound of the parameter.
#' - ub: The numerical upper bound of the parameter.
#' - dist: The distribution describing the prior distribution from standard R disributions (e.g. unif, norm)
#' - dist_par1: The first argument of the prior distribution, 
#' - dist_par2: The second argument of the prior distribution, (NA if only one parameter is taken like in exp)
#' @param par_name The name of the parameter as as string.
#' @param lb The numerical lower bound of the parameter.
#' @param ub The numerical upper bound of the parameter.
#' @param dist The name of the parameter as as string.
#' @param dist_par1  The first argument of the prior distribution, 
#' @param dist_par2 The second argument of the prior distribution, (NA if only one parameter is taken like in exp)
#' @return a data.frame of the prior information
#' @export
addPrior <- function(par_name, lb, ub, dist, dist_par1, dist_par2) {
    # Convert the ellipsis arguments to a list
    args_list <- list(par_name, lb, ub, dist, dist_par1, dist_par2)
    names(args_list) <- c("par_name", "lb", "ub", "dist", "dist_par1", "dist_par2")
    df <- as.data.frame(args_list)
    df$dist_par1 <- as.character(df$dist_par1)
    df$dist_par2 <- as.character(df$dist_par2)
    df$part_type <- "prior"

    # Check the distribution exists
    distributions <- ls("package:stats")
    # Filter for distribution functions
    distribution_functions <- distributions[grepl("^d", distributions)]
    if (!paste0("d", dist) %in% distribution_functions) {
        stop(paste0("Error: `dist` = ", dist," does not correspond to a probability density function in the stats package."))
    }
    if (lb >= ub) {
        stop(paste0("Error: `lower bound`, ", lb, ",  is greater than `upper bound`, ", ub))
    }

    aschar_dist <- paste0("r", df$dist)
    my_dist_name <- get(aschar_dist)

    if (df$dist == "unif"){
        if (dist_par1 >= dist_par2) {
            stop("Invalid arguments: for a uniform distribution, lower must be less than upper.")
        }
    }

    if (df$dist != "exp") {
        temp <- do.call(my_dist_name, list(1, as.numeric(df$dist_par1), as.numeric(df$dist_par2)) )
    } else {
        temp <- do.call(my_dist_name, list(1, as.numeric(df$dist_par1)) )
    }
    if (is.nan(temp) | is.na(temp)) {
        stop("Error: cannot sample a random variable, check inputs")
    }
    return(df)
}

addPriorHeir <- function(par_name, lb, ub, dist, dist_par1, dist_par2, dist_sd, dist_sd_par1, dist_sd_par2, dim) {

    
    mean_df <- addPrior(par_name, lb, ub, dist, dist_par1, dist_par2)
    # Convert the ellipsis arguments to a list

    args_list <- list(par_name, dim, dist_sd, dist_sd_par1, dist_sd_par2)
    names(args_list) <- c("par_name", "dim", "dist_sd", "dist_sd_par1", "dist_sd_par2")
    
    list_par_pool <- list(mode = "list", length = args_list$dim)
    for (i in 1:args_list$dim) {
        args_list_temp <- list(par_name = paste0("z_", args_list$par_name, "_", i), lb = -10, ub = 10, dist = "norm",
            dist_par1 = "0", dist_par2 = "1")
        list_par_pool[[i]] <- as.data.frame(args_list_temp)
        list_par_pool[[i]]$part_type <- "hprior"
    }
    list_par_reg <- list(paste0("sigma_", args_list$par_name), 0, 5, args_list$dist_sd, args_list$dist_sd_par1, args_list$dist_sd_par2)
    names(list_par_reg) <-  c("par_name", "lb", "ub", "dist", "dist_par1", "dist_par2")
    list_par_reg$dist_par1 <- as.character(list_par_reg$dist_par1)
    list_par_reg$dist_par2 <- as.character(list_par_reg$dist_par2)
    list_par_reg$part_type <- "hprior"

    df <- mean_df %>% bind_rows(bind_rows(list_par_pool)) %>% bind_rows( as.data.frame(list_par_reg))

    
    # Check the distribution exists
    distributions <- ls("package:stats")
    # Filter for distribution functions
    distribution_functions <- distributions[grepl("^d", distributions)]
    if (!paste0("d", list_par_reg$dist) %in% distribution_functions) {
        stop(paste0("Error: `dist_sd` = ", list_par_reg$dist," does not correspond to a probability density function in the stats package."))
    }

    aschar_dist <- paste0("r", list_par_reg$dist)
    my_dist_name <- get(aschar_dist)

    if (list_par_reg$dist != "exp") {
        temp <- do.call(my_dist_name, list(1, as.numeric(args_list$dist_sd_par1), as.numeric(args_list$dist_sd_par2)) )
    } else {
        temp <- do.call(my_dist_name, list(1, as.numeric(args_list$dist_sd_par1)) )
    }
    if (is.nan(temp) | is.na(temp)) {
        stop("Error: cannot sample a random variable, check inputs")
    }

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


# Can add heirarchical priors with the functions below 
# add_par_pool_df_non_centered("boost_naive_pvnt", length = 5, 0, "boost_naive_sigma_pvnt", "exp", 5, NA)


get_sample_non_centered <- function(par_tab, seed = -1) {
    if (seed > -1) {
        set.seed(seed)
    }
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

#' @title check_input_data
#' @details This function checks the data_sero input for the serojump model.
#' @param data_sero A data.frame with columns id, time and IgG
#' @param known_inf A data.frame with columns id, time and exposure_type
#' @return A warning if there are single entries
#' @export 
check_input_data  <- function(data_sero, known_inf) {

    check_sero_no_single_entries(data_sero)
    check_sero_timings(data_sero)
}

#' @title check_no_single_entries
#' @details This function checks that the data_sero has more than one observation per person.
#' @param data_sero A data.frame with columns id, time and IgG
#' @return A warning if there are single entries
#'
#' @note This function can be used to check the data_sero before running the serojump model.
#' 
#' @export
check_sero_no_single_entries <- function(data_sero) {
      data_sero_single_obs <- data_sero %>% group_by(id) %>% mutate(r = row_number()) %>% filter(max(r) == 1)
    if(nrow(data_sero_single_obs) != 0) {
        stop("ERROR in data_sero: must have more than one observation be person, check ids: ", 
        paste(data_sero_single_obs$id, collapse = ", "), 
        ". Please remove")
    } else {
        cat("No single entries in data_sero!, \n")
    }
}

#' @title check_sero_timings
#' @details This function checks that the data_sero has more than one observation per person.
#' @param data_sero A data.frame with columns id, time and IgG
#' @return A warning if there are single entries
#' @note This function can be used to check the data_sero before running the serojump model.
#' 
#' @export
check_sero_timings <- function(data_sero) {

    initialTitreTime <- data_sero %>% group_by(id) %>% filter(time == min(time)) %>% unique %>% .[["time"]]
    endTitreTime <- data_sero %>% group_by(id) %>% filter(time == max(time)) %>% unique %>% .[["time"]]

    small_time <- vector()
    for (i in 1:length(initialTitreTime)) {
        if ((endTitreTime[i] - 15) - initialTitreTime[i] < 0) {
            small_time <- c(small_time, i)
        }   
    }
    if (length(small_time) != 0) {
        stop("There are individuals which are in the study for less than 15 days, ids: ",  paste(small_time, collapse = ", "), ". It is hard to perform inference on these people, so we recommend you remove them.")
    } else {
        cat("No individuals with less than 15 days in the study! \n")
    }
}   


#know_inf_df_i <- data.frame(
#    id = c(1, 2, 3, 4),
#    time = c(-1, 40, 30, -1),
#    start = c(1, 2, 3, 4),
#    type = c("in_sample", "in_sample", "in_sample", "in_sample"),
#    end = c(50, 70, 80, 90)
#)
#T <- 100
#N <- 4
#exp_prior_list_i <- lapply(1:N, function(x) rep(1, 100))


#update_exp_prior(exp_prior_list_i, know_inf_df_i, T, N)

#known_exp <- known_inf_w2
#data_sero <- sero_data_w2
#exposureType <- "vac"
#initialTitreTime <- data_sero %>% group_by(id) %>% filter(time == min(time)) %>% unique %>% .[["time"]]
#endTitreTime <- data_sero %>% group_by(id) %>% filter(time == max(time)) %>% unique %>% .[["time"]]

##know_inf_df_i <- check_oos_sample_df(known_exp, exposureType, initialTitreTime, endTitreTime)
#update_exp_prior(exp_prior_list, know_inf_df_i, T, N)

check_oos_sample_df <- function(known_exp, exposureType, initialTitreTime, endTitreTime) {

    N <- length(endTitreTime)
    known_exp %>% filter(exposure_type == exposureType) %>%
        complete(id = 1:N, fill = list(time = -1)) %>%
        mutate(start = initialTitreTime, end = endTitreTime) %>%
        mutate(type = case_when(time >= start & time <= end ~ "in_sample", TRUE ~ "oos_sample"))
}

update_exp_prior <- function(exp_prior_list_par, know_inf_df_i, T, N) {

    for (i in 1:N) {
        start_t <- know_inf_df_i[i, ]$start
        end_t <- know_inf_df_i[i, ]$end
        exp_prior_list_par[[i]][seq_len(start_t)] <- 0
        exp_prior_list_par[[i]][(end_t - 7):T] <- 0

        if(know_inf_df_i[i, ]$time != -1 && know_inf_df_i[i, ]$type == "in_sample") {
            timeexp <- know_inf_df_i[i, ]$time
            exp_prior_list_par[[i]][(timeexp - 14):(timeexp + 14)] <- 0
        } 
        # normalise 
        exp_prior_list_par[[i]] <- exp_prior_list_par[[i]] / sum(exp_prior_list_par[[i]])
    }
    exp_prior_list_par
}

#' @title check_exposures_times
#' @details This function checks that the data_sero has exposure times within the range the bleeds taken for an individual
#' @param data_sero A data.frame with columns id, time and `biomarkers`
#' @param known_exp A data.frame with columns id, time and exposure_type
#' @param exposure_types A vector of the exposure types
#' @param fitted_exp A string of the exposure type that is fitted
#' @param exposurePriorTime A data.frame describing either a functional or empirical prior
#' @param exposurePriorTimeType define a type of exposure prior, `func` for functional and 'empirical' for empirical. If blank it will assume a uniform distribution between 1 and T.
#' @return A warning if there are values outside this range
#' @note This function can be used to check the data_sero before running the serojump model.
#' 
#' @export
check_exposures_times <- function(data_sero, known_exp, exposure_types, fitted_exp, exposurePriorTime = NULL, exposurePriorTimeType = NULL) {


    initialTitreTime <- data_sero %>% group_by(id) %>% filter(time == min(time)) %>% unique %>% .[["time"]]
    endTitreTime <- data_sero %>% group_by(id) %>% filter(time == max(time)) %>% unique %>% .[["time"]]
    T <- max(endTitreTime)
    N <- length(initialTitreTime)

    known_exp_i <- known_exp %>% rename(exp_time = time) 

    outside_obs <- vector()
    inside_obs <- vector()
    key_exposure <- vector()
    
    exp_prior <- convertExpPriorEmpirical(exposurePriorTime, T, exposurePriorTimeType)
    exp_prior_list <- lapply(1:N, function(x) exp_prior$prob)
    
    T <- length(exp_prior$prob)
    know_inf_df_i <- list()
    j <- 1
    for (exposureType in exposure_types) { 
            # Get the known infections
        know_inf_df_i[[j]] <- check_oos_sample_df(known_exp, exposureType, initialTitreTime, endTitreTime)
        exp_prior_list <- update_exp_prior(exp_prior_list, know_inf_df_i[[j]], T, N)
        j <- j + 1  
    }
    know_inf_df <- know_inf_df_i %>% bind_rows

    # Get information on people 
    outside_obs <- know_inf_df %>% filter(!is.na(exposure_type), type == "oos_sample") %>% pull(id)
    for (i in 1:N) {
        if(sum(exp_prior_list[[i]]) == 0 ) {
            inside_obs <- c(inside_obs, i)
        }
    }    

    if (length(outside_obs) != 0) {
        cat("There are individuals which have exposure before their first bleed or before their last bleed, ids: ",  paste(outside_obs, collapse = ", "), ". It is hard to perform inference on these people, so we recommend you remove the exposures from the known inf, but, if included, the model ignores these exposures.\n")
    }
    if (length(inside_obs) != 0) {
        cat("There are individuals which no possible exposure times after the known infection assumptions have been included, ids:",  paste(outside_obs, collapse = ", "), ". These people cannot be exposure to fitted exposured. \n")
    }  
    
    if (length(outside_obs) == 0 && length(inside_obs) == 0) {
        cat("No individuals with exposure times outside the bleed times! All inf priors look good, \n")
    }
   # known_exp$key_exposure <- key_exposure
    #known_exp
}


generate_data <- function(data_sero, biomarkers, known_exp_bool = NULL) {
    
    # Check the entries!    
    check_sero_no_single_entries(data_sero)
    check_sero_timings(data_sero)

    data_titre_model <- data_sero
    N <- data_titre_model$id %>% unique %>% length  
    N_data <- nrow(data_titre_model)

    titre_true <-  data_titre_model %>% select(all_of(biomarkers)) %>% as.matrix
    times_full <- data_titre_model$time
    id_full <- data_titre_model$id
    pid_full <- data_titre_model$pid

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
    if (!is.null(known_exp_bool)) {
        knownExpVec <- rep(-1, N)
    } else {
        knownExpVec <- NA
    }

    max_titre <- titre_true %>% apply(2, max)

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
        max_titre = max_titre,
        pid_full = pid_full,
        id_full = id_full,
        knownExpVec = knownExpVec
    )
    data_t
}


convertExpPriorEmpirical <- function(exposurePriorTime, T_max, exposurePriorTimeType = NULL) {

    if (is.null(exposurePriorTimeType)) {
        cat("Exposure rate is not defined over the time period. Defaulting to uniform distribution between 1 and ", T_max, ". \n")

        exp_prior <- data.frame(
            time = 1:T_max,
            prob = dunif(1:T_max, 1, T_max)
        )
    } else if (exposurePriorTimeType == "func") {
        addExposurePrior_checkfunction(exposurePriorTime)
        dist_name <- paste0("d", exposurePriorTime[1, 3])
        my_dist_name <- get(dist_name)

        s <- do.call(my_dist_name, list(1:T_max, as.numeric(exposurePriorTime[1, 4]),  as.numeric(exposurePriorTime[1, 5])) )

        exp_prior <- data.frame(
            time = 1:T_max,
            prob = s
        )
    } else if (exposurePriorTimeType == "empirical") {
        exp_prior <- exposurePriorTime
        addExposurePrior_checkempirical(exp_prior, T_max)
    } else {
        cat("'exposurePriorTimeType' argument must be either NULL, 'func' or 'empirical'. \n")
    }
    exp_prior
}



calculateIndExposure <- function(model_type, data_t, exp_prior_i, type = NULL) {

    if (is.null(type)) {
        cat("Exposure rate is not defined over the time period. Defaulting to uniform distribution between 1 and ", data_t$T, ". \n")

        exp_prior <- data.frame(
            time = 1:data_t$T,
            prob = dunif(1:data_t$T, 1, data_t$T)
        )
    } else if (type == "func") {
        addExposurePrior_checkfunction(exp_prior_i)
        dist_name <- paste0("d", exp_prior_i[1, 3])
        my_dist_name <- get(dist_name)

        s <- do.call(my_dist_name, list(1:data_t$T, as.numeric(exp_prior_i[1, 4]),  as.numeric(exp_prior_i[1, 5])) )

        exp_prior <- data.frame(
            time = 1:data_t$T,
            prob = s
        )
    } else if (type == "empirical") {
        exp_prior <- exp_prior_i
        addExposurePrior_checkempirical(exp_prior, data_t)
    } else {
        cat("'type' argument must be either NULL, 'func' or 'empirical'. \n")
    }

    # Get the inidividual level exposure probabilities 
    T <- exp_prior$time %>% length
    UnknownInfsVec <- rep(-1, data_t$N)
    exp_list <- list()
    for (i in 1:data_t$N) {
        exp_i <- exp_prior$prob
        start_t <- data_t$initialTitreTime[i]
        end_t <- data_t$endTitreTime[i]
        exp_i[seq_len(start_t + 7)] <- 0
        exp_i[(end_t - 7):T] <- 0

        # Impossible times due to known infection
        exposureInfo <- model_type$infoModel$exposureInfo
        known_vec_ind <- vector()
        for (j in 1:length(exposureInfo)) {
            known_exp <- exposureInfo[[j]]$known_inf[i]
            if (known_exp > -1) {
                known_vec_ind[j] <- known_exp
                exp_i[(known_exp):(known_exp + 21)] <- 0
            }
        }

        # Add other known exposures into the mix 
        if (data_t$knownInfsVec[i] == 1) {
            exp_i[data_t$knownInfsTimeVec[i]] <- 1
        }

        if ( sum(exp_i) == 0) {
            exp_list[[i]] <- rep(0, T)
            UnknownInfsVec[i] <- 1
        } else {
            exp_list[[i]] <- exp_i / sum(exp_i)
        }
    }

    data_t$UnknownInfsVec <- UnknownInfsVec
    data_t$knownInfsN <- sum(data_t$knownInfsVec)

    data_t$exp_list <- exp_list
    data_t$exp_prior <- exp_prior

    data_t
}

#' @title This function adds the exposure prior to the rjmc model.
#' 
#' @description This function takes two types of priors. Either a functional prior or an empirical prior.
#' 
#' @param model_type rjmc_model created from createModelRJCMCFull()
#' @param data_t data included in the model created from generate_data_t()
#' @param exp_prior a data.frame describing either a functional or empirical prior
#' @param type define a type of exposure prior, `func` for functional and 'empirical' for empirical. If blank it will assume a uniform distribution between 1 and T.
#' @return returns a rjmc model
#' 
#' @seealso createModelRJCMCFull() and generate_data_t()
#' @importFrom stats runif
addExposurePrior <- function(model_type, data_t, exp_prior, type = NULL) {

    #model_type <- modelSeroJump
    #data_t
    ##exp_prior <- modeldefinition$exposurePrior
    #type <- modeldefinition$exposurePriorType


    if (is.null(type)) {
        cat("Exposure rate is not defined over the time period. Defaulting to uniform distribution between 1 and ", data_t$T, ". \n")
        model_type$exposureFunctionSample <- function() {
            s <- runif(1, 1, data_t$T)
        }
        model_type$exposureFunctionDensity <- function(jump_i,  i) {
            d <- log(1/data_t$T)
            d
        }
    } else if (type == "func") {
        # Code to check form of exp_prior
        addExposurePrior_checkfunction(exp_prior)

        model_type$exposureFunctionSample <- function(i) {
            dist_name <- paste0("r", exp_prior[1, 3])
            my_dist_name <- get(dist_name)
            s <- do.call(my_dist_name, list(1, as.numeric(exp_prior[1, 4]),  as.numeric(exp_prior[1, 5])) )
            s
        }

        model_type$exposureFunctionDensity <- function(jump_i, i) {
            dist_name <- paste0("d", exp_prior[1, 3])
            my_dist_name <- get(dist_name)
            d <- do.call(my_dist_name, list(as.numeric(jump_i), as.numeric(exp_prior[1, 4]), as.numeric(exp_prior[1, 5])) )
            d
        }
    } else if (type == "empirical") {
            
        addExposurePrior_checkempirical(exp_prior, data_t)

        T <- exp_prior$time %>% length
        exp_prior$prob
        exp_list <- list()
        for (i in 1:data_t$N) {
            exp_i <- exp_prior$prob
            start_t <- data_t$initialTitreTime[i]
            end_t <- data_t$endTitreTime[i]
            exp_i[seq_len(start_t + 7)] <- 0
            exp_i[(end_t - 7):T] <- 0

            # Impossible times 
            known_exp <- data_t$knownInfsTimeVec[i]
            if (known_exp > -1) {
                exp_i[known_exp:(known_exp + 21)] <- 0
            }
            exp_list[[i]] <- exp_i / sum(exp_i)
        }

        # Code to check form of exp_prior
        model_type$exposureFunctionSample <- function(i) {
            sample(exp_list[[i]]$time, 1, prob = exp_prior[[i]]$prob) 
        }

        model_type$exposureFunctionDensity <- function(jump_i, i) {
            exp_prior[[i]]$prob[max(round(jump_i, 0), 1)] %>% log
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

clean_simulated_rjmcmc <- function(modelname_sim, obs_er) {
    modeli <- readRDS(here::here("outputs", "sim_data", modelname_sim, "inputs.RDS"))
    res <- readRDS(file = here::here("outputs", "sim_data", modelname_sim, paste0("sim_data_", obs_er, ".rds")))
    data_titre_model <- res$observed_biomarker_states %>% select(i, t, value) %>% rename(id = i, time = t, titre = value) 
    data_titre_model
}
