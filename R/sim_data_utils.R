antibody_model_exp_asymp <-  function(i, t1, b, immune_histories, biomarker_states, kinetics_parameters, biomarker_map, biomarker_quantity_init, ...){
    biomarker_quantity <- biomarker_quantity_init
    
    ## Get kinetics parameters for this individual
    if(!is.null(kinetics_parameters[[i]])){
    tmp_kinetics_parameters <- kinetics_parameters[[i]]
    } else {
    return(biomarker_quantity)
    }
    
    ## Get only exposures relevant to this biomarker ID and time
    tmp_kinetics_parameters <- tmp_kinetics_parameters[tmp_kinetics_parameters$b == b & tmp_kinetics_parameters$t <= t1,]
    
    ## Only continue if there are relevant exposures to calculate kinetics for
    if(nrow(tmp_kinetics_parameters) > 0){
    a <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "a",]$realized_value ## Boosts
    b <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "b",]$realized_value ## Boosts
    c <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "c",]$realized_value ## Waning rate
    alpha <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "biomarker_ceiling_gradient",]$realized_value ## Waning rate


    t_infs <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "a",]$t ## Time of infections
    ## Sum up contribution of each boost, with waning
    for(j in seq_along(t_infs)){
        if ((t1-t_infs[j]) < 14) {
            biomarker_quantity <- biomarker_quantity + max(0, (1 - alpha * biomarker_quantity)) * (log(exp(a) +  exp(c)) * (t1-t_infs[j]) / 14);
        } else {
            biomarker_quantity <- biomarker_quantity + max(0, (1 - alpha * biomarker_quantity)) * (log(exp(a) * exp(-b/10 * ((t1-t_infs[j]) - 14)) + exp(c)));
        }
    }
    }
    biomarker_quantity
}



runserosim_edit <- function(
    ## SIMULATION SETTINGS
    simulation_settings, ## List of parameters governing the simulation settings
    demography=NULL, ## tibble of demographic information for each individual
    observation_times=NULL, ## tibble of observation times and biomarkers measured for each individual
    foe_pars, ## any object (usually a 3D array) giving force of infection for each exposure ID, groups and time
    biomarker_map, ## Object determining relationship between exposure IDs and biomarkers
    model_pars,
    
    ## FUNCTIONS
    exposure_model, ## calculates the probability of infection given the FOE array, foe_pars
    immunity_model, ## function determining probability of infection conditional on foe_pars and individual's immune state
    antibody_model, ## function determining biomarker quantity, commonly antibody state, as a function of exposure history and kinetics parameters (model_pars)
    observation_model, ## function generating observed biomarker quantity as a function of latent biomarker quantities and model_pars
    draw_parameters, ## function to simulate biomarker (e.g, antibody)  kinetics parameters
    
    ## Pre-specified parameters/events (optional)
    immune_histories_fixed=NULL,
    start_titre = NULL,
    ## UPDATE MESSAGE
    VERBOSE=NULL,
    attempt_precomputation=TRUE,
    parallel=FALSE,
    n_cores=4,
    parallel_packages=NULL,
    ...
                    ){
    
    ## Simulation settings
    t_start <- simulation_settings[["t_start"]]
    t_end <- simulation_settings[["t_end"]]
    times <- seq(simulation_settings[["t_start"]],simulation_settings[["t_end"]],1)
    
    ## Extract key demographic information
    indivs <- unique(demography$i)
    N <- length(indivs)
    
    ## Note "birth" refers to first time point in the population and "removal" refers to time point of removal from population
    if(!("birth" %in% colnames(demography))) {
        demography$birth <- t_start
    } 
    birth_times <- demography %>% dplyr::select(i, birth) %>% distinct()
    if(!("removal" %in% colnames(demography))) {
        demography$removal <- t_end
    } 
    removal_times <- demography %>% dplyr::select(i, removal) %>% distinct() 
    
    ## If no groups information provided, assume 1 group
    ## ...
    if(!("group" %in% colnames(demography))) {
        demography$group <- 1
    }
    groups <- demography %>% dplyr::select(group) %>% distinct() %>% pull(group)
    N_groups <- length(groups)
    
    ## Make sure demography is arranged by i and times (if there is a time element). 
    ## This is important for subsetting below. 
    ## Also pre-compute which group the individual is in for each time point to speed up
    ## indexing below
    if("times" %in% colnames(demography)){
        use_time <- TRUE
        demography <- demography %>% arrange(i, times)
        ## Create groups matrix
        all_groups <- matrix(NA, nrow=N,ncol=length(times))
        indices <- convert_indices_matrix_to_vector(demography$i, demography$times, N)
        all_groups[indices] <- demography$group
        if(any(is.na(all_groups))) message(cat("Warning - could not find group membership for all individuals at all times. This may lead to an error or unexpected output"))
    } else {
        use_time <- FALSE
        demography <- demography %>% arrange(i)
        all_groups <- demography$group
        all_groups <- matrix(rep(all_groups,each=length(times)),nrow=N,ncol=length(times),byrow=TRUE)
    }
    ## Extract information on number of exposure types
    exposure_ids <- unique(biomarker_map$exposure_id)
    biomarker_ids <- unique(biomarker_map$biomarker_id)
    N_exposure_ids <- length(exposure_ids)
    N_biomarker_ids <- length(biomarker_ids)
    demography<-data.table::data.table(demography)
    
    ########################################################################
    ## Checking if pre-compuation of exposure probabilities is possible
    if(attempt_precomputation){
        if(!is.null(VERBOSE)) message(cat("Checking for possible pre-computation to save time...\n"))
        
        precomputations <- precomputation_checks(N, times,exposure_ids, groups,
                                                exposure_model, foe_pars, demography, 
                                                VERBOSE, ...)
        if(precomputations$flag == TRUE){
            foe_pars <- precomputations$foe
            exposure_model <- exposure_model_indiv_fixed
        } 
        if(!is.null(VERBOSE)) message(cat("Using pre-computed exposure probabilities\n"))
        
    } 
    ## If successful precomputation, change the exposure model to just read directly from foe_pars
    ########################################################################
    
    ## Parallelize here
    serosim_internal <- function(tmp_indivs, ...){
        ## Create empty arrays to store immune histories
        immune_histories <- array(NA, dim=c(N, length(times), N_exposure_ids))
        exposure_probabilities <- array(NA, dim=c(N, length(times), N_exposure_ids))
        exposure_force <-array(NA, dim=c(N, length(times), N_exposure_ids))
        biomarker_states <- array(0, dim=c(N, length(times), N_biomarker_ids))
        kinetics_parameters <- vector(mode="list",length=N)
        
        ## Merge in any pre-specified immune history information
        ## ...
        if(!is.null(immune_histories_fixed)){
            immune_histories <- ifelse(!is.na(immune_histories_fixed), immune_histories_fixed, immune_histories) 
        }
        
        if(!is.null(VERBOSE)) message(cat("Beginning simulation\n"))
    
        ## For each individual
        for(i in tmp_indivs){
            ## Print update message
            update(VERBOSE,i)
            ## Pull birth time for this individual
            birth_time <- birth_times$birth[i]
            removal_time <- ifelse(is.na(removal_times$removal[i]), simulation_settings[["t_end"]], removal_times$removal[i])
            
            ## Only consider times that the individual was alive for
            simulation_times_tmp <- times[times >= birth_time & times <= removal_time]
            ## Go through all times relevant to this individual
            for(t in simulation_times_tmp){
                ## Pull group for this individual at this time 
                g <- all_groups[i, t]
                ## Work out antibody state for each biomarker
                ## The reason we nest this at the same level as the immune history generation is
                ## that immune histories may be conditional on antibody state
                if(immune_histories[i, t, 1] == 1){
                    kinetics_parameters[[i]] <- bind_rows(kinetics_parameters[[i]],
                                                    draw_parameters(i, t, x, demography, biomarker_states, model_pars, ...))
                }
                if (t == 1) {
                        biomarker_states[i, t, 1] <- start_titre[i]
                } else {
                for(b in biomarker_ids){
                    biomarker_states[i,t,b] <- antibody_model(i, t, b, immune_histories, 
                                                                biomarker_states, kinetics_parameters, biomarker_map, start_titre[i], ...)

                }
                }
                
                ## Work out exposure result for each exposure ID
                for(x in exposure_ids){
                    ## Only update if immune history entry is NA here. If not NA, then pre-specified
                    if(is.na(immune_histories[i,t,x])){
                    
                        ## What is the probability that exposure occurred?
            

                        prob_exposed <- exposure_model(i, t, x, g, foe_pars, demography, ...)

                        ## If an exposure event occurred, what's the probability 
                        ## of successful exposure?
                        prob_success <- immunity_model(i, t, x, immune_histories, 
                                                    biomarker_states, demography, 
                                                    biomarker_map, model_pars, ...)
                        ## Randomly assign success of exposure event based on immune state
                        successful_exposure <- as.integer(runif(1) < prob_success*prob_exposed)
                        successful_exposure <- successful_exposure[1]
                        
                        ## Simulate kinetics parameters for this exposure event
                        ## Each successful exposure event will create a tibble with parameters
                        ## for this event, drawn from information given in model_pars

                        immune_histories[i,t,x] <- successful_exposure
                        exposure_probabilities[i,t,x] <- prob_success*prob_exposed
                        exposure_force[i,t,x] <- prob_exposed
                        if(successful_exposure == 1){
                            for(b in biomarker_ids){
                                biomarker_states[i, t, b] <- antibody_model(i, t, b, immune_histories, 
                                                                        biomarker_states, kinetics_parameters, biomarker_map, start_titre[i],...)
                            }
                        }
                    }
                }
            }
        }
        return(list(array(biomarker_states[tmp_indivs,,],dim=c(length(tmp_indivs), length(times),N_biomarker_ids)), 
                    kinetics_parameters[tmp_indivs], 
                    immune_histories[tmp_indivs,,], 
                    exposure_probabilities[tmp_indivs,,], 
                    exposure_force[tmp_indivs,,]))
    }
    ## Run either the entire simulation in one go, or split into n_cores jobs and run in parallel
    if(!parallel){
        res <- serosim_internal(indivs,...)
        biomarker_states <- res[[1]]
        kinetics_parameters <- res[[2]]
        immune_histories <- res[[3]]
        exposure_probabilities <- res[[4]]
        exposure_force <- res[[5]]
    } else {

        ## Set up socket cluster to split the population across n_cores processes
        cluster <- makeCluster(n_cores) 
        registerDoParallel(cluster)
        
        if(!is.null(VERBOSE)) message(cat("Running jobs in parallel across ", n_cores, " cores\n"))
        
        n_indivs_per_block <- ceiling(length(indivs)/n_cores)
        indiv_blocks <- split(indivs, ceiling(seq_along(indivs)/n_indivs_per_block))
        
        ## Submit job. Note this is likely where any object or package exports may become an issue.
        res <- foreach(block = 1:n_cores, .packages=c("abind",parallel_packages)) %dopar% {
        serosim_internal(indiv_blocks[[block]],...)
        }
        
        ## Close socket cluster
        stopCluster(cluster)
        
        ## Combine outputs from each process. This could probably be improved without requiring abind.
        biomarker_states <- do.call("abind",args=list(lapply(res, function(x) x[[1]]), along=1))
        immune_histories <- do.call("abind",args=list(lapply(res, function(x) x[[3]]), along=1))
        exposure_probabilities <- do.call("abind",args=list(lapply(res, function(x) x[[4]]), along=1))
        exposure_force <- do.call("abind",args=list(lapply(res, function(x) x[[5]]), along=1))
        kinetics_parameters <- do.call("bind_rows",lapply(res, function(x) x[[2]]))
    }
    
    if(!is.null(VERBOSE)) message(cat("Simulation complete! Cleaning up...\n"))
    
    all_kinetics_parameters <- do.call("bind_rows", kinetics_parameters)
    
    ## Reshape antibody states
    biomarker_states <- reshape2::melt(biomarker_states)
    colnames(biomarker_states) <- c("i","t","b","value")
    biomarker_states <- biomarker_states %>% arrange(i, t, b)

    ## Reshape immune histories
    immune_histories_long <- NULL
    cat(str(immune_histories))
    cat(str(reshape2::melt(immune_histories)))

    if(sum(immune_histories, na.rm = TRUE) > 0){
        immune_histories_long <- reshape2::melt(immune_histories)
        colnames(immune_histories_long) <- c("i","t","value")
        # immune_histories_long <- immune_histories_long %>% filter(value != 0) %>% select(-value)
        immune_histories_long <- immune_histories_long %>% arrange(i, t)
    }

    ## Reshape probability of a successful exposure event
    exposure_probabilities_long <- reshape2::melt(exposure_probabilities)
    colnames(exposure_probabilities_long) <- c("i","t","value")
    exposure_probabilities_long <- exposure_probabilities_long %>% arrange(i, t)

    ## Reshape exposure probabilities
    exposure_force_long <- reshape2::melt(exposure_force)
    colnames(exposure_force_long) <- c("i","t","value")
    exposure_force_long <- exposure_force_long %>% arrange(i, t)

    ## Observation process
    if(!is.null(observation_times)){
        observed_biomarker_states <- observation_model(left_join(observation_times,biomarker_states), model_pars, ...)
    } else {
        observed_biomarker_states <- observation_model(biomarker_states, model_pars, ...)
    }
    
    ## Remove latent states before individual was born and after they left the study
    biomarker_states <- biomarker_states %>% 
        left_join(demography %>% dplyr::select(c(i, birth,removal)) %>% distinct()) %>%
        mutate(value=ifelse(t >= birth & t <= removal, value, NA))%>% 
        dplyr::select(-c(birth,removal))

    ## Remove observations before individual was born and after they left the study
    observed_biomarker_states <- observed_biomarker_states %>% 
        left_join(demography %>% dplyr::select(c(i, birth,removal)) %>% distinct()) %>%
        mutate(value=ifelse(t >= birth & t <= removal, value, NA))%>% 
        mutate(observed=ifelse(t >= birth & t <= removal, observed, NA))%>% 
        dplyr::select(-c(birth,removal))
    
    return(list("immune_histories"=immune_histories,
                "start_titre" = start_titre,
                "immune_histories_long"=immune_histories_long,
                "exposure_probabilities"=exposure_probabilities,
                "exposure_probabilities_long"=exposure_probabilities_long,
                "exposure_force"=exposure_force,
                "exposure_force_long"=exposure_force_long,
                "biomarker_states"=biomarker_states,
                "observed_biomarker_states"=observed_biomarker_states,
                "kinetics_parameters"=all_kinetics_parameters,
                "demography"=demography))
}


draw_parameters_random_fx_biomarker_dep <- function(i, t, x, demography, biomarker_states, model_pars, ...){
  ## Filter for only exposure stimulated 
    model_pars_tmp <- as.data.frame(model_pars[model_pars$exposure_id == x & !is.na(model_pars$exposure_id),])
    pars <- numeric(nrow(model_pars_tmp))
  realized <- numeric(nrow(model_pars_tmp))
  par_names <- character(nrow(model_pars_tmp))
  ## For each parameter; randomly sample from the distribution given the mean and sd 
  for(par_name in seq_len(nrow(model_pars_tmp))){
      if(is.na(model_pars_tmp$distribution[par_name]) | model_pars_tmp$distribution[par_name]==""){
          pars[par_name] <- model_pars_tmp$mean[par_name]
          par_names[par_name] <- model_pars_tmp$name[par_name]
      } else if(model_pars_tmp$distribution[par_name] == "log-normal"){ #Convert the normal distributions to log-normal distributions 
      
      ## Create functions to convert normal distributions to log-normal distributions
      pars[par_name] <- rlnorm(1, normal_to_lognormal_mean(model_pars_tmp$mean[par_name], model_pars_tmp$sd[par_name]), normal_to_lognormal_sd(model_pars_tmp$mean[par_name], model_pars_tmp$sd[par_name]))
      par_names[par_name] <- model_pars_tmp$name[par_name]
    } else if(model_pars_tmp$distribution[par_name]=="normal"){
      pars[par_name] <- rnorm(1, model_pars_tmp$mean[par_name], model_pars_tmp$sd[par_name])
      par_names[par_name] <- model_pars_tmp$name[par_name]
    } else {
        pars[par_name] <- model_pars_tmp$mean[par_name]
        par_names[par_name] <- model_pars_tmp$name[par_name]
    }
    if(par_names[par_name] %in% c("boost_short","boost_long","boost")){
      ## Pull out all biomarker
      biomarker<-model_pars_tmp$biomarker_id[par_name]
      t1<-t-1
      biomarker_threshold <- min(biomarker_states[i,t1,biomarker], model_pars_tmp[model_pars_tmp$name=="biomarker_ceiling_threshold" & model_pars_tmp$biomarker_id==biomarker, "mean"])
      realized[par_name] <- pars[par_name]*(1-model_pars_tmp[model_pars_tmp$name=="biomarker_ceiling_gradient" & model_pars_tmp$biomarker_id==biomarker, "mean"]*biomarker_threshold)
    }
    if(!(par_names[par_name] %in% c("boost_short","boost_long","boost"))){
      ## Non-boost realized parameters don't get affected by biomarker quantity ceiling
      realized[par_name] <- pars[par_name]
    }
  }
  all_pars <- tibble(i=rep(i,nrow(model_pars_tmp)),t=rep(t,nrow(model_pars_tmp)), x=rep(x,nrow(model_pars_tmp)), b=model_pars_tmp$biomarker_id, name=par_names, value=pars, realized_value=realized) 
  return(all_pars)
}




clean_simulated_rjmcmc_alt <- function(modelname_sim, obs_er, prob_known, known_exp = FALSE) {

    modeli <- readRDS(here::here("outputs", "sim_data", modelname_sim, "inputs.RDS"))
    res <- readRDS(file = here::here("outputs", "sim_data", modelname_sim, paste0("sim_data_", obs_er, ".rds")))

    name <- modeli$name
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