#' @title addObservationalModel
#' @description This function adds an observational model to the model definition.
#' @param biomarker The name of the biomarker.
#' @param pars The parameters of the model.
#' @param logLikelihood The log likelihood function.
#' @return A list with the biomarker name, the parameters and the log likelihood function.
#' @export
addObservationalModel <- function(biomarker, pars, logLikelihood) {
    if (is.null(biomarker) || is.null(pars) || is.null(logLikelihood) ) {
        stop("One or more arguments of `addObservationalModel` are NULL.")
    }
    list(
        biomarker = biomarker,
        pars = pars,
        logLikelihood = logLikelihood
    )
}

#' @title addAbkineticsModel
#' @description This function adds an antibody kinetics model to the model definition.
#' @param id The name of the biomarker.
#' @param biomarker The name of the biomarker.
#' @param exposureType The name of the exposure type.
#' @param pars The parameters of the model.
#' @param funcForm The antibody kinetics function.
#' @return A list with the biomarker name, the exposure name, whether the exposure is inferred, the parameters and the antibody kinetics function.
#' @export
addAbkineticsModel <- function(id, biomarker, exposureType, pars, funcForm) {
    if (is.null(id) || is.null(biomarker) || is.null(exposureType) || is.null(pars) || is.null(funcForm)) {
        stop("One or more arguments of `addAbkineticsModel` are NULL.")
    }
    list(
        id = id,
        biomarker = biomarker,
        exposureType = exposureType,
        pars = pars,
        funcForm = funcForm,
        hierFlag = FALSE
    )
}

#' @title addAbkineticsModelHier
#' @description This function adds an antibody kinetics model to the model definition.
#' @param id The name of the biomarker.
#' @param biomarker The name of the biomarker.
#' @param exposureType The name of the exposure type.
#' @param pars The parameters of the model.
#' @param parsHier The parameters of the Hierarchical model.
#' @param dataHier The data for the Hierarchical model.
#' @param funcForm The antibody kinetics function.
#' @return A list with the biomarker name, the exposure name, whether the exposure is inferred, the parameters and the antibody kinetics function.
#' @export
addAbkineticsModelHier <- function(id, biomarker, exposureType, pars, parsHier, dataHier, funcForm) {
    if (is.null(id) || is.null(biomarker) || is.null(exposureType) || 
        is.null(pars) || is.null(parsHier) || is.null(dataHier) || is.null(funcForm)) {
        stop("One or more arguments of `addAbkineticsModelHier` are NULL.")
    }
    
    # Check that elements of parsHier are a subset of pars
    if (!all(parsHier %in% pars)) {
        stop("All elements of parsHier must be a subset of pars.")
    }
    
    # Check that dataHier is a numeric vector
    if (!is.numeric(dataHier) || !is.vector(dataHier)) {
        stop("dataHier must be a numeric vector.")
    }
    
    # Check that funcForm is a function and its formals are (titre, time, pars)
    if (!is.function(funcForm)) {
        stop("funcForm must be a function.")
    }

    M <- length(unique(dataHier))
    add_pars <- c()
    k <- 1
    for (j in 1:length(pars)) {
        add_pars <- c(add_pars, pars[j])
        if (pars[j] %in% parsHier) {
            for (i in 1:M) {
                add_pars <- c(add_pars, paste0("z_" , parsHier[k], "_", i))
            }
            add_pars <- c(add_pars, paste0("sigma_", parsHier[k]))
            k <- k + 1
        }
    }


  
    list(
        id = id,
        biomarker = biomarker,
        exposureType = exposureType,
        pars = add_pars,
        parsBase = pars,
        parsHier = parsHier,
        dataHier = dataHier,
        dataHierN = M,
        funcForm = funcForm,
        hierFlag = TRUE
    )
}


# testing hierarchical model
logit_inverse <- function(x) {
    exp(x) / (1 + exp(x))
}


#' @title addCopModel
#' @description This function adds a cop model to the model definition.
#' @param biomarker The name of the biomarker.
#' @param exposureType The name of the exposure type.
#' @param pars The parameters of the model.
#' @param funcForm The functional form of the observational model.
#' @param logLikelihood The log likelihood function.
#' @return A list with the biomarker name, the exposure name, the parameters and the log likelihood function.
#' @export
addCopModel <- function(biomarker, exposureType, pars, funcForm, logLikelihood) {
    if (is.null(biomarker) || is.null(exposureType) || is.null(pars) || is.null(funcForm) || is.null(logLikelihood) ) {
        stop("One or more arguments of `addObservationalModel` are NULL.")
    }
    list(
        biomarker = biomarker,
        exposureType = exposureType,
        pars = pars,
        funcForm = funcForm,
        logLikelihood = logLikelihood
    )
}


#depth <- function(this,thisdepth=0){
 #   return(thisdepth)
##  if(!is.list(this)){
#  }else{
#  }
#    return(max(unlist(lapply(this,depth,thisdepth=thisdepth+1))))    
#}

#' @title makeModel
#' @description This function creates a model list for observationalModel, copModel, and abKineticsModel.
#' @param ... The models to include in the definition.
#' @return A list with the models.
#' @export
makeModel <- function(...) {
    if (length(list(...)) == 1) {
        models <- list(...)
    } else {
        models <- list(...)
    }
    models
}

#console_update <- function(data_t, modelSeroJump) {
#    name_bio <- modelSeroJump$observationalModel$name
#    cat("There are ", length(name_bio), " measured biomarkers: ", name_bio, "\n")
#    name_exp <- modelSeroJump$abkineticsModel$model$exposureNames
#    cat("There are ", length(name_exp), " exposure types in the study period: ", name_exp, "\n")
#    cat("The fitted exposure type is ", modelSeroJump$abkineticsModel$model$exposureNameInf)
#}


inf_prior_base <- function(N, E, I, K) {
    N_adj <- N - K
    E_adj <- E - K 
    logPriorExpInf <- lfactorial(E_adj) + lfactorial(N_adj - E_adj) - lfactorial(N_adj ) 
    logPriorExpInf
}


#' @title createSeroJumpModel
#' @description This function creates a model for the serology jump model.
#' @param data_sero The serology data.
#' @param data_known The known exposure data.
#' @param biomarkers A vector of the biomarkers 
#' @param exposureTypes A vector of the exposure types
#' @param exposureFitted The exposure type that is fitted
#' @param observationalModel The observational model
#' @param abkineticsModel The antibody kinetics model
#' @param exposurePriorTime The prior time of the exposure
#' @param exposurePriorTimeType The type of the exposure prior
#' @param exposurePriorPop The population prior of the exposure
#' @param known_exp_bool Boolean describing if the exposures are known in the model
#' @return A list with the data and the model.
#' @export
createSeroJumpModel <- function(
    data_sero, 
    data_known,
    biomarkers,
    exposureTypes,
    exposureFitted,
    observationalModel,
    abkineticsModel,
    exposurePriorTime = NULL,
    exposurePriorTimeType = NULL,
    exposurePriorPop = NULL,
    known_exp_bool = NULL, 
    seed = -1) {

    cat("OUTLINE OF INPUTTED MODEL\n")
    check_inputs(data_sero, data_known, biomarkers, exposureTypes, exposureFitted, observationalModel, abkineticsModel, exposurePriorTime, exposurePriorTimeType)
    check_priors(observationalModel, abkineticsModel)

    modelSeroJump <- list()

    # Generate data for model
    # 1. Extract priors and defined define distributions for model 
    priors <- bind_rows(
        observationalModel$prior,
        abkineticsModel$prior
    ) %>% filter(part_type == "prior")


    modelSeroJump$lowerParSupport_fitted <- priors$lb
    modelSeroJump$upperParSupport_fitted <- priors$ub
    modelSeroJump$namesOfParameters <- priors$par_name

    # Add in custom functions
    modelSeroJump$samplePriorDistributions <- function(datalist) {
        get_sample_non_centered(priors, seed)
    }

    modelSeroJump$evaluateLogPrior <- function(params, jump, datalist) {
        cal_lprior_non_centered(priors, params)
    }

    modelSeroJump$initialiseJump <- function(datalist) {
    }

    if(is.null(exposurePriorPop)) {
        modelSeroJump$evaluateLogPriorInfExp <- inf_prior_base
    } else {
        modelSeroJump$evaluateLogPriorInfExp <- exposurePriorPop
    } 


    # 2. Define the models for the likelihoods
    modelSeroJump$observationalModel <- observationalModel$model
    modelSeroJump$abkineticsModel <- abkineticsModel$model

    if (is.null(modelSeroJump$copModel)) {
        modelSeroJump$copModel <- list()
    }

    known_exp_bool <- NULL

    data_t <- generate_data(data_sero, biomarkers, known_exp_bool) # look into this later
    data_t$par_names <- priors[, 1]
    data_t$priorPredFlag <- FALSE
    T <- data_t$T

    # Add known infections to the model
    modelSeroJump$infoModel$biomarkers <- biomarkers
    modelSeroJump$infoModel$exposureFitted <- exposureFitted
    modelSeroJump$exposureTypes <- exposureTypes
    modelSeroJump$infoModel$exposureInfo <- list()
    modelSeroJump$infoModel$logitBoundaries <- abkineticsModel$prior %>% filter(part_type == "logit_boundary") %>% .[1:3]


    exp_prior <- convertExpPriorEmpirical(exposurePriorTime, T, exposurePriorTimeType = exposurePriorTimeType)
    exp_prior_list <- lapply(1:data_t$N, function(x) exp_prior$prob)

    # Extrac the information from the known exposures to the seorjump model
    know_inf <- list()
    for (i in 1:length(exposureTypes)) {
        # for each exposure type, extract the known exposure times
        modelSeroJump$infoModel$exposureInfo[[i]] <- list()
        exposureType <- exposureTypes[i]
        # The first exposure is the none
        if (i == 1 & (exposureType != exposureFitted)) {
            know_inf <- data_t$initialTitreTime
        } else {
            if (is.null(data_known)) {
                know_inf <- rep(-1, data_t$N)
            } else{
                initialTitreTime <- data_t$initialTitreTime
                endTitreTime <- data_t$endTitreTime

                # Trimmed info on known exposure
                know_inf_df <- check_oos_sample_df(data_known, exposureType, initialTitreTime, endTitreTime) 
                know_inf <- know_inf_df %>% 
                    mutate(time = case_when(time >= initialTitreTime & time <= endTitreTime ~ time, TRUE ~ -1)) %>%
                    pull(time)

                exp_prior_list <- update_exp_prior(exp_prior_list, know_inf_df, data_t$T, data_t$N)
            }
        }

        modelSeroJump$infoModel$exposureInfo[[i]]$exposureType <- exposureType
        modelSeroJump$infoModel$exposureInfo[[i]]$known_inf <- know_inf

        if (exposureType == exposureFitted) {
            data_t$knownInfsTimeVec = know_inf
            data_t$knownInfsVec = as.numeric(know_inf > -1)
            data_t$knownInfsN = length(know_inf[know_inf > -1])
        }
    }
    knownNoneInfsVec <- vector(mode = "numeric", length = data_t$N)
    for (i in 1:data_t$N) {
        if(sum(exp_prior_list[[i]]) == 0 && data_t$knownInfsVec[i] == -1 ) {
            knownNoneInfsVec <- 0
        } else {
            knownNoneInfsVec[i] <- 1
        }
    }    
    data_t$knownNoneInfsVec <- knownNoneInfsVec

    # then add information on 
    if (!is.null(known_exp_bool)) {
        data_t$knownExpVec <- data_t$knownInfsTimeVec
    }

    # Code to check form of exp_prior
    modelSeroJump$exposureFunctionSample <- function(i) {
        sample(1:length(exp_prior_list[[i]]), 1, prob = exp_prior_list[[i]]) 
    }

    if (!is.null(data_known)) {
        check_exposures_times(data_sero, data_known, exposureTypes, exposureFitted, exposurePriorTime, exposurePriorTimeType)
    }
    
    data_t$exp_list <- exp_prior_list
    data_t$exp_prior <- exp_prior

    modelSeroJump$exposureFunctionDensity <- function(jump_i, i) {
        if ((exp_prior_list[[i]][max(round(jump_i, 0), 1)] %>% log) < -100) {

        }
        exp_prior_list[[i]][max(round(jump_i, 0), 1)] %>% log
    }

    # Add help with exposure prior
    data_t$raw_sero <- data_sero
    data_t$raw_exp <- data_known
   # data_t$logitBoundaries <- modelSeroJump$logitBoundaries

    list(
        data = data_t,
        model = modelSeroJump
    )
}


#' @title run the RJMCMC algorithm
#' @description This function runs the RJMCMC algorithm given a defined seroModel, settings and filepaths for output
#' @param seroModel The seroModel previously defined.
#' @param settings Settings used for the calibration
#' @param priorPred Boolean option on whether to run the prior predictive model
#' @param save_info Filepath of where the outputs are saved
#' @param seed Seed for the random number generator 
#' @return A list with the posterior samples, the model and the data.
#' @export
runSeroJump <- function(seroModel, settings, priorPred = FALSE, save_info = NULL, seed = -1) {
   

    settings <- settings
    settings$numberFittedPar <- seroModel$model$namesOfParameters %>% length
    settings$lowerParBounds <- seroModel$model$lowerParSupport_fitted
    settings$upperParBounds <- seroModel$model$upperParSupport_fitted
    settings$lengthJumpVec <- seroModel$data$N
    settings <- check_settings_sero(settings)

    if (!is.null(save_info)) {
        dir.create(here::here("outputs", "fits"), recursive = TRUE, showWarnings = FALSE)
        check_save_info(save_info)
    }
   # dir.create(here::here("outputs", "fits", filename, "figs", modelname), recursive = TRUE, showWarnings = FALSE)

    if (priorPred) {
        seroModel$data$priorPredFlag <- TRUE
    }
    
    if(is.null(settings$runParallel)) {
        settings$runParallel <- TRUE
         cat("`runParallel` not specified in settings.  Default value TRUE. \n")
    }

    if(settings$runParallel) {
        out_pp_full <- mclapply(list(seroModel), 
        function(i) { 
            rjmc_sero_func(model = i$model, data = i$data, settings = settings, seed = seed)
        },
        mc.cores = settings$numberCores
        )
    } else {
        out_pp_full <- lapply(list(seroModel), 
            function(i) { 
                rjmc_sero_func(model = i$model, data = i$data, settings = settings, seed = seed)
            }
        )
    }

    model_fit <- list(post = out_pp_full[[1]],  model = seroModel$model, data_t = seroModel$data, settings = settings)
    model_post <- postprocess_fit(model_fit)

    model_summary <- list(fit = model_fit, post = model_post)

    if (!is.null(save_info)) {
        dir.create(here::here("outputs", "fits", save_info$file_name, save_info$model_name), recursive = TRUE, showWarnings = FALSE)
        if (seroModel$data$priorPredFlag) {
            saveRDS(model_summary, here::here("outputs", "fits", save_info$file_name, save_info$model_name, paste0("model_summary_pp.RDS")))
        } else {
            saveRDS(model_summary, here::here("outputs", "fits", save_info$file_name, save_info$model_name, paste0("model_summary.RDS")))
        }
    }
    model_summary
}

check_save_info <- function(save_info) {
    cat("--CHECKING SAVE_INFO DATA--\n")
    if (is.null(save_info[["file_name"]])) {
        stop("`file_name` must be specified in `save_info`.\n")
    } 
    if (is.null(save_info[["model_name"]])) {
        stop("`model_name` must be specified in `model_name`.\n")
    } 
    cat("Saving model outputs to: ", here::here("outputs", "fits", save_info[["file_name"]] , save_info[["model_name"]]), "\n")
}


#' @title Define the Log Likelihood function for the observational model
#' @name obsLogLikelihood
#'
#' @description This function is defined by the user and describes the observational likelihood of the data given the model.
#'
#' @param titre_val The titre value from the data
#' @param titre_est The model-estimate titre value.
#' @param pars The fitted parameters needed to calculate the log likelihood. These are defined in the \code{prior} entry of the observationalModel list.
#' @return A functoion that returns the log likelihood value for the COP model.
#' 
#' @details Add information here.
#' 
#' @examples
#' # Example usage:
#' obsLogLikelihood = function(titre_val, titre_est, pars) {
#'    if (titre_val <= log10(40)) {
#'        ll <- pcauchy(log10(40), titre_est, pars[1], log.p = TRUE)
#'    } else {
#'        ll <- dcauchy(titre_val, titre_est, pars[1], log = TRUE)
#'    }
#'    ll
#' }
#' 
#' 
#' @author Your Name
#' @keywords data format function
#' export
NULL



#' @title Define the Log Likelihood function for the COP model
#' @name copLogLikelihood
#'
#' @description This function is defined by the user and describes the COP likelihood of the data given the model.
#'
#' @param inf_status The infection status of an individual binary (0 or 1).
#' @param esttitreExp The model-estimated titre value at the time of exposure
#' @param pars The fitted parameters needed to calculate the log likelihood. These are defined in the \code{prior} entry of the copModel list.
#' @return A function that returns the log likelihood value for the COP model.
#' 
#' @details Add information here.
#' 
#' @examples
#' # Example usage:
#' copLogLikelihood <- function(inf_status, esttitreExp, params) {
#'    # COP parameters
#'    beta0 <- params[1]
#'    beta1 <- params[2]
#'    p <- 1.0 / (1.0 + exp(- (beta0 + beta1 * esttitreExp) ) )
#'    ll <- inf_status * log(p) + (1 - inf_status) * log(1 - p)
#'    ll
#'}
#' 
#' 
#' @author Your Name
#' @keywords data format function
#' export
NULL


#' @title Define the function form of th antibody kinetics
#' @name abkineticsFunction
#'
#' @description This function is defined by the user and describes antibody kinetics given the model.
#'
#' @param titre_est The model-estimated titre value at the observational time
#' @param inf_status How long since the last exposur event occured
#' @param pars The fitted parameters needed to calculate the estimate titre value. These are defined in the \code{prior} entry of the abkineticsModel list.
#' @return A function that returns estimate titre value 
#' 
#' @details Add information here.
#' 
#' @examples
#' # Example usage: we define two functional form examples, one describing waning until infection, and then the other describing kinetics post infection. 
#'noInfSerumKinetics <- function(titre_est, timeSince, pars) {
#'    titre_est <- titre_est - pars[1] * (timeSince)
#'    titre_est
#'}
#'
#'infSerumKinetics <- function(titre_est, timeSince, pars) {
#'    a <- pars[1] 
#'    b <- pars[2] 
#'    c <- pars[3]
#'    if (timeSince < 14) {
#'        titre_est <- titre_est + (log(exp(a) + exp(c)) * (timeSince) / 14) * 1;
#'    } else {
#'        titre_est <- titre_est + (log(exp(a) * exp(-b/10 * (timeSince - 14))) + exp(c)) * 1;
#'    }
#'    titre_est
#'}
#'

#' @author Your Name
#' @keywords data format function
#' export
NULL