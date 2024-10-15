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
#' @param is The name of the biomarker.
#' @param biomarker The name of the biomarker.
#' @param exposureName The name of the exposure type.
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
        funcForm = funcForm
    )
}

#' @title addCopModel
#' @description This function adds a cop model to the model definition.
#' @param biomarker The name of the biomarker.
#' @param exposureName The name of the exposure type.
#' @param pars The parameters of the model.
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


depth <- function(this,thisdepth=0){
  if(!is.list(this)){
    return(thisdepth)
  }else{
    return(max(unlist(lapply(this,depth,thisdepth=thisdepth+1))))    
  }
}

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

console_update <- function(data_t, modelSeroJump) {
    name_bio <- modelSeroJump$observationalModel$name
    cat("There are ", length(name_bio), " measured biomarkers: ", name_bio, "\n")
    name_exp <- modelSeroJump$abkineticsModel$model$exposureNames
    cat("There are ", length(name_exp), " exposure types in the study period: ", name_exp, "\n")
    cat("The fitted exposure type is ", modelSeroJump$abkineticsModel$model$exposureNameInf)
}

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
#' @param modeldefinition The model definition.
#' @return A list with the data and the model.
#' @export
createSeroJumpModel <- function(data_sero, data_known, modeldefinition, known_exp_bool = NULL) {
    cat("OUTLINE OF INPUTTED MODEL\n")
    check_inputs(data_sero, data_known, modeldefinition)
    check_priors(modeldefinition)

    modelSeroJump <- list()

    # Generate data for model
    # 1. Extract priors and defined define distributions for model 
    priors <- bind_rows(
        modeldefinition$observationalModel$prior,
        modeldefinition$abkineticsModel$prior,
        modeldefinition$copModel$prior
    )

    modelSeroJump$lowerParSupport_fitted <- priors$lb
    modelSeroJump$upperParSupport_fitted <- priors$ub
    modelSeroJump$namesOfParameters <- priors$par_name

    # Add in custom functions
    modelSeroJump$samplePriorDistributions <- function(datalist) {
        get_sample_non_centered(priors)
    }

    modelSeroJump$evaluateLogPrior <- function(params, jump, datalist) {
        cal_lprior_non_centered(priors, params)
    }

    modelSeroJump$initialiseJump <- function(datalist) {
    }

    if(is.null(modeldefinition$expInfPrior)) {
        modelSeroJump$evaluateLogPriorInfExp <- inf_prior_base
    } else {
        modelSeroJump$evaluateLogPriorInfExp <- modeldefinition$expInfPrior
    } 


    # 2. Define the models for the likelihoods
    modelSeroJump$observationalModel <- modeldefinition$observationalModel$model
    modelSeroJump$abkineticsModel <- modeldefinition$abkineticsModel$model
    modelSeroJump$copModel <- modeldefinition$copModel$model

    if (is.null(modelSeroJump$copModel)) {
        modelSeroJump$copModel <- list()
    }

   # id_ab <- modeldefinition$abkineticsModel$model %>% map(~.x$id) %>% unlist
   # names(modelSeroJump$abkineticsModel) <- id_ab
    known_exp_bool <- NULL
    data_t <- generate_data_alt(data_sero, modeldefinition$biomarkers, known_exp_bool)
    data_t$par_names <- priors[, 1]
    data_t$priorPredFlag <- FALSE

    # Add known infections to the model
    modelSeroJump$infoModel$biomarkers <- modeldefinition$biomarkers
    modelSeroJump$infoModel$exposureFitted <- modeldefinition$exposureFitted

    modelSeroJump$exposureTypes <- modeldefinition$exposureTypes

    modelSeroJump$infoModel$exposureInfo <- list()


        know_inf <- list()
        for (i in 1:length(modeldefinition$exposureTypes)) {
            modelSeroJump$infoModel$exposureInfo[[i]] <- list()
            exposureType <- modeldefinition$exposureTypes[i]

            if (i == 1 & (exposureType != modeldefinition$exposureFitted)) {
                know_inf <- data_t$initialTitreTime
            } else {
                if (is.null(data_known)) {
                    know_inf <- rep(-1, data_t$N)
                } else{
                    know_inf <- data_known %>% filter(exposure_type == exposureType) %>%
                        complete(id = 1:data_t$N, fill = list(time = -1)) %>% 
                        mutate(start = data_t$initialTitreTime, end = data_t$endTitreTime) %>%
                        mutate(time = case_when(time >= start & time <= end ~ time, TRUE ~ -1)) %>% pull(time)
                }
            }
            modelSeroJump$infoModel$exposureInfo[[i]]$exposureType <- exposureType
            modelSeroJump$infoModel$exposureInfo[[i]]$known_inf <- know_inf

            if (exposureType == modeldefinition$exposureFitted) {
                data_t$knownInfsTimeVec = know_inf
                data_t$knownInfsVec = as.numeric(know_inf > -1)
                data_t$knownInfsN = length(know_inf[know_inf > -1])
            }
        }

    if (!is.null(known_exp_bool)) {
        data_t$knownExpVec <- data_t$knownInfsTimeVec
    }

    data_t <- calculateIndExposure(modelSeroJump, data_t, modeldefinition$exposurePrior, type = modeldefinition$exposurePriorType)

    # Code to check form of exp_prior
    modelSeroJump$exposureFunctionSample <- function(i) {
        sample(1:length(data_t$exp_list[[i]]), 1, prob = data_t$exp_list[[i]]) 
    }

    modelSeroJump$exposureFunctionDensity <- function(jump_i, i) {
        if ((data_t$exp_list[[i]][max(round(jump_i, 0), 1)] %>% log) < -100) {

        }
        data_t$exp_list[[i]][max(round(jump_i, 0), 1)] %>% log
    }

    # Add help with exposure prior
    data_t$raw_sero <- data_sero
    data_t$raw_exp <- data_known

    list(
        data = data_t,
        model = modelSeroJump
    )
}


#' @title run the RJMCMC algorithm
#' @name runSeroJump
#' @description This function runs the RJMCMC algorithm given a defined seroModel, settings and filepaths for output
#' @param seroModel The seroModel previously defined.
#' @param settings Settings used for the calibration
#' @param filename Filepath of where the outputs are saved (outputs/fits/filename)
#' @param modelname Name of the model outputs (in outputs/fits/filename)
#' @return A list with the posterior samples, the model and the data.
#' @export
runSeroJump <- function(seroModel, settings, priorPred = FALSE, save_info = NULL) {
   

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
            rjmc_sero_func(model = i$model, data = i$data, settings = settings)
        },
        mc.cores = 8
        )
    } else {
        out_pp_full <- lapply(list(seroModel), 
            function(i) { 
                rjmc_sero_func(model = i$model, data = i$data, settings = settings)
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
#' @param pars The fitted parameters needed to calculate the log likelihood. These are defined in the \code{prior} entry of the \ref{observationalModel} list.
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
#' @param pars The fitted parameters needed to calculate the log likelihood. These are defined in the \code{prior} entry of the \ref{copModel} list.
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
#' @param pars The fitted parameters needed to calculate the estimate titre value. These are defined in the \code{prior} entry of the \ref{abkineticsModel} list.
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