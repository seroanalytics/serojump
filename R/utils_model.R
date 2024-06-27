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

prior_predictive_ab <- function(seroModel) {
    abModels <- seroModel$model$abkineticsModel 
    M <- length(abModels)
    seroModel$model$samplePriorDistributions() 

    full_pp_ab <- map_df(
        1:4,
        function(i) {
            abModel <- abModels[[i]]
            pars <- abModel$pars
            funcForm <- abModel$funcForm
            id <- abModel$id

            map_df(1:1000,
                function(j) {
                    samples <- as.numeric(seroModel$model$samplePriorDistributions()[pars] )
                    prior_samples <- map_dbl(1:250, ~funcForm(0, .x, samples) )
                    data.frame(
                        t = 1:250,
                        titre = prior_samples,
                        sample = j,
                        name = id
                    )
                }
            )
        }
    )
    full_pp_ab_sum <- full_pp_ab %>% group_by(t, name) %>% mean_qi(titre)

    full_pp_ab_sum %>%
        ggplot() + 
            geom_ribbon(aes(x = t, ymin = .lower, ymax = .upper), alpha = 0.3) + 
            geom_line(aes(x = t, titre), size = 3) + 
            facet_wrap(vars(name)) + 
            theme_minimal()
}



check_inputs <- function(data_sero, data_known, modeldefinition) {
    # CHECK inputs of modeldefinition are present
    if(is.null(modeldefinition$biomarkers)) {
        stop("Please define the `biomarkers` variable in `modeldefinition`")
    }
    if(is.null(modeldefinition$exposureTypes)) {
        stop("Please define the `exposureTypes` variable in `modeldefinition`")
    }
    if(is.null(modeldefinition$exposureFitted)) {
        stop("Please define the `exposureFitted` variable in `modeldefinition`")
    }
    if(is.null(modeldefinition$observationalModel)) {
        stop("Please define the `observationalModel` structure in `modeldefinition`")
    }
    if(is.null(modeldefinition$abkineticsModel)) {
        stop("Please define the `abkineticsModel` structure in `modeldefinition`")
    }
   # if(is.null(modeldefinition$copModel)) {
  #      stop("Please define the `copModel` structure in `modeldefinition`")
  #  }
    if(is.null(modeldefinition$exposurePrior)) {
        stop("Please define the `exposurePrior` data in `modeldefinition`")
    } 
    if(is.null(modeldefinition$exposurePriorType)) {
        stop("Please define the `exposurePriorType` flag (func or empirical) in `modeldefinition`")
    } 
 
    # CHECK BIOMARKERS ARE WELL DEFINED
    # Check columns of data_sero match model definition 
    data_sero_name <- data_sero %>% names
    biomarkers_md <- modeldefinition$biomarkers
    biomarkers_cop <- modeldefinition$copModel$model %>% map(~.x$biomarker) %>% unlist
    biomarkers_obs <- modeldefinition$observationalModel$model %>% map(~.x$biomarker) %>% unlist
    biomarkers_abkin <- modeldefinition$abkineticsModel$model %>% map(~.x$biomarker) %>% unlist %>% unique

    for(b in biomarkers_md) {
        if(!b %in% data_sero_name) {
            stop("Biomarker, ", b, ", in `modeldefinition$biomarkers` is not a column of serological data; `",
                paste(data_sero_name, collapse = ", "), "`")
        }
    }
    if(!identical(biomarkers_md, biomarkers_cop) ) {
        stop("Biomarkers in copModel (", paste(biomarkers_cop, collapse = ", "), ") do not match biomarkers in `modeldefinition$biomarkers` (", paste(biomarkers_md, collapse = ", "), ")")
    }
     if(!identical(biomarkers_md, biomarkers_obs) ) {
        stop("Biomarkers in observationalModel (", paste(biomarkers_obs, collapse = ", "), ") do not match biomarkers in `modeldefinition$biomarkers` (", paste(biomarkers_md, collapse = ", "), ")")
    }
    if(!identical(biomarkers_md, biomarkers_abkin) ) {
        stop("Biomarkers in abkineticsModel (", paste(biomarkers_abkin, collapse = ", "), ") do not match biomarkers in `modeldefinition$biomarkers` (", paste(biomarkers_md, collapse = ", "), ")")
    }  

    # CHECK EXPSURETYPES ARE WELL DEFINED
    exposure_type_names <- data_known$exposure_type %>% unique
    exposures_md <- modeldefinition$exposureTypes
    exposures_obs <- modeldefinition$abkineticsModel$model %>% map(~.x$exposureType) %>% unlist %>% unique
    if (!is.null(data_known)) {
        exposure_type_names <- data_known$exposure_type %>% unique
        for(e in exposure_type_names) {
            if(!e %in% exposures_md) {
                stop("Exposure type, ", e, ", in known exposure data.frame column 'exposure_type' (", paste(exposure_type_names, collapse = ", "),
                    "is not defined in `modeldefinition$exposureTypes` (", paste(exposures_md, collapse = ", "), ")")
            }
        }
    }
    if(!identical(exposures_md, exposures_obs) ) {
        stop("Exposure types in abkineticsModel (", paste(exposures_obs, collapse = ", "), ") do not match exposure types in `modeldefinition$exposureTypes` (", paste(exposures_md, collapse = ", "), ")")
    }
    if(is.null(modeldefinition$exposureFitted)) {
        stop("`modeldefinition$exposureFitted` is NULL, please define a biomarker to fit.")
    }
    exposure_fitted <- modeldefinition$exposureFitted
    if(!exposure_fitted %in% exposures_md) {
        stop("The fitted exposure type, ", exposure_fitted, ", is not defined in, `modeldefinition$exposureTypes`: ", paste(exposures_md, collapse = ", "))
    }

    names_cop <- modeldefinition$copModel$model %>% map(~.x$name) %>% unlist
    names_obs <- modeldefinition$observationalModel$model %>% map(~.x$name) %>% unlist
    names_abkin <- modeldefinition$abkineticsModel$model %>% map(~.x$name) %>% unlist %>% unique
    # Read out into console
    cat("There are ", length(biomarkers_md), " measured biomarkers: ", paste(biomarkers_md, collapse = ", "), "\n")
    cat("There are ", length(exposures_md), " exposure types in the study period: ", paste(exposures_md, collapse = ", "), "\n")
    cat("The fitted exposure type is ", modeldefinition$exposureFitted, "\n")

}

check_priors <- function(modeldefinition) {
    priors <- bind_rows(
        modeldefinition$observationalModel$prior,
        modeldefinition$abkineticsModel$prior,
        modeldefinition$copModel$prior
    )
    if(any(duplicated(priors$par_name))) {
        stop("Priors: ", paste0(priors$par_name[duplicated(priors$par_name)], collapse = ", "), " are duplicated, please assign original names to each prior")
    }
    if(any(priors$lb >= priors$ub)) {
        stop("Priors: ", paste(priors$par_name[any(priors$lb < priors$ub)], collapse = ", "), " have their lower bound greater than or equal to upper bound, please change.")
    }
    for (i in 1:nrow(priors)) {
        func <- paste("r", priors$dist[[i]], sep = "") 
        if(!exists( func ) ){
            stop("Prior function `",  priors$dist[[i]], "` for `", priors$par_name[[i]], "` is not defined in R environment.")
        }
    }

    cat("PRIOR DISTRIBUTIONS", "\n")
    cat("Prior parameters of observationalModel are: ", paste(modeldefinition$observationalModel$prior$par_name, collapse = ", "), "\n")
    cat("Prior parameters of abkineticsModel are: ", paste(modeldefinition$abkineticsModel$prior$par_name, collapse = ", "), "\n")
    cat("Prior parameters of copModel are: ", paste(modeldefinition$copModel$prior$par_name, collapse = ", "), "\n")

}

#' @title createSeroJumpModel
#' @description This function creates a model for the serology jump model.
#' @param data_sero The serology data.
#' @param data_known The known exposure data.
#' @param modeldefinition The model definition.
#' @return A list with the data and the model.
#' @export
createSeroJumpModel <- function(data_sero, data_known, modeldefinition, known_exp = NULL) {
    cat("OUTLINE OF INPUTTED MODEL\n")
   # check_inputs(data_sero, data_known, modeldefinition)
    check_priors(modeldefinition)

    #data_sero <- gambia_pvnt_w2
#    data_known <- known_exposure

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

    modelSeroJump$evaluateLogPriorInfExp <- modeldefinition$expInfPrior



    # 2. Define the models for the likelihoods
    modelSeroJump$observationalModel <- modeldefinition$observationalModel$model
    modelSeroJump$abkineticsModel <- modeldefinition$abkineticsModel$model
    modelSeroJump$copModel <- modeldefinition$copModel$model
   # id_ab <- modeldefinition$abkineticsModel$model %>% map(~.x$id) %>% unlist
   # names(modelSeroJump$abkineticsModel) <- id_ab

    data_t <- generate_data_alt(data_sero, modeldefinition$biomarkers, known_exp)
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

    modelSeroJump <- addExposurePrior(modelSeroJump, data_t, modeldefinition$exposurePrior, type = modeldefinition$exposurePriorType)

    
    # Add help with exposure prior

    list(
        data = data_t,
        model = modelSeroJump
    )
}



#' @title run the RJMCMC algorithm
#' @name runRJMCMC
#' @description This function runs the RJMCMC algorithm given a defined seroModel, settings and filepaths for output
#' @param seroModel The seroModel previously defined.
#' @param settings Settings used for the calibration
#' @param filename Filepath of where the outputs are saved (outputs/fits/filename)
#' @param modelname Name of the model outputs (in outputs/fits/filename)
#' @return A list with the posterior samples, the model and the data.
#' @export
runRJMCMC <- function(seroModel, settings, filename, modelname) {
    settings <- settings
    settings$numberFittedPar <- seroModel$model$namesOfParameters %>% length
    settings$lowerParBounds <- seroModel$model$lowerParSupport_fitted
    settings$upperParBounds <- seroModel$model$upperParSupport_fitted
    settings$lengthJumpVec <- seroModel$data$N
    
    post <- rjmc_full_func(model = seroModel$model, data = seroModel$data, settings = settings)
    fitfull <- list(post = post,  model = seroModel$model, data_t = seroModel$data)

    dir.create(here::here("outputs", "fits", filename, modelname, "figs"), recursive = TRUE, showWarnings = FALSE)
    saveRDS(fitfull, here::here("outputs", "fits", filename, modelname, paste0("fit_", modelname, ".RDS")))
}

#' @title run the RJMCMC algorithm
#' @name runRJMCMC
#' @description This function runs the RJMCMC algorithm given a defined seroModel, settings and filepaths for output
#' @param seroModel The seroModel previously defined.
#' @param settings Settings used for the calibration
#' @param filename Filepath of where the outputs are saved (outputs/fits/filename)
#' @param modelname Name of the model outputs (in outputs/fits/filename)
#' @return A list with the posterior samples, the model and the data.
#' @export
runInfRJMCMC <- function(seroModel, settings, filename, modelname, priorPred = TRUE) {

    settings <- settings
    settings$numberFittedPar <- seroModel$model$namesOfParameters %>% length
    settings$lowerParBounds <- seroModel$model$lowerParSupport_fitted
    settings$upperParBounds <- seroModel$model$upperParSupport_fitted
    settings$lengthJumpVec <- seroModel$data$N

    dir.create(here::here("outputs", "fits", filename, "figs", modelname), recursive = TRUE, showWarnings = FALSE)

    if (!priorPred) {
        out_pp <- rjmc_sero_func(model = seroModel$model, data = seroModel$data, settings = settings)
        fitfull <- list(post = out_pp,  model = seroModel$model, data_t = seroModel$data)
        saveRDS(fitfull, here::here("outputs", "fits", filename, paste0("fit_", modelname, ".RDS")))
    } else {
        seroModel_pp <- seroModel
        seroModel_pp$data$priorPredFlag <- TRUE

        if(settings$runParallel) {
            out_pp_full <- mclapply(list(seroModel_pp, seroModel), 
            function(i) { 
                rjmc_sero_func(model = i$model, data = i$data, settings = settings)
            },
            mc.cores = 8
            )
        } else {
            out_pp_full <- lapply(list(seroModel_pp, seroModel), 
            function(i) { 
                rjmc_sero_func(model = i$model, data = i$data, settings = settings)
            }
            )
        }

        fitfull_pp <- list(post = out_pp_full[[1]],  model = seroModel_pp$model, data_t = seroModel_pp$data)
        fitfull <- list(post = out_pp_full[[2]],  model = seroModel$model, data_t = seroModel$data)

        dir.create(here::here("outputs", "fits", filename, "figs", modelname), recursive = TRUE, showWarnings = FALSE)
        saveRDS(fitfull_pp, here::here("outputs", "fits", filename, paste0("fit_prior_", modelname, ".RDS")))
        saveRDS(fitfull, here::here("outputs", "fits", filename, paste0("fit_", modelname, ".RDS")))
    }
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