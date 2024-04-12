#' @title addObservationalModel
#' @description This function adds an observational model to the model definition.
#' @param biomarker The name of the biomarker.
#' @param pars The parameters of the model.
#' @param logLikelihood The log likelihood function.
#' @return A list with the biomarker name, the parameters and the log likelihood function.
#' @export
addObservationalModel <- function(biomarker, pars, logLikelihood) {
    list(
        biomarker = biomarker,
        pars = pars,
        logLikelihood = logLikelihood
    )
}

#' @title addAbkineticsModel
#' @description This function adds an antibody kinetics model to the model definition.
#' @param biomarker The name of the biomarker.
#' @param exposureName The name of the exposure type.
#' @param inferred Whether the exposure is inferred with in the RJMCMC algorithm (can only do this for ONE exposure type).
#' @param pars The parameters of the model.
#' @param funcForm The antibody kinetics function.
#' @return A list with the biomarker name, the exposure name, whether the exposure is inferred, the parameters and the antibody kinetics function.
#' @export
addAbkineticsModel <- function(biomarker, exposureName, inferred, pars, funcForm) {
    list(
        biomarker = biomarker,
        exposureName = exposureName,
        inferred = inferred,
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
addCopModel <- function(biomarker, exposureName, pars, logLikelihood) {
    list(
        biomarker = biomarker,
        exposureName = exposureName,
        pars = pars,
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


#' @title createSeroJumpModel
#' @description This function creates a model for the serology jump model.
#' @param data_sero The serology data.
#' @param data_known The known exposure data.
#' @param modeldefinition The model definition.
#' @return A list with the data and the model.
#' @export
createSeroJumpModel <- function(data_sero, data_known, modeldefinition) {
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

    # Checks on this data.frame, i) no duplicated names, ii) all parameters have a priors defined, iii) lb < ub, iv) dist is a valid distribution, v) draw all priors

    modelSeroJump$lowerParSupport_fitted <- priors$lb
    modelSeroJump$upperParSupport_fitted <- priors$ub
    modelSeroJump$namesOfParameters <- priors$par_name

    modelSeroJump$samplePriorDistributions <- function(datalist) {
        get_sample_non_centered(priors)
    }

    modelSeroJump$evaluateLogPrior <- function(params, jump, datalist) {
        cal_lprior_non_centered(priors, params)
    }

    for (i in 1:length( abkineticsModel$model)) {
        if(abkineticsModel$model[[i]]$inferred) {
            inferredExpType = abkineticsModel$model[[i]]$exposureName
        }
    }

    # 2. Define the exposure prior 

    # 3. Define the models for the likelihoods
    modelSeroJump$observationalModel <- modeldefinition$observationalModel
    modelSeroJump$abkineticsModel <- modeldefinition$abkineticsModel
    modelSeroJump$copModel <- modeldefinition$copModel
    names_here <- modeldefinition$abkineticsModel$names
    names(modelSeroJump$abkineticsModel$model) <- names_here


    data_t <- generate_data_alt(data_sero, modeldefinition$biomarkers)
    data_t$par_names <- priors[, 1]

    know_inf <- list()
    for (i in 1:length(modeldefinition$abkineticsModel$model)) {
        if (i == 1) {
            know_inf[[names_here[i]]] <- data_t$initialTitreValue
        } else {
            if (is.null(data_known)) {
                know_inf[[names_here[i]]] <- rep(-1, data_t$N)
            } else{
                know_inf[[names_here[i]]] = data_known %>% filter(exposure_type == names_here[i]) %>%
                    complete(id = 1:data_t$N, fill = list(time = -1)) %>% 
                    mutate(start = data_t$initialTitreTime, end = data_t$endTitreTime) %>%
                    mutate(time = case_when(time >= start & time <= end ~ time, TRUE ~ -1)) %>% pull(time)
            }
        }
        modelSeroJump$abkineticsModel$model[[names_here[i]]]$known_inf <- know_inf[[names_here[i]]]
        if (modelSeroJump$abkineticsModel$model[[names_here[i]]]$inferred) {
            data_t$knownInfsTimeVec = know_inf[[names_here[i]]]
            data_t$knownInfsVec = as.numeric(know_inf[[names_here[i]]] > -1)
            data_t$knownInfsN = sum(know_inf[[names_here[i]]] > -1)
        }
    }

    modelSeroJump$initialiseJump <- function(datalist) {
        init_exposure <- data_t$knownInfsTimeVec
        init_exposure
    }
    

    modelSeroJump <- addExposurePrior(modelSeroJump, data_t, modeldefinition$exposurePrior, type = modeldefinition$exposurePriorType)
    modelSeroJump$abkineticsModel$model$exposurePeriods <- abkineticsModel$names
    modelSeroJump$abkineticsModel$model$exposureNames <- modeldefinition$exposureType
    modelSeroJump$abkineticsModel$model$exposureNameInf <- inferredExpType

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

    dir.create(here::here("outputs"), showWarnings = FALSE)
    dir.create(here::here("outputs", "fits", filename), showWarnings = FALSE)
    dir.create(here::here("outputs", "fits", filename, "figs"), showWarnings = FALSE)

    saveRDS(fitfull, here::here("outputs", "fits", filename, paste0("fit_", modelname, ".RDS")))
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