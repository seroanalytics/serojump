#' @useDynLib serojump
#' @importFrom Rcpp sourceCpp
#' @import coda
#' @import parallel
#' @import tidyr
#' @import dplyr
#' @import foreach
#' @import doParallel
#' @importFrom magrittr %>% %<>%
NULL

#' Create a model
#'
#' @param model The model to run ptmc.
#' @param data Data used in the calibration process
#' @param settings settings
#' @param par NULL
#' @return Returns a list with the fist element being the mcmc samples formatting for analysis and plottig with the CODA package. The second is the log posterior value at each time step.
#'
#' @export
rjmc_sero_func <- function(model, data, settings, par = NULL, seed = -1) {
  settings <- check_settings_sero(settings)
  if (length(par) == 0) {
    par <- rep(list(list(type = "None")), settings[["numberChainRuns"]])
    output <- get_output_sero(model, data, settings, FALSE, par, seed)
  } else {
    output <- get_output_sero(model, data, settings, TRUE, par, seed)
  }
  output
}

get_output_sero <- function(model, data_list, settings, update_ind, par, seed = -1) {

  outPTpost <- vector(mode = "list", length = settings[["numberChainRuns"]])
  outPTjump <- vector(mode = "list", length = settings[["numberChainRuns"]])
  outPTtitreexp <- vector(mode = "list", length = settings[["numberChainRuns"]])
  outPTobstitre <- vector(mode = "list", length = settings[["numberChainRuns"]])
  outPTobsloglik <- vector(mode = "list", length = settings[["numberChainRuns"]])
  outPTlp <- vector(mode = "list", length = settings[["numberChainRuns"]])
  outPTacc <- vector(mode = "list", length = settings[["numberChainRuns"]])
  outPTpar <- vector(mode = "list", length = settings[["numberChainRuns"]])
  out_raw <- vector(mode = "list", length = settings[["numberChainRuns"]])


  # Run the chains in parallel
  if (settings[["runParallel"]]) {
    out_raw <- mclapply(1:settings[["numberChainRuns"]], 
      function(i) {
        run_rjmc_sero(model, data_list, settings, update_ind, par[[i]], i, seed)
      },
      mc.cores = settings[["numberCores"]]
    )
  } else {
    for (i in 1:settings[["numberChainRuns"]]) {
      cat(i)
      out_raw[[i]] <- run_rjmc_sero(model, data_list, settings, update_ind, par[[i]], i)
    }
  }

  for(i in 1:settings[["numberChainRuns"]]) {
    out_post <- out_raw[[i]][["output"]][, 1:settings$numberFittedPar]
    outPTpar[[i]] <- out_raw[[i]][["PTMCpar"]]
    if (settings$numberFittedPar > 1){
        colnames(out_post) <- model[["namesOfParameters"]]
    }
    outPTpost[[i]] <- mcmc(out_post)
    outPTjump[[i]] <- out_raw[[i]][["jump"]]
    outPTtitreexp[[i]] <- out_raw[[i]][["titreexp"]]
    outPTobstitre[[i]] <- out_raw[[i]][["obstitre"]]
    outPTobsloglik[[i]] <- out_raw[[i]][["obsloglik"]]
    outPTlp[[i]] <- out_raw[[i]][["output"]][, settings$numberFittedPar + 1]
    outPTacc[[i]] <- out_raw[[i]][["output"]][, settings$numberFittedPar + 2]
  }

  outlpv <- data.frame(matrix(unlist(outPTlp), nrow = length(outPTlp[[1]])))
  colnames(outlpv) <- c(1:settings[["numberChainRuns"]])
  outlpv <- outlpv %>% gather(colnames(outlpv), key="chain_no",value="lpost")
  outlpv$sample_no <- rep(1:length(outPTlp[[1]]), settings[["numberChainRuns"]])

  outlaccv <- data.frame(matrix(unlist(outPTacc), nrow=length(outPTacc[[1]])))
  colnames(outlaccv) <- c(1:settings[["numberChainRuns"]])
  outlaccv <- outlaccv %>% gather(colnames(outlaccv), key="chain_no", value="acceptance rate")
  outlaccv$sample_no <- rep(1:length(outPTacc[[1]]), settings[["numberChainRuns"]])

  output <- list(
    mcmc = as.mcmc.list(outPTpost),
    jump = outPTjump,
    titreexp = outPTtitreexp,
    obstitre = outPTobstitre,
    obsloglik = outPTobsloglik,
    lpost = outlpv,
    acc = outlaccv,
    outPTpar = outPTpar
  )
  output
}

check_settings_sero <- function(settings) {
    cat("--CHECKING SETTINGS LIST--\n")

  if (is.null(settings[["numberChainRuns"]])) {
    settings[["numberChainRuns"]] <- 4
    cat("`numberChainRuns` not specified in settings. Default value 4. \n")
  }

  if (is.null(settings[["numberCores"]])) {
    settings[["numberCores"]] <- settings[["numberChainRuns"]]
    cat("`numberCores` not specified in settings. Default value equal to `numberChainRuns`. \n")
  }

  if (is.null(settings[["iterations"]])) {
    settings[["iterations"]] <- 20000
    cat("`iterations` not specified in settings. Default value 20,000. \n")
  }  
  if (is.null(settings[["burninPosterior"]])) {
    settings[["burninPosterior"]] <- 10000
    cat("`numberChainRuns` not specified in settings. Default value 10,000. \n")
  }  
  if (is.null(settings[["thin"]])) {
    settings[["thin"]] <- 100
    cat("`thin` not specified in settings. Default value 100. \n")
  }
  if (is.null(settings[["consoleUpdates"]])) {
    settings[["consoleUpdates"]] <- 100
    cat("`consoleUpdates` not specified in settings. Default value 100. \n")
  }
  if (is.null(settings[["numberFittedPar"]])) {
    stop("`numberFittedPar` not specified in settings. MUST be specified. \n")
  }
  if (is.null(settings[["onAdaptiveCov"]])) {
        settings[["onAdaptiveCov"]] <- TRUE
    cat("`onAdaptiveCov` not specified in settings. Default value TRUE. \n")
  }
  if (is.null(settings[["updatesAdaptiveCov"]])) {
        settings[["updatesAdaptiveCov"]] <- 100
    cat("`updatesAdaptiveCov` not specified in settings. Default value 100. \n")
  }
  if (is.null(settings[["burninAdaptiveCov"]])) {
        settings[["burninAdaptiveCov"]] <- 1000
    cat("`burninAdaptiveCov` not specified in settings. Default value 1000. \n")
  }
  if (is.null(settings[["onAdaptiveTemp"]])) {
        settings[["onAdaptiveTemp"]] <- TRUE
    cat("`onAdaptiveTemp` not specified in settings.  Default value TRUE. \n")
  }
  if (is.null(settings[["updatesAdaptiveTemp"]])) {
        settings[["updatesAdaptiveTemp"]] <- 10
    cat("`updatesAdaptiveTemp` not specified in settings.  Default value 10. \n")
  }
  if (is.null(settings[["onDebug"]])) {
        settings[["onDebug"]] <- FALSE
  }
  if (is.null(settings[["lowerParBounds"]])) {
    stop("`lowerParBounds` not specified in settings. MUST be specified. \n")
  }
  if (is.null(settings[["upperParBounds"]])) {
    stop("`upperParBounds` not specified in settings. MUST be specified. \n")
  }
  if (is.null(settings[["covarInitVal"]])) {
        settings[["covarInitVal"]] <- 1e-2
    cat("`covarInitVal` not specified in settings.  Default value 1e-10. \n")
  }
  if (is.null(settings[["covarInitValAdapt"]])) {
        settings[["covarInitValAdapt"]] <- 1e-2
    cat("`covarInitValAdapt` not specified in settings.  Default value 1e-10. \n")
  }
  if (is.null(settings[["covarMaxVal"]])) {
        settings[["covarMaxVal"]] <- 1
    cat("`covarMaxVal` not specified in settings. Default value 1. \n")
  }
  if (is.null(settings[["noGibbsSteps"]])) {
    settings[["noGibbsSteps"]] <- 1
    cat("`noGibbsSteps` not specified in settings. Default value 1. \n")
  }  
  if (is.null(settings[["runParallel"]])) {
        settings[["runParallel"]] <- TRUE
    cat("`runParallel` not specified in settings. Default value TRUE. \n")
  }

  settings
}