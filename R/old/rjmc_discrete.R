#' @useDynLib ptmc
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
#' @return Returns a list with the fist element being the mcmc samples formatting for analysis and plottig with the CODA package. The second is the log posterior value at each time step.
#'
#' @export
rjmc_alt_func <- function(model, data, settings, par = NULL) {
    settings <- check_settings_discete(settings)

  if (length(par) == 0) {
    par <- rep(list(list(type = "None")), settings[["numberChainRuns"]])
    output <- get_alt_output(model, data, settings, FALSE, par)
  } else {
    output <- get_alt_output(model, data, settings, TRUE, par)
  }
  output
}

get_alt_output <- function(model, data_list, settings, update_ind, par) {

  outPTpost <- vector(mode = "list", length = settings[["numberChainRuns"]])
  outPTmissed <- vector(mode = "list", length = settings[["numberChainRuns"]])
  outPTlp <- vector(mode = "list", length = settings[["numberChainRuns"]])
  outPTtemp <- vector(mode = "list", length = settings[["numberChainRuns"]])
  outPTacc <- vector(mode = "list", length = settings[["numberChainRuns"]])
  outPTpar <- vector(mode = "list", length = settings[["numberChainRuns"]])
  out_raw <- list()

  # Run the chains in parallel
  if (settings[["runParallel"]]) {
    out_raw <- mclapply(1:settings[["numberChainRuns"]], 
      function(i) {
        run_ptmc_discrete(model, data_list, settings, update_ind, par[[i]], i)
      },
      mc.cores = settings[["numberCores"]]
    )
  } else {
    for (i in 1:settings[["numberChainRuns"]]) {
      out_raw[[i]] <- run_ptmc_discrete(model, data_list, settings, update_ind, par[[i]], i)
    }
  }

  for(i in 1:settings[["numberChainRuns"]]) {
    out_post <- out_raw[[i]][["output"]][, 1:settings$numberFittedPar]
    outPTpar[[i]] <- out_raw[[i]][["PTMCpar"]]
    if (settings$numberFittedPar > 1){
        colnames(out_post) <- model[["namesOfParameters"]]
    }
    outPTpost[[i]] <- mcmc(out_post)
    outPTmissed[[i]] <- out_raw[[i]][["missed"]]
    outPTlp[[i]] <- out_raw[[i]][["output"]][, settings$numberFittedPar + 1]
    outPTtemp[[i]] <- out_raw[[i]][["output"]][, settings$numberFittedPar + 2]
    outPTacc[[i]] <- out_raw[[i]][["output"]][, settings$numberFittedPar + 3]
  }

  outlpv <- data.frame(matrix(unlist(outPTlp), nrow = length(outPTlp[[1]])))
  colnames(outlpv) <- c(1:settings[["numberChainRuns"]])
  outlpv <- outlpv %>% gather(colnames(outlpv), key="chain_no",value="lpost")
  outlpv$sample_no <-rep(1:length(outPTlp[[1]]), settings[["numberChainRuns"]])

  outltempv <- data.frame(matrix(unlist(outPTtemp), nrow=length(outPTtemp[[1]])))
  colnames(outltempv) <- c(1:settings[["numberChainRuns"]])
  outltempv <- outltempv %>% gather(colnames(outltempv), key="chain_no", value="temperature")
  outltempv$sample_no <- rep(1:length(outPTtemp[[1]]), settings[["numberChainRuns"]])
  
  outlaccv <- data.frame(matrix(unlist(outPTacc), nrow=length(outPTacc[[1]])))
  colnames(outlaccv) <- c(1:settings[["numberChainRuns"]])
  outlaccv <- outlaccv %>% gather(colnames(outlaccv), key="chain_no", value="acceptance rate")
  outlaccv$sample_no <- rep(1:length(outPTacc[[1]]), settings[["numberChainRuns"]])

  output <- list(
    mcmc = as.mcmc.list(outPTpost),
    missed = outPTmissed,
    lpost = outlpv,
    temp = outltempv,
    acc = outlaccv,
    outPTpar = outPTpar
  )
  output
}

check_settings_discete <- function(settings) {
  if (is.null(settings[["numberChainRuns"]])) {
    settings[["numberChainRuns"]] <- 4
    cat("`numberChainRuns` not specified in settings. Default value 4. \n")
  }

  if (is.null(settings[["numberCores"]])) {
    settings[["numberCores"]] <- settings[["numberChainRuns"]]
    cat("`numberCores` not specified in settings. Default value equal to `numberChainRuns`. \n")
  }

  if (is.null(settings[["numberTempChains"]])) {
    settings[["numberTempChains"]] <- 10
    cat("`numberTempChains` not specified in settings. Default value 10. \n")
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
        settings[["burninAdaptiveCov"]] <- 2000
    cat("`burninAdaptiveCov` not specified in settings. Default value 2000. \n")
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
        settings[["covarInitVal"]] <- 1e-10
    cat("`covarInitVal` not specified in settings.  Default value 1e-10. \n")
  }
  if (is.null(settings[["covarInitValAdapt"]])) {
        settings[["covarInitValAdapt"]] <- 1e-10
    cat("`covarInitValAdapt` not specified in settings.  Default value 1e-10. \n")
  }
  if (is.null(settings[["covarMaxVal"]])) {
        settings[["covarMaxVal"]] <- 1
    cat("`covarMaxVal` not specified in settings. Default value 1. \n")
  }
  if (is.null(settings[["runParallel"]])) {
        settings[["runParallel"]] <- TRUE
    cat("`runParallel` not specified in settings. Default value TRUE. \n")
  }

  if (is.null(settings[["lengthDiscreteVec"]])) {
    stop("`lengthDiscreteVec` not specified in settings. MUST be specified. \n")
  }
  if (is.null(settings[["updateDiscreteFreq"]])) {
        settings[["updateDiscreteFreq"]] <- 0
    cat("`updateDiscreteFreq` not specified in settings. Default value 0. \n")
  }

  settings
}