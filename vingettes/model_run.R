library(furrr)
run_rjmc <- function(model_info, save_info,  rj_settings, runmodel = TRUE) {
    
    if(runmodel) {
        model_summary <- runInfRJMCMC(model_info, rj_settings, save_info = save_info)
    } else {
        model_summary <- readRDS(here::here("outputs", "fits", 
                save_info$file_name, save_info$model_name, paste0("model_summary.RDS")))
    }

    # Need to have save_info model_summary to run these
    plotMCMCDiagnosis(model_summary, save_info = save_info)
    plotPostFigs(model_summary, save_info = save_info)
}

models <- readRDS(here::here("data", "fudan", "model_A.RDS"))


plan(multisession, workers = 4)

rj_settings <- list(
        numberChainRuns = 4, 
        iterations = 5000,
        burninPosterior = 2500,
        thin = 10
    )


#run_rjmc(models$model[[1]], models$save[[1]], rj_settings)
run_rjmc(models$model[[2]], models$save[[2]], rj_settings)


