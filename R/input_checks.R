addExposurePrior_checkfunction <- function(exp_prior) {
    columns <- c("lb", "ub", "func", "par1", "par2")
    existing_columns <- names(exp_prior)

    non_existing_columns <- setdiff(columns, existing_columns)
    if (length(non_existing_columns) > 0) {
        stop(paste("Error: Column(s)", paste(non_existing_columns, collapse = ", "), "do not exist in the exposure prior dataframe"))
    }
}

check_titre <- function(data_titre) {
    check_time(data_titre)
}

check_time <- function(data_titre) {
    
    time_diff_too_small <- data.frame(
        id = 1:(data_titre$id %>% max),
        start_time = data_titre %>% group_by(id) %>% filter(time == min(time)) %>% unique %>% pull(time),
        end_time = data_titre %>% group_by(id) %>% filter(time == max(time)) %>% unique %>% pull(time)
    ) %>% mutate(time_diff = end_time - start_time) %>% filter(time_diff <= 7)

    if(nrow(time_diff_too_small) > 0) {
        stop(paste("Error: Time difference between first and last titre measurement is less than 7 days for the following id(s):", paste(time_diff_too_small$id, collapse = ", "), ". Please remove these values and recode the data inputs."))
    }
}


addExposurePrior_checkempirical <- function(exp_prior, data_t) {
    columns <- c("time", "prob")
    existing_columns <- names(exp_prior)

    non_existing_columns <- setdiff(columns, existing_columns)
    if (length(non_existing_columns) > 0) {
        stop(paste("Error: Column(s)", paste(non_existing_columns, collapse = ", "), "do not exist in the exposure prior dataframe"))
    }
    column_types <- sapply(exp_prior[, columns, drop = FALSE], class)
    if (!all(column_types %in% c("numeric", "integer"))) {
        stop(paste("Error: Column(s)", paste(names(column_types[!column_types %in% c("numeric", "integer")]), collapse = ", "), "are not numeric or integers in the exposure prior dataframe"))
    }

    if (data_t$T != nrow(exp_prior) ) {
        warning("Warning: The number of rows in the exposure prior dataframe does not match the number of time points in the study\n")
    }
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

 
    # CHECK BIOMARKERS ARE WELL DEFINED
    # Check columns of data_sero match model definition 
    data_sero_name <- data_sero %>% names
    biomarkers_md <- modeldefinition$biomarkers
    biomarkers_obs <- modeldefinition$observationalModel$model %>% map(~.x$biomarker) %>% unlist
    biomarkers_abkin <- modeldefinition$abkineticsModel$model %>% map(~.x$biomarker) %>% unlist %>% unique

    for(b in biomarkers_md) {
        if(!b %in% data_sero_name) {
            stop("Biomarker, ", b, ", in `modeldefinition$biomarkers` is not a column of serological data; `",
                paste(data_sero_name, collapse = ", "), "`")
        }
    }
     if(!identical(biomarkers_md, biomarkers_obs) ) {
        stop("Biomarkers in observationalModel (", paste(biomarkers_obs, collapse = ", "), ") do not match biomarkers in `modeldefinition$biomarkers` (", paste(biomarkers_md, collapse = ", "), ")")
    }
    if(!identical(biomarkers_md, biomarkers_abkin) ) {
        stop("Biomarkers in abkineticsModel (", paste(biomarkers_abkin, collapse = ", "), ") do not match biomarkers in `modeldefinition$biomarkers` (", paste(biomarkers_md, collapse = ", "), ")")
    }  

    exposures_md <- modeldefinition$exposureTypes
    exposures_obs <- modeldefinition$abkineticsModel$model %>% map(~.x$exposureType) %>% unlist %>% unique

    # CHECK EXPSURETYPES ARE WELL DEFINED
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

}
