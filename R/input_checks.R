addExposurePrior_checkfunction <- function(exp_prior) {
    columns <- c("lb", "ub", "func", "par1", "par2")
    existing_columns <- names(exp_prior)

    non_existing_columns <- setdiff(columns, existing_columns)
    if (length(non_existing_columns) > 0) {
        stop(paste("Error: Column(s)", paste(non_existing_columns, collapse = ", "), "do not exist in the exposure prior dataframe"))
    }
}

addExposurePrior_checkempirical <- function(exp_prior, data_t) {
    columns <- c("day", "prob")
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