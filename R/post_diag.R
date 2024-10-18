#' @useDynLib serojump
#' @importFrom Rcpp sourceCpp
#' @import ggplot2
#' @import purrr
#' @import tidyr
#' @import dplyr
#' @import ggdist
#' @import patchwork
#' @import tidybayes
#' @import furrr
#' @import utils
#' @importFrom magrittr %>% %<>%
NULL



#' Postprocess the posteriors of the simulation study
#' @param model_summary The fitted model from serojump
#' @param save_info Information on how to save the diagnosis figures
#' @export 
plotMCMCDiagnosis <- function(model_summary, save_info) {
    
    check_save_info(save_info)

    if (model_summary$fit$data$priorPredFlag) {
        dir.create(here::here("outputs", "fits", save_info$file_name, save_info$model_name,  "figs_pp", "diag"), showWarnings = FALSE, recursive = TRUE)
        file_path <- here::here("outputs", "fits", save_info$file_name, save_info$model_name,  "figs_pp", "diag")
    } else{
        dir.create(here::here("outputs", "fits", save_info$file_name, save_info$model_name,  "figs", "diag"), showWarnings = FALSE, recursive = TRUE)
        file_path <- here::here("outputs", "fits", save_info$file_name, save_info$model_name,  "figs", "diag")
    }

    df_smi_df <- calcScaleModelIndicator(model_summary)
    transDimConvPlot(df_smi_df, file_path) 
    invariantParamConvPlot(model_summary, file_path) 
    invariantParamConvPriorPlot(model_summary, file_path)
    plotRhatTime(model_summary, file_path)
}


plotRhatTime <- function(model_summary, file_path) {

    outputfull <- model_summary$post

    model_outline <- model_summary$fit$model
    bio_all <- model_outline$infoModel$biomarkers

    fit_states_dt <- as.data.table(outputfull$fit_states)
    S <- fit_states_dt %>% filter(id == 1) %>% nrow

    ids <- fit_states_dt %>% group_by(id) %>% summarise(prob = sum(inf_ind) / S) %>% filter(prob > 0.5) %>% pull(id) %>% unique

    if (length(ids) == 0) {
        cat("No individuals have posterior prob of infection > 0.5")
    } else {
        # extract values here
        df_mcmc_time <- fit_states_dt %>% filter(id %in% ids) %>% filter(inf_ind == 1) %>% 
            select(id, chain_no, sample, inf_time, !!bio_all) %>% rename(chain = chain_no) 

        df_mcmc_time_wide <- df_mcmc_time %>% 
            select(id, sample, chain, inf_time) %>% unique %>%
            pivot_wider(!chain, names_from = "id", values_from = "inf_time") 

        cols <- ncol(df_mcmc_time_wide)

        df_summary_disc <- 
                map_df(2:cols,
            ~df_mcmc_time_wide %>% select(sample, .x) %>% drop_na %>% summarise_draws() %>% .[2, ]
        )

        p1 <- df_mcmc_time %>% 
            ggplot() +
                stat_pointinterval(aes(x = inf_time, y = as.character(id), color = as.character(chain)), 
                    position = position_dodge(0.4)) + theme_bw() + 
                    labs(x = "Time in study", y = "ID", color = "Chain number") 

        p2 <- df_summary_disc %>% ggplot() + geom_col(aes(x = rhat, y = as.character(variable))) + theme_bw() + 
            geom_vline(xintercept = 1.1, color = "red", linetype = "dashed") + 
            labs(x = "Rhat", y = "ID")

        p1 + p2
        ggsave(here::here(file_path, "timing_convergence.png"), height = 10, width = 10)

    }

}

#https://www.tandfonline.com/doi/abs/10.1198/1061860031347
#https://www.semanticscholar.org/reader/c5a24c1fcafcc80ec6fd22489585b13a0c3643de#
calcScaleModelIndicator <- function(model_summary) {

    fit <- model_summary$fit

    C <- length(fit$post$jump)
    dims <- dim(fit$post$jump[[1]])
    cat(str(fit$post$jump[[1]]))
    M <- dims[1]
    N <- dims[2]

    # Initializing the scaleModel
    log2_scaleModel <- 0
    df_smi_df <- map_df(1:C, 
        function(c) 
        {
            sMI <- map_dbl(1:M, 
                function(i) {
                    theta_i <- fit$post$jump[[c]][i, ]
                    theta_i_D <- sum(theta_i > -1)
                    theta_i_trim <- theta_i[theta_i > -1]
                    # Compute the terms vector
                    terms <- theta_i_trim * 2 ^ (0:(theta_i_D - 1))
                    cumulative_sums <- cumsum(terms)
                    log2_scaleModel <- log2(1 + cumulative_sums[theta_i_D])
                }
            )
        dims <- map_dbl(1:M, 
                function(i) {
                    theta_i <- fit$post$jump[[c]][i, ]
                    theta_i_D <- sum(theta_i > -1)
                }
            )
            data.frame(
                .chain = c,
                .iteration = 1:M,
                sMI = sMI,
                dims = dims
            )
        }
    )
    df_smi_df
}


transDimConvPlot <- function(df_smi_df, file_path) {

    pdims_trace <- df_smi_df %>% 
        ggplot() + 
            geom_line(aes(x = .iteration, y = dims, color = as.character(.chain))) + 
            labs(x = "Iteration", y = "Model dimension", color = "Chain") 


    pdims_hist <- df_smi_df %>% 
        ggplot() + 
            geom_histogram(aes(x = dims, fill = as.character(.chain))) + 
            labs(x = "Model dimension", y = "Count", fill = "Chain") 

    p1 <- pdims_trace + pdims_hist + plot_layout(guides = "collect") +
        plot_annotation(title = "Transdimensional convergence: dimensions of model") & theme_minimal() & theme(legend.position = "top")



    psmi_trace <- df_smi_df %>% 
        ggplot() + 
            geom_line(aes(x = .iteration, y = sMI, color = as.character(.chain))) + 
            labs(x = "Iteration", y = "Log2 of SMI", color = "Chain") 


    psmi_hist <- df_smi_df %>% 
        ggplot() + 
            geom_histogram(aes(x = sMI, fill = as.character(.chain))) + 
            labs(x = "Log2 of SMI", y = "Count", fill = "Chain") 

    p2 <- psmi_trace + psmi_hist  + plot_layout(guides = "collect") & theme_minimal() & theme(legend.position = "top")

    p3 <- summarise_draws(df_smi_df) %>% select(variable, rhat, ess_bulk, ess_tail) %>% pivot_longer(!variable, names_to = "stat", values_to = "value") %>%
        filter(variable == "sMI") %>% 
        ggplot()+ 
            geom_col(aes(y = stat, x = value)) + 
            facet_grid(rows = vars(stat), scales = "free") + theme_minimal()


    p1 / p2 / p3 +
        plot_annotation(title = "Transdimensional convergence: dimensions/SMI of chain")

    ggsave(here::here(file_path, "transdim_conv.png"), height = 20)
}

invariantParamConvPlot <- function(model_summary, file_path) {

    fit <- model_summary$fit
    p1 <- (fit$post$mcmc  %>% bayesplot::mcmc_trace()) + theme_minimal() + theme(legend.position = "top")
    p2 <- fit$post$lpost %>% ggplot() + geom_line(aes(x = sample_no, y = lpost, color = chain_no))  + theme_minimal() + theme(legend.position = "top")

    p3 <- df_conver_stat <- summarise_draws(fit$post$mcmc ) %>% select(variable, rhat, ess_bulk, ess_tail) %>% 
        pivot_longer(!variable, names_to = "stat", values_to = "value") %>%
        ggplot() + 
            geom_col(aes(y = variable, x = value)) +
            facet_wrap(~stat, scales = "free") + theme_minimal()

    p1 / p2 / p3
    ggsave(here::here(file_path, "invariant_param_conv.png"), height = 20)
}

invariantParamConvPriorPlot <- function(model_summary, file_path) {
    fit <- model_summary$fit 
    model_outline <- fit$model
    post <- fit$post

    mcmc_combined <- post$mcmc %>% lapply(as.data.frame) %>% do.call(rbind, .) %>% as.data.frame 
    n_post <- nrow(mcmc_combined)

    priorpost <- bind_rows(
        mcmc_combined %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>%
                        mutate(type = "Posterior distribution") ,
        purrr::map_df(1:n_post,
            ~model_outline$samplePriorDistributions(par_tab)
        )  %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>%  mutate(type = "Prior distribution") 
    )

    priorpost %>% 
        ggplot() + geom_histogram(aes(x = value, fill = type), position = "identity", alpha = 0.5) + facet_wrap(~param, scales = "free") +
        theme_bw() + theme(legend.position = "top")
    
    ggsave(here::here(file_path, "invariant_param_priorpost.png"), height = 15, width = 15)

}
           