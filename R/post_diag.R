#' @useDynLib rjmc
#' @importFrom Rcpp sourceCpp
#' @import ggplot2
#' @import purrr
#' @import tidyr
#' @import dplyr
#' @import ggdist
#' @import patchwork
#' @import tidybayes
#' @import furrr
#' @importFrom magrittr %>% %<>%
NULL



#' Postprocess the posteriors of the simulation study
#' @param modelname The nams of the model
#' @param modelname_sim The name of the simulation model
#' @param obs_er The observation error
#' @param n_chains The number of chains
#' @export 
seroJumpPostDaig <- function(fitfull, filepathfig) {
    df_smi_df <- calcScaleModelIndicator(fitfull)
    transDimConvPlot(df_smi_df, filepathfig) 
    invariantParamConvPlot(fitfull, filepathfig) 
}

#https://www.tandfonline.com/doi/abs/10.1198/1061860031347
#https://www.semanticscholar.org/reader/c5a24c1fcafcc80ec6fd22489585b13a0c3643de#
calcScaleModelIndicator <- function(fitfull, filepathfig) {

    C <- length(fitfull$post$jump)
    cat(C)
    dims <- dim(fitfull$post$jump[[1]])
    cat(str(fitfull$post$jump[[1]]))
    M <- dims[1]
    N <- dims[2]

    # Initializing the scaleModel
    log2_scaleModel <- 0
    df_smi_df <- map_df(1:C, 
        function(c) 
        {
            sMI <- map_dbl(1:M, 
                function(i) {
                    theta_i <- fitfull$post$jump[[c]][i, ]
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
                    theta_i <- fitfull$post$jump[[c]][i, ]
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


transDimConvPlot <- function(df_smi_df, filepathfig) {


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
    ggsave(here::here(filepathfig, "diag", "transdim_conv.png"), height = 20)
}

invariantParamConvPlot <- function(fitfull, filepathfig) {

    p1 <- fitfull$post$mcmc  %>% mcmc_trace + theme_minimal() + theme(legend.position = "top")
    p2 <- fitfull$post$lpost %>% ggplot() + geom_line(aes(x = sample_no, y = lpost, color = chain_no))  + theme_minimal() + theme(legend.position = "top")

    p3 <- df_conver_stat <- summarise_draws(fitfull$post$mcmc ) %>% select(variable, rhat, ess_bulk, ess_tail) %>% 
        pivot_longer(!variable, names_to = "stat", values_to = "value") %>%
        ggplot() + 
            geom_col(aes(y = variable, x = value)) +
            facet_wrap(~stat, scales = "free") + theme_minimal()

    p1 / p2 / p3
    ggsave(here::here(filepathfig, "diag", "invariant_param_conv.png"), height = 20)
}
