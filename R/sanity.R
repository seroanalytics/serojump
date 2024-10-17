#' @title A plotting function for the serological data
#' @param seroModel A serological model object
#' @return A ggplot object with the plotted data
#' @description This function takes a serological model object and plots the serological data. The plot gives an overview of the data across individuals
#' @export
plotSero <- function(seroModel, N_max = 100){
    N <- seroModel$data$N

    if (N > N_max) {
        sero_plot <- seroModel$data$raw_sero %>% filter(id <= N_max)  %>% mutate(id = factor(id, levels = unique(id))) 
    } else {
        sero_plot <- seroModel$data$raw_sero %>% mutate(id = factor(id, levels = unique(id))) 
    }

    T <- sero_plot %>% pull(time) %>% max

    start_bleed <- sero_plot %>% group_by(id) %>% filter(time == min(time)) 
    end_bleed <- sero_plot %>% group_by(id)   %>% filter(time == max(time)) %>%
      mutate(time = time - 7) 

    sero_plot %>% 
        ggplot2::ggplot() + 
        geom_tile(data = start_bleed, aes(y = id, x = time / 2, width = time ), alpha = 0.9, fill = "gray90") +
        geom_tile(data = end_bleed, aes(y = id, x = time + (T - time) / 2, width = (T - time) ), alpha = 0.9, fill = "gray90") +
        geom_point(aes(x = time, y = id), shape = "|", size = 2.5) + 
      #  geom_point(data = epi_times_sims, aes(x = t, y = i, color = exp_type) ,  shape = 4) 
        labs(x = "Time in study (days)", y = "Id of individual", color = "Exposure type") + 
        geom_hline(yintercept = seq(0.5, N + 0.5, 1), color = "gray90") +
        theme(panel.grid.major.y = element_blank()) + theme_bw()
}

#' @title A plotting function for the prior predictive distribution of the antibody kinetics
#' @param seroModel A serological model object
#' @return A ggplot object with the plotted data
#' @description This function takes a serological model object and plots the prior predictive distribution of the antibody kinetics. The plot gives an overview of the prior predictive distribution of the antibody kinetics
#' @export
plotPriorPredictive <- function(seroModel) {
    abModels <- seroModel$model$abkineticsModel 

    M <- length(abModels)

    full_pp_ab <- map_df(
        1:M,
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
    full_pp_ab_sum <- full_pp_ab %>% group_by(t, name) %>% mean_qi(titre, .width = c(0.5, 0.95))

    full_pp_ab_sum %>%
        ggplot2::ggplot() + 
            geom_ribbon(data = full_pp_ab_sum %>% filter(.width == 0.95), aes(x = t, ymin = .lower, ymax = .upper), alpha = 0.3) + 
            geom_ribbon(data = full_pp_ab_sum %>% filter(.width == 0.5), aes(x = t, ymin = .lower, ymax = .upper), alpha = 0.7) + 
            geom_line(aes(x = t, titre), size = 3) + 
            facet_wrap(vars(name)) + 
            theme_minimal() + theme(text = element_text(size = 20)) + 
            labs(x = "Time post exposure", y = "Biomarker value")
}

#' @title A plotting function for the prior on the infection rate
#' @param seroModel A serological model object
#' @return A ggplot object with the plotted data
#' @description This function takes a serological model object and plots the prior distribution of the exposure rate both on the population-level and individual-level.
#' @export
plotPriorInfection <- function(seroModel) {

    p1 <- seroModel$data$exp_prior %>% 
        ggplot2::ggplot() + 
            geom_line(aes(x = time, prob), size = 3) + 
            theme_minimal() + theme(text = element_text(size = 20)) + 
            labs(x = "Time in study", y = "PDF of infection") + 
            ylim(0, NA)

    exp_full_df <- map_df(1:length(seroModel$data$exp_list), 
        function(i) {
            data.frame(
                time = 1:length(seroModel$data$exp_list[[i]]),
                prob = seroModel$data$exp_list[[i]],
                id = i
            )
        }
    ) 
    p2 <- exp_full_df %>% 
        ggplot() + 
            geom_tile(aes(x = time, y = id, fill = prob), alpha = 0.9) 

    library(patchwork)
    p1 / p2
}