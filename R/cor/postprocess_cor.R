
plot_cop_recInfstan <- function(outputfull, fitfull, fig_folder, scale_ab = NULL) {
    
    require(rstan)

    fit_states <- outputfull$fit_states
    filename <- outputfull$filename
    modelname <- outputfull$modelname

    post <- fitfull$post
    data_t <- fitfull$data_t

    n_chains <- outputfull$n_chains
    n_post <- outputfull$n_post
    n_length <- n_chains * n_post
    chain_samples <- 1:n_chains %>% map(~c(rep(.x, n_post))) %>% unlist

    model_outline <- fitfull$model

    post_fit <- post$mcmc %>% lapply(as.data.frame) %>% do.call(rbind, .) %>% as.data.frame %>% mutate(chain = as.character(chain_samples ))

    n_post <- outputfull$n_post
    n_length <- n_chains * n_post
  
    cop_exp_sum_plot_all <- map_df(1:length(model_outline$observationalModel), 
        function(i) {
        biomarker <- model_outline$observationalModel[[i]]$biomarker

        df_full_info <- fit_states %>% mutate(titre_type = case_when(
            inf_time == -1 ~ "No inf",
            inf_time != -1 ~ "Inf")
        ) %>%  group_by(id, titre_type) %>% add_count(name = "n_rows") %>% 
        group_by(id, titre_type, n_rows) %>% mutate(prop = n_rows / (n_length)) %>%
         group_by(id, titre_type, prop) %>%
        mean_qi(!!sym(biomarker) ) %>% 
        complete(id = 1:data_t$N, titre_type, fill = list(prop = 0)) %>%
        arrange(!!biomarker) %>% as.data.frame 

        # Get response for each individual
        df_full_info_res <- df_full_info %>% filter(titre_type == "Inf") %>% select(id, prop)

        # Get covariate for each individual
        df_full_info_x <- df_full_info %>% group_by(id) %>% mutate(!!sym(biomarker) := weighted.mean(x = !!sym(biomarker), w = prop)) %>% 
         select(id, !!sym(biomarker))

        df_data <- data.frame(
            id = 1:data_t$N,
            start_titre = data_t$initialTitreValue %>% as.data.frame %>% pull(!!biomarker),
            known_inf = data_t$knownInfsVec
        ) %>% mutate(known_inf = ifelse(known_inf == 0, "Not known", "Known"))


        df_cor_info <- df_full_info_res %>% left_join(df_full_info_x) %>% left_join(df_data)  %>%
            rename(titre_val = !!(biomarker)) %>% mutate(biomarker = !!biomarker) 


        ## here
        #titre_cop_sum <- fit_states %>% filter(!!sym(biomarker) != -1) %>%
        #    group_by(id) %>% mean_qi(!!sym(biomarker) ) %>% arrange(!!biomarker) %>%
        #    rename(titre_val = !!(biomarker))
        #u_ids <- titre_cop_sum$id

        #df_data_post <- fit_states %>% filter(inf_ind == 1) %>% select(id, inf_ind) %>%
        #         summarise(inf_post = mean(inf_ind), n = n(), .by = "id") %>%
        #        left_join(df_data) %>% mutate(prop = n / n_length)

        #cop_exp_sum_plot <- titre_cop_sum %>% left_join(df_data_post) %>% mutate(id = factor(id, levels = u_ids)) %>%
         #   filter(!is.na(n)) %>% mutate(biomarker = !!biomarker)
        # here
        })


        # Compile the model
    model_code <- here::here("src", "stan", "normal.stan")
    sm <- stan_model(file = model_code)
    
    # split data up by biomarker
    cop_exp_sum_plot_all_split <- map(model_outline$infoModel$biomarker,
        ~cop_exp_sum_plot_all %>% 
        filter(biomarker == .x)
    )

    # Get split draws
    spread <- map(cop_exp_sum_plot_all_split, 
        ~ seq(min(.x$titre_val), max(.x$titre_val), length.out = 100) 
    )

    make_data <- function(data_i, spread) {
         stan_data <- list(
        N = data_i %>% nrow,
        X = data_i$titre_val,
        y = data_i$prop,
        N_val = 100,
        X_val = spread,
        lb_alpha = min(spread)

        )
    }

    stan_data <- map2(cop_exp_sum_plot_all_split, spread, make_data)

    # Fit the model
    fit_cor <- map(stan_data, ~sampling(sm, data = .x, iter = 2000, chains = 4, cores = 4) )
   # fig_folder <- "post"

    # PLot the mcmc traces
    biomarker_file <- model_outline$infoModel$biomarker %>% map(~gsub("/", "_", .x)) %>% unlist

    map(1:length(fit_cor), 
        function(i) {
        p1 <- fit_cor[[i]] %>% as.array(c("beta_i", "alpha", "gamma", "p_base", "k", "sigma")) %>% mcmc_trace
        p2 <- fit_cor[[i]] %>% as.array(c("beta_i", "alpha", "gamma", "p_base", "k", "sigma")) %>% mcmc_dens_overlay
        p1 / p2
        dir.create(here::here("outputs", "fits", filename,  "figs", modelname, fig_folder, "cor"), showWarnings = FALSE)
        ggsave(here::here("outputs", "fits", filename,  "figs", modelname, fig_folder, "cor", paste0("fitted_trace_", biomarker_file[i], ".png")), height = 10, width = 10)
        }
    )


    # Extract the gradient and posterior predictions
    y_predictions <- 
        map_df(1:length(fit_cor), 
        function(i) {
            rstan::extract(  fit_cor[[i]])$y_pred_mean %>% as.data.frame %>% rename_with(
                .fn = ~ paste0(spread[[i]], ""),  # Create new names
                .cols = starts_with("V")  # Apply to all columns starting with "V"
            ) %>% mutate(sample = 1:4000) %>% pivot_longer(!sample, names_to = "id", values_to = "prop") %>%
            mutate(biomarker = model_outline$infoModel$biomarker[i])
        }
        ) %>% mutate(id = as.numeric(id)) %>% group_by(sample, biomarker) %>% mutate(prop_rel = prop / max(prop))

    require(ggrepel)

    y_predictions_plot <- y_predictions %>%
        left_join(cop_exp_sum_plot_all %>% select(id, titre_val, biomarker) %>% distinct) %>% 
        group_by(biomarker) %>%
        mutate(label = if_else(id == min(id), as.character(biomarker), NA_character_)) 

    p1 <- y_predictions %>%
        ggplot() +
        stat_lineribbon(aes(x = id, prop), .width = 0.5, alpha = 0.5, fill = "gray") + 
        geom_point(data = cop_exp_sum_plot_all, aes(x = titre_val, y = prop), shape = 1, alpha = 0.4) +
        facet_grid(cols = vars(biomarker)) + 
        ylim(0, 1) + theme_bw() + 
        labs(x = "Titre at exposure", y = "Posterior probability of infection") 

    y_predictions_plot_label <- y_predictions %>% 
        group_by(id, biomarker) %>%
        summarise(prop = mean(prop), prop_rel = mean(prop_rel)) %>%
        ungroup %>% 
        group_by(biomarker) %>%
        mutate(label = if_else(id == min(id), as.character(biomarker), NA_character_)) 


    p2 <- y_predictions_plot_label %>%
         ggplot() +
        stat_lineribbon(aes(x = id, prop, fill = biomarker, color = biomarker), .width = 0.01, size = 6, alpha = 0.9) + 
        theme_bw() + 
        geom_label_repel(data = y_predictions_plot_label, aes(x = id, y = prop, color = biomarker, label = label), size = 7,
                  nudge_x = 0,
                  na.rm = TRUE) + 
        guides(color = FALSE, fill = FALSE) +
        ylim(0, NA) + 
        labs(x = "Titre at exposure", y = "Posterior probability of infection", color = "Biomarker", fill = "Biomarker")  

    p1 / p2

    ggsave(here::here("outputs", "fits", filename,  "figs", modelname, fig_folder, "cor", paste0("fitted_biomarker_", "combined", ".png")), height = 10, width = 10)


}