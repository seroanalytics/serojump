filename <- "hpc/transvir_w2_inf/p3"
n_chains <- 4
require(readxl)
require(lubridate)


gambia_pvnt_w2 <- get_data_titre_model_wave2_pvnt_iga()# This is the empirical prior for the exposure time
gambia_exp_w2 <- get_exposures_wave2() # This is the empirical prior for the exposure time
exp_prior_w2 <- get_exp_prior_wave2() # This is the empirical prior for the exposure time

fitfull_w2 <- readRDS(here::here("outputs", "fits", filename, paste0("fit_", "w2", ".RDS")))
outputfull_w2 <- readRDS(file = here::here("outputs", "fits", filename, paste0("pp_", "w2", ".RDS")))

gambia_pvnt_w3 <- get_data_titre_model_wave3_pvnt_iga()# This is the empirical prior for the exposure time
gambia_exp_w3 <- get_exposures_wave3() # This is the empirical prior for the exposure time
exp_prior_w3 <- get_exp_prior_wave3() # This is the empirical prior for the exposure time

fitfull_w3 <- readRDS(here::here("outputs", "fits", filename, paste0("fit_", "w3", ".RDS")))
outputfull_w3 <- readRDS(file = here::here("outputs", "fits", filename, paste0("pp_", "w3", ".RDS")))

########################################################
################## FIGURES FOR FIGURE 1 ##################
########################################################

fitfull <- fitfull_w2
outputfull <- outputfull_w2
wave_name <- "Delta"

get_infer_df <- function(fitfull, outputfull, wave_name) {

    gambia_meta <- read.csv(file = here::here("data", "transvir","full_demog_data.csv") )
    gambia_meta_short <- gambia_meta %>% dplyr::select(pid = Participant_ID, age_cat, sex = Sex, vaccinated) 

    N <- fitfull$data_t$N
    known_inf <- fitfull$data_t$knownInfsN
    infer_inf <- outputfull$fit_states %>% group_by(sample) %>% filter(inf_ind == 1) %>% summarise(n_inf = n()) %>% pull(n_inf) %>%
        mean_qi %>% mutate(wave = wave_name) %>% mutate(across(c(y, ymin, ymax), ~ . / N))

    data_mean <- data.frame(
    wave = rep(wave_name, 2),
    Value = c(infer_inf$y, known_inf / N),
    Group = rep(c('Inferred infections', 'Known infections'), each = 1)
    ) %>% mutate(Group = factor(Group, levels = c('Inferred infections', 'Known infections')))
    list(infer_inf, data_mean)
}

get_infer_df_cov <- function(fitfull, outputfull, data_i, wave_name) {

    gambia_meta <- read.csv(file = here::here("data", "transvir","full_demog_data.csv") )
    gambia_meta_short <- gambia_meta %>% dplyr::select(pid = Participant_ID, age_cat, sex = Sex, vaccinated) 

    meta_df <- data_i %>% left_join(gambia_meta_short) %>% mutate(pid = factor(pid, levels = unique(pid))) %>%
        mutate(id = as.numeric(pid)) %>% select(id, age_cat, sex, vaccinated) %>% unique

    N <- fitfull$data_t$N
    meta_df$known_inf <- fitfull$data_t$knownInfsVec

    map_df(c("age_cat", "sex", "vaccinated"),
        function(covar) {
            known_inf <- meta_df %>% group_by(!!sym(covar)) %>% filter(known_inf == 1) %>% summarise(n_known = n())

            infer_inf <- outputfull$fit_states %>% group_by(sample) %>% filter(inf_ind == 1) %>% 
                left_join(meta_df) %>% group_by(sample, !!sym(covar)) %>% summarise(n_inf = n()) %>% 
                group_by(!!sym(covar)) %>% mean_qi(n_inf) %>% mutate(wave = wave_name) 
            

            infer_inf %>% left_join(known_inf) %>% mutate(prop = n_known / n_inf) %>% 
                mutate(error = sqrt((prop * (1 - prop)) / n_inf)) %>% mutate(covariate = covar) %>% 
                rename(level = !!sym(covar))
        }
    ) %>% mutate(wave = wave_name)
}
    


df_delta <- get_infer_df(fitfull_w2, outputfull_w2, "Delta")
df_omicron <- get_infer_df(fitfull_w3, outputfull_w3, "Omicron")

infer_inf <- bind_rows(df_delta[[1]], df_omicron[[1]])
data_mean <- bind_rows(df_delta[[2]], df_omicron[[2]])

df_delta_cov <- get_infer_df_cov(fitfull_w2, outputfull_w2, gambia_pvnt_w2, "Delta")
df_omicron_cov <- get_infer_df_cov(fitfull_w3, outputfull_w3,gambia_pvnt_w3,  "Omicron")

df_both <- bind_rows(df_delta_cov, df_omicron_cov)

data_mean %>% 
    ggplot() +
  geom_bar(aes(x = wave, y = Value, fill = Group), stat = "identity", position = position_dodge(width = 0.0), alpha = 1) +
  geom_errorbar(data = infer_inf, aes(x = wave, ymin = ymin, ymax = ymax), width = 0.4, size = 2) + 
  scale_fill_manual(values = c("#cf725e", "#62656e")) +
  labs(title = "Known and inferred attack rates",
       x = "SARS-CoV-2 wave",
       y = "Attack rate", fill = "") +
    theme_minimal() +
    theme(text = element_text(size = 15, color = "black")) 

df_both %>% ggplot() + 
    geom_col(aes(x = level, y = 1 - prop)) + 
    facet_grid(cols = vars(covariate), rows = vars(wave), space = "free_x", scale = "free_x") + 
    theme_bw()



########################################################
################## FIGURES FOR FIGURE 2 ##################
########################################################
# Prior pred
gamma <- rnorm(10000, 0, 1.68)
z_gamma <- rnorm(10000, 0, 1)
gamma_sigma <- rexp(10000, 3)
inv.logit(gamma + z_gamma * gamma_sigma)%>% hist#+ z_gamma[X_g1[n]] * gamma_sigma


# COP STUFF
require(rstan)
run_all_mcmc <- function(fitfull, outputfull, data_i, wave_name) {
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
            select(id, !!sym(biomarker)) %>% unique

        df_data <- data.frame(
            id = 1:data_t$N,
            start_titre = data_t$initialTitreValue %>% as.data.frame %>% pull(!!biomarker),
            known_inf = data_t$knownInfsVec
        ) %>% mutate(known_inf = ifelse(known_inf == 0, "Not known", "Known"))


        df_cor_info <- df_full_info_res %>% left_join(df_full_info_x) %>% left_join(df_data)  %>%
            rename(titre_val = !!(biomarker)) %>% mutate(biomarker = !!biomarker) 

    })

    gambia_meta_short <- gambia_meta %>% dplyr::select(pid = Participant_ID, age_cat, sex = Sex, vaccinated) 

    meta_df <- data_i %>% left_join(gambia_meta_short) %>% mutate(pid = factor(pid, levels = unique(pid))) %>%
        mutate(id = as.numeric(pid)) %>% select(id, age_cat, sex, vaccinated) %>% unique

    df_fit_all <- cop_exp_sum_plot_all %>% left_join(meta_df)

        # Compile the model
    model_code <- here::here("src", "stan", "normal.stan")
    sm <- stan_model(file = model_code)

    model_code <- here::here("src", "stan", "normal_cov.stan")
    sm_cov <- stan_model(file = model_code)

    # split data up by biomarker
    cop_exp_sum_plot_all_split <- map(model_outline$infoModel$biomarker,
        ~df_fit_all %>% 
        filter(biomarker == .x)
    )

    run_mcmc_cor <- function(curr_data) {

        # Get split draws
        spread <- seq(min(curr_data$titre_val), max(curr_data$titre_val), length.out = 100) 
        
        X_g1 <- curr_data %>% mutate(age_cat_i = 
            case_when(age_cat == "<5"~1,
            age_cat == "5-17"~2,
            age_cat == "18-49"~3,
            age_cat == ">50"~4
            )) %>% pull(age_cat_i)

        N_g1 <- unique(X_g1) %>% length

        X_g1 <- curr_data %>% mutate(vaccinated_i = 
            case_when(vaccinated == "No"~1,
            vaccinated == "Yes"~2
            )) %>% pull(vaccinated_i)

        N_g1 <- unique(X_g1) %>% length
        prop_g1 <- X_g1 %>% table %>% `/`(., sum(.)) %>% as.numeric

        make_data <- function(data_i, spread) {
            stan_data <- list(
                N = curr_data %>% nrow,
                X = curr_data$titre_val,
                y = curr_data$prop,
                N_val = 100,
                X_val = spread,
                lb_alpha = min(spread),
                ub_alpha = max(spread),
                X_g1 = X_g1,
                N_g1 = N_g1,
                prop_g1 = prop_g1
            )
        }

        stan_data <- make_data(curr_data, spread)
        # Fit the model
        fit_cor <- sampling(sm_cov, data = stan_data, iter = 2000, chains = 4, cores = 4) 

    }

    curr_data_spike <- cop_exp_sum_plot_all_split[[1]]
    curr_data_iga <- cop_exp_sum_plot_all_split[[2]]

    fit_spike <- run_mcmc_cor(curr_data_spike)
    fit_iga <- run_mcmc_cor(curr_data_iga)


    fit_cor <- list(fit_spike, fit_iga)
    # fig_folder <- "post"

    all_mean_traj <- map_df(
        1:2,
        function(i) {
            y_predictions_mean <- rstan::extract(fit_cor[[i]])$y_pred_mean %>% as.array %>% as.data.frame %>% rename_with(
                        .fn = ~ paste0(spread, ""),  # Create new names
                        .cols = starts_with("V")  # Apply to all columns starting with "V"
                    ) %>% mutate(sample = 1:4000) %>% pivot_longer(!sample, names_to = "id", values_to = "prop")  %>%
                    mutate(id = as.numeric(id)) %>% group_by(sample) %>% mutate(prop_rel = prop / max(prop)) %>%
                    mutate(biomarker = model_outline$infoModel$biomarker[i])
            y_predictions_mean
        }

    ) %>% mutate(wave = wave_name)


    curr_data <- bind_rows(curr_data_spike, curr_data_iga) %>%  mutate(wave = wave_name)

    list(fit = fit_cor, model = all_mean_traj, data = curr_data)
}


output_w2 <- run_all_mcmc(fitfull_w2, outputfull_w2, gambia_pvnt_w2, "Delta")
output_w3 <- run_all_mcmc(fitfull_w3, outputfull_w3, gambia_pvnt_w3, "Omicron")

k <- 1
alpha <- 1
beta <- -2
esttitreExp <- seq(0, 4, 0.1)
r <- beta * (esttitreExp - alpha);

(r / (1 + abs(r)^k)^(1 / k)) * 0.5 + 0.5

data.frame(
    value = c(rstan::extract(output_w2$fit[[1]])$beta_i_ref_mean_out, rstan::extract(output_w2$fit[[2]])$beta_i_ref_mean_out),
    wave = c(rep("spike", 4000), rep("IgA", 4000))
) %>% 
    ggplot() + geom_boxplot(aes(x = value, y = wave))

data.frame(
    value = c(rstan::extract(output_w2$fit[[1]])$alpha_ref_mean_out, rstan::extract(output_w2$fit[[2]])$alpha_ref_mean_out),
    wave = c(rep("spike", 4000), rep("IgA", 4000))
) %>% 
    ggplot() + geom_boxplot(aes(x = value, y = wave))


output_traj_combined <- bind_rows(
    output_w2$model, output_w3$model
)

y_predictions_plot_label <- output_traj_combined %>% 
    group_by(id, biomarker, wave) %>%
    summarise(prop = mean(prop), prop_rel = mean(prop_rel)) %>%
    ungroup %>% 
    group_by(biomarker) %>%
    mutate(label = if_else(id == min(id), as.character(biomarker), NA_character_)) 


y_predictions_plot_label %>%
        ggplot() +
    stat_lineribbon(aes(x = id, prop, fill = biomarker, color = biomarker), .width = 0.01, size = 6, alpha = 0.9) + 
    theme_bw() + 
    geom_label_repel(data = y_predictions_plot_label, aes(x = id, y = prop, color = biomarker, label = label), size = 7,
                nudge_x = 0,
                na.rm = TRUE) + 
    guides(color = FALSE, fill = FALSE) +
    ylim(0, NA) + 
    facet_grid(rows = vars(wave)) + 
    labs(x = "Titre at exposure", y = "Posterior probability of infection", color = "Biomarker", fill = "Biomarker")  

# PLot the mcmc traces
biomarker_file <- model_outline$infoModel$biomarker %>% map(~gsub("/", "_", .x)) %>% unlist

age_cats <- c("<5", "5-17", "18-49", ">50")
vaccination_val <- c( "No", "Yes")

y_predictions <- map_df(1:N_g1,
    function(i) {
    rstan::extract(fit_cor)$y_pred_mean_cov %>% as.array %>% .[, ,i] %>% as.data.frame %>% rename_with(
                .fn = ~ paste0(spread, ""),  # Create new names
                .cols = starts_with("V")  # Apply to all columns starting with "V"
            ) %>% mutate(sample = 1:4000) %>% pivot_longer(!sample, names_to = "id", values_to = "prop")  %>%
            mutate(vaccinated = vaccination_val[i])
    }
) %>% mutate(id = as.numeric(id)) %>% group_by(sample, vaccinated) %>% mutate(prop_rel = prop / max(prop))



require(ggrepel)

y_predictions_plot <- y_predictions %>%
    left_join(curr_data %>% select(id, titre_val, vaccinated) %>% distinct) %>% 
    group_by(vaccinated) %>%
    mutate(label = if_else(id == min(id), as.character(vaccinated), NA_character_)) 

p1 <- y_predictions %>%
    ggplot() +
    stat_lineribbon(aes(x = id, prop), .width = 0.5, alpha = 0.5, fill = "gray") + 
    geom_point(data = curr_data, aes(x = titre_val, y = prop), shape = 1, alpha = 0.4) +
    facet_grid(cols = vars(vaccinated)) + 
    ylim(0, 1) + theme_bw() + 
    labs(x = "Titre at exposure", y = "Posterior probability of infection") 

y_predictions_plot_label <- y_predictions %>% 
    group_by(id, vaccinated) %>%
    summarise(prop = mean(prop), prop_rel = mean(prop_rel)) %>%
    ungroup %>% 
    group_by(vaccinated) %>%
    mutate(label = if_else(id == min(id), as.character(vaccinated), NA_character_)) 


p2 <- y_predictions_plot_label %>%
        ggplot() +
    stat_lineribbon(aes(x = id, prop, fill = vaccinated, color = vaccinated), .width = 0.01, size = 6, alpha = 0.9) + 
    theme_bw() + 
    geom_label_repel(data = y_predictions_plot_label, aes(x = id, y = prop, color = vaccinated, label = label), size = 7,
                nudge_x = 0,
                na.rm = TRUE) + 
    guides(color = FALSE, fill = FALSE) +
    ylim(0, NA) + 
    labs(x = "Titre at exposure", y = "Posterior probability of infection", color = "Biomarker", fill = "Biomarker")  

p1 / p2
