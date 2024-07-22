filename <- "hpc/nih_2024_inf/p3"
n_chains <- 4
require(readxl)
require(lubridate)

nih_inf_raw <- read_excel(path = here::here("data", "nih_2024_XX", "2022_2023_Flu_Swabs.xlsx") )

data_exp_h1 <- get_data_titre_nih_2023_h1()
data_titre_h1 <- data_exp_h1[[1]] #%>% check_titre 
data_exp_h3 <- get_data_titre_nih_2023_h3()
data_titre_h3 <- data_exp_h3[[1]] #%>% check_titre 


fitfull_h1 <- readRDS(here::here("outputs", "fits", filename, paste0("fit_", "h1", ".RDS")))
outputfull_h1 <- readRDS(file = here::here("outputs", "fits", filename, paste0("pp_", "h1", ".RDS")))

n_post <- outputfull_h1$n_post
n_length <- n_chains * n_post
model_outline <- fitfull_h1$model
data_t <- fitfull_h1$data_t
fit_states <- outputfull_h1$fit_states
summarise_inf_prop_h1 <- map_df(1:length(model_outline$observationalModel), 
    function(i) {

    biomarker <- model_outline$observationalModel[[i]]$biomarker
    fit_states %>% filter(id == "194")
    titre_cop_sum <- fit_states %>% filter(!!sym(biomarker) != -1) %>%
        group_by(id) %>% mean_qi(!!sym(biomarker) ) %>% arrange(!!biomarker) %>%
        rename(titre_val = !!(biomarker))
    u_ids <- titre_cop_sum$id

    fit_states %>% group_by(id) %>% filter(inf_ind == 1) %>% summarise(n = n()) %>% as.data.frame


    df_data <- data.frame(
        id = 1:data_t$N,
        start_titre = data_t$initialTitreValue %>% as.data.frame %>% pull(!!biomarker),
        known_inf = data_t$knownInfsVec
    ) %>% mutate(known_inf = ifelse(known_inf == 0, "Not known", "Known"))

    df_data_post <- fit_states %>% filter(inf_ind == 1) %>% select(id, inf_ind) %>%
                summarise(inf_post = mean(inf_ind), n = n(), .by = "id") %>%
            left_join(df_data) %>% mutate(prop = n / n_length)

    cop_exp_sum_plot <- titre_cop_sum %>% left_join(df_data_post) %>% mutate(id = factor(id, levels = u_ids)) %>%
        filter(!is.na(n)) %>% mutate(biomarker = !!biomarker)

    }) %>% mutate(subtype = "h1")



fitfull_h3 <- readRDS(here::here("outputs", "fits", filename, paste0("fit_", "h3", ".RDS"))) 
outputfull_h3 <- readRDS(file = here::here("outputs", "fits", filename, paste0("pp_", "h3", ".RDS")))


n_post <- outputfull_h3$n_post
n_length <- n_chains * n_post
model_outline <- fitfull_h3$model
data_t <- fitfull_h3$data_t
fit_states <- outputfull_h3$fit_states
summarise_inf_prop_h3 <- map_df(1:length(model_outline$observationalModel), 
    function(i) {
    biomarker <- model_outline$observationalModel[[i]]$biomarker

    titre_cop_sum <- fit_states %>% filter(!!sym(biomarker) != -1) %>%
        group_by(id) %>% mean_qi(!!sym(biomarker) ) %>% arrange(!!biomarker) %>%
        rename(titre_val = !!(biomarker))
    u_ids <- titre_cop_sum$id

    df_data <- data.frame(
        id = 1:data_t$N,
        start_titre = data_t$initialTitreValue %>% as.data.frame %>% pull(!!biomarker),
        known_inf = data_t$knownInfsVec
    ) %>% mutate(known_inf = ifelse(known_inf == 0, "Not known", "Known"))

    df_data_post <- fit_states %>% filter(inf_ind == 1) %>% select(id, inf_ind) %>%
                summarise(inf_post = mean(inf_ind), n = n(), .by = "id") %>%
            left_join(df_data) %>% mutate(prop = n / n_length)

    cop_exp_sum_plot <- titre_cop_sum %>% left_join(df_data_post) %>% mutate(id = factor(id, levels = u_ids)) %>%
        filter(!is.na(n)) %>% mutate(biomarker = !!biomarker)

    }
) %>% mutate(subtype = "h3")



# summarise_inf_prop_h1

summarise_inf_prop_h1 %>% ggplot() + 
    geom_point(aes(x = titre_val, y = prop, alpha = prop)) + 
    facet_grid(vars(biomarker)) + 
    theme_bw()

pid_cols <- data_titre_h3 %>% select(pid, id) %>% unique
swab_info <- nih_inf_raw %>% filter(year == 2023) %>% left_join(pid_cols) %>% select(pid, id, swab_virus) %>% unique

df_post_compare <- bind_rows(
    swab_info %>% left_join(summarise_inf_prop_h1 %>% mutate(id = as.numeric(as.character(id))), by = "id" ) %>% mutate(subtype = "h1"),
    swab_info %>% left_join(summarise_inf_prop_h3 %>% mutate(id = as.numeric(as.character(id))), by = "id" ) %>% mutate(subtype = "h3")
) %>% mutate(prop = ifelse(is.na(prop), 0, prop))

summary_plot <- df_post_compare %>% select(pid, id, swab_virus, prop, subtype ) %>% unique %>% pivot_wider(names_from = "subtype", values_from = "prop") %>%
            mutate(h_diff = h1 - h3) 
            
summary_plot %>%
          ggplot() + 
            geom_point( 
            aes(x = h1, y = h3, color = swab_virus), size = 5, alpha = 0.7) + 
            theme_bw() 

summary_plot %>%
          ggplot() + 
            geom_point( 
            aes(x = pid, y = h_diff, color = swab_virus), size = 5, alpha = 0.7) + 
            theme_bw() + facet_wrap(vars(swab_virus)) + 
            labs(x = "PID", y = "Differece in probability of infection (h1 - h3)")