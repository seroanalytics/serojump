filename <- "hpc/transvir_w2_inf/p3"
n_chains <- 4
require(readxl)
require(lubridate)


gambia_pvnt_w2 <- get_data_titre_model_wave2_pvnt_iga()# This is the empirical prior for the exposure time
gambia_exp_w2 <- get_exposures_wave2() # This is the empirical prior for the exposure time
exp_prior_w2 <- get_exp_prior_wave2() # This is the empirical prior for the exposure time

fitfull_w2 <- readRDS(here::here("outputs", "fits", filename, paste0("fit_", "w2", ".RDS")))
outputfull_w2 <- readRDS(file = here::here("outputs", "fits", filename, paste0("pp_", "w2", ".RDS")))


ids_infer <- outputfull_w2$fit_states %>% group_by(id) %>% summarise(prop = sum(inf_ind) / 800 ) %>% filter(prop > 0.75) %>% pull(id)

df_temp <- outputfull_w2$fit_states %>% filter(id %in% ids_infer) %>% group_by(id) %>% filter(inf_ind == 1) %>% 
    summarise(sVNT = mean(sVNT), IgA = mean(IgA))

df_temp


# Get priors for peak value
joiner <- gambia_exp_w2 %>% filter(exposure_type == "delta") %>% rename(inf_time = time) %>% select(!exposure_type)
summarise_post_inf <- gambia_pvnt_w2 %>% left_join(joiner) %>% filter(!is.na(inf_time) ) %>% mutate(time_diff = time - inf_time) %>% 
    filter(time_diff > 0) 

summarise_post_inf %>%
        ggplot() + 
            geom_point(aes(time_diff, sVNT)) + 
            geom_smooth(aes(time_diff, sVNT))

summarise_post_inf %>%
        ggplot() + 
            geom_point(aes(time_diff, IgA)) + 
            geom_smooth(aes(time_diff, IgA))
