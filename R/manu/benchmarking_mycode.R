

# Get TRANSVIR data run timings

cop_label <-c(rep("With COP", 11), rep("No COP", 11))
uncert_label <- rep(c(0.01, seq(0.05, 0.5, 0.05)), 2)

df_sim_status <- map_df( 1:22, function(i) {
        readRDS(here::here("hpc", "sim", paste0("bm_t_", i, ".RData")) ) %>%
            mutate(labels_cop = cop_label[i]) %>%
            mutate(uncert_label = uncert_label[i])
    }
) %>% mutate(total_time = as.numeric(total_time))

p1 <- df_sim_status %>% 
    ggplot() + 
        geom_point(aes(x = uncert_label, y = total_time / 60, color = labels_cop), size = 3, alpha = 0.9) + 
        theme_bw() + 
        labs(x = "Uncertainty in the data", y = "Run time (minutes)", color = "Model", title = "Total run time") + theme_bw()

p2 <- df_sim_status %>% 
    ggplot() + 
        geom_point(aes(x = uncert_label, y = total_time / 400000, color = labels_cop), size = 3, alpha = 0.9) + 
        theme_bw() + 
        labs(x = "Uncertainty in the data", y = "Run time (seconds)", color = "Model", title = "Run time per MCMC iteration (seconds)") + theme_bw()
p1 / p2 +
    plot_annotation(tag_levels = 'A') + plot_layout(guides = "collect") & 
    theme(plot.tag = element_text(size = 12, face = "bold")) & 
    scale_x_continuous(breaks = seq(0, 0.5, 0.05)) 
ggsave(here::here("outputs", "figs", "supp", "bm_sim.pdf"), width = 10, height = 5)

# Get TRANSVIR data run timings
transvir1 <- readRDS(here::here("hpc", "transvir", "bm_t_1.RData"))
transvir2 <- readRDS(here::here("hpc", "transvir", "bm_t_2.RData"))

bind_rows(
    transvir1 %>% mutate(labs = "TRANSVIR, NO PCR"), 
    transvir2 %>% mutate(labs = "TRANSVIR, PCR")
) %>% mutate(total_time = as.numeric(total_time)) %>% 
    mutate(total_time_mins = total_time / 60) %>%
    mutate(total_time_mc_secs = total_time / 400000) 



## Probably do a run with more data than other

files <- list.files("hpc/sim_comp", pattern = "\\.RData$", full.names = TRUE)
data_list <- list()

df_full_times <- map_df (files, 
    function(f) {
      #  cat(f, "\n")
        matches <- regmatches(f, regexec("bm_t_(\\d+)_(\\d+)\\.RData", f))[[1]]
       # cat(as.integer(matches[2]), "\n")
     #   cat(as.integer(matches[3]), "\n")
        data_perm <- readRDS(f)
      #  print(data_perm)
        data.frame(
            model_num = as.integer(matches[2]),
            model_size = c(10, 20, 50, 100, 200, 500, 1000)[(as.integer(matches[2]) - 1) %% 7 + 1],
            model_type = c("With COP", "No COP")[as.integer(matches[2]) %/% 7 + 1],
            sample = as.integer(matches[3]),
            time = as.numeric(readRDS(f)$total_time),
            time_ind = as.numeric(readRDS(f)$total_time) / (c(10, 20, 50, 100, 200, 500, 1000)[(as.integer(matches[2]) - 1) %% 7 + 1])
        )
    }
)

p1 <- df_full_times %>% filter(!model_num %in% c(7, 14)) %>% 
    ggplot() + 
        geom_boxplot(aes(x = factor(model_size), y = time / 60, fill = factor(model_type)),  alpha = 0.5)  + 
        labs(x = "Number of individuals in seroepidemiological study", y = "Run time (minutes)", fill = "Model", title = "Total run time") + theme_bw()
p2 <- df_full_times %>% filter(!model_num %in% c(7, 14)) %>% 
    ggplot() + 
        geom_boxplot(aes(x = factor(model_size), y = time / 400000, fill = factor(model_type)),  alpha = 0.5)  + 
        labs(x = "Number of individuals in seroepidemiological study", y = "Run time (seconds)", fill = "Model", title = "Run time per MCMC iteration (seconds)") + theme_bw()
p3 <- df_full_times %>% filter(!model_num %in% c(7, 14)) %>% 
    ggplot() + 
        geom_boxplot(aes(x = factor(model_size), y = time_ind, fill = factor(model_type)),  alpha = 0.5)  + 
        labs(x = "Number of individuals in seroepidemiological study", y = "Run time (seconds)", fill = "Model", title = "Run time per MCMC iteration (seconds)") + theme_bw()
p1 / p2 / p3
ggsave(here::here("outputs", "figs", "supp", "bm_full_times.pdf"), width = 10, height = 5)
