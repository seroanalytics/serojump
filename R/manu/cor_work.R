# 1. Load the relevant data 

devtools::load_all()
model_summary_w2 <- readRDS(file = "./outputs/fits/transvir_data/wave2_base/model_summary.RDS")

save_info_w2 <- list(
    file_name = "transvir_data",
    model_name = "wave2_no_pcr"
)

### COP work 

model_summary <- model_summary_w2

fitfull <- model_summary$fit    
outputfull <- model_summary$post

fit_states <- outputfull$fit_states
model_outline <- fitfull$model

posterior_binary <- fit_states %>% mutate(titre_type = case_when(
    inf_time == -1 ~ "No inf",
    inf_time != -1 ~ "Inf")
) %>% group_by(id) %>% mutate(sample = row_number()) %>% ungroup %>% 
mutate(
    spike_s = (spike - min(spike))/ (max(spike) - min(spike)), 
    NCP_s = (NCP - min(NCP))/ (max(NCP) - min(NCP)), 
)

posterior_binary_sum <- posterior_binary %>% group_by(id) %>% summarise(
    inf_ind = mean(inf_ind), spike_s = mean(spike_s), NCP_s = mean(NCP_s)
)



P <- 20

summary_discrete <- map_df(1:length(model_outline$observationalModel), 
    function(i) {
    biomarker <- model_outline$observationalModel[[i]]$biomarker

    posterior_binary %>% mutate(titre_cut = cut(!!sym(biomarker), breaks = P)) %>% 
        group_by(titre_cut) %>% summarise(inf = sum(inf_ind), n = n()) %>% mutate(inf_risk = inf/n) %>% 
        complete(titre_cut, fill = list(inf = 0, n = 0, inf_risk = 0)) %>%
        mutate(biomarker = biomarker, decile = 1:P)  %>% mutate(dens = n / sum(n))
        }
    )
    


# FT-Style Theme
theme_ft <- function() {
  theme_minimal(base_family = "Arial") +  # Simple sans-serif font
    theme(
      panel.grid.major.y = element_line(color = "grey80", size = 0.3),  # Light grid lines
      panel.grid.major.x = element_blank(),  # No vertical grid lines
      panel.grid.minor = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 14, color = "grey40", margin = margin(b = 15)),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      plot.margin = margin(20, 20, 20, 20)
    )
}


#### NOTES 
##### posterior_binary is the full posterior binaries and critical if you wanna fit a logistic model 
##### summary_discrete is the summary of the discrtise data and ctitical if you wanna fit a continuous curve 

# 2. Fitting to infection risk (or clinicial risk) per titre group 

## 2.1 For continuous titre vaues 

posterior_binary
posterior_binary_sum
summary_discrete



## use stan model
require(cmdstanr)
data_list_A <- list(
    N = nrow(posterior_binary),
    x = posterior_binary$spike_s,
    y = posterior_binary$inf_ind
)

data_list_B <- list(
    N = nrow(posterior_binary_sum),
    x = posterior_binary_sum$spike_s,
    y = posterior_binary_sum$inf_ind
)

log_curve_fit <- cmdstan_model( here::here("src", "stan", "cops", "logit_curve_fit.stan") )
post_A <- log_curve_fit$sample(data = data_list_A, parallel_chains = 4, cores = 4, iter_warmup = 1000, iter_sampling = 1000)

y_hat_new_post <- post_A$draws() %>% gather_draws(y_hat_new[i]) 
post_under <- post_A$draws() %>% spread_draws(L, k, x0) 

df_xy <- data.frame(
    x = data_list_B$x,
    y = data_list_B$y
)

y_hat_new_post %>% 
    ggplot() + 
        geom_point(aes(x = x * 100, y = y), data = df_xy, color = "black") +
        stat_lineribbon(aes(x = i, y = .value), .width = 0.95, fill = "red", alpha = 0.3) + 
        theme_bw()


## use stan model
require(cmdstanr)
summary_discrete %>% filter(biomarker == "spike") -> infection_risk_bio_n
data_list_C <- list(
    N = nrow(infection_risk_bio_n),
    x = infection_risk_bio_n$decile,
    y = infection_risk_bio_n$inf_risk,
    weights = infection_risk_bio_n$n
)
log_curve_fit_w <- cmdstan_model( here::here("src", "stan", "cops", "logit_curve_fit_weighted.stan") )
post_C <- log_curve_fit_w$sample(data = data_list_C, parallel_chains = 4, iter_warmup = 1000, iter_sampling = 1000)

y_hat_new_post <- post_C$draws() %>% gather_draws(y_hat_new[i]) 
y_hat_post <- post_C$draws() %>% gather_draws(y_hat[i]) 

post_under <- post_C$draws() %>% spread_draws(L, k, x0) 


df_xy <- data.frame(
    w = data_list_C$weights,
    x = data_list_C$x,
    y = data_list_C$y
)

p1 <- y_hat_post %>% 
    ggplot() + 
        geom_point(aes(x = x, y = y, size = w), data = df_xy, color = "black") +
        stat_lineribbon(aes(x = i, y = .value), .width = 0.95, fill = "red", alpha = 0.3) + 
        theme_bw()

p2 <- y_hat_new_post %>% 
    ggplot() + 
        stat_lineribbon(aes(x = i, y = .value), .width = 0.95, fill = "red", alpha = 0.3) + 
        theme_bw()
p1 + p2