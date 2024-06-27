create_data <- function(df_titre, df_inf) {

  #df_titre %>% head
  #df_inf %>% head
    # in function
    df_titre_key <- df_titre %>% select(pid) %>% unique %>% mutate(id = 1:nrow(.))
    ids_titre_key_id <- df_titre_key %>% pull(pid)

    df_inf_clean <- df_inf %>% filter(pid %in% ids_titre_key_id )

    N <- nrow( df_titre %>% select(pid) %>% unique)
    ids_pcr <- df_inf_clean %>% filter(inf_type == "pcr") %>% left_join(df_titre_key) %>% pull(id)
    N_pcr <- nrow(df_inf_clean %>% filter(inf_type == "pcr"))
    ids_sero <- df_inf_clean %>% filter(inf_type == "sero") %>% left_join(df_titre_key) %>% pull(id)
    N_sero <- nrow(df_inf_clean %>% filter(inf_type == "sero"))

    if (sum(!df_titre %>% pivot_wider(names_from = "virus", values_from = "titre") %>% as.data.frame %>% complete.cases) > 0) {
      stop("ERROR! Some NAs in titre values")
    }
    x <- df_titre %>% pivot_wider(names_from = "virus", values_from = "titre") %>% select(3:ncol(.)) %>% as.matrix
    max_t <- apply(x, 2, max) %>% as.numeric

    y <- rep(0, N)
    y[unique(c(ids_pcr, ids_sero))] <- 1

    data_l <- list(
        obs_y = y,
        obs_x = x,
        biomarkers = colnames(x),
        N_known = N_pcr + N_sero,
        ids_sero = ids_sero,
        N_sero = N_sero,
        ids_pcr = ids_pcr,
        N_pcr = N_pcr,
        N_data = N,
        max_t = max_t
    )

}

create_model <- function(exp_prior, data_l) {
  model_null <- create_model_i(exp_prior, 0, data_l$max_t, data_l$biomarkers)
  model <- create_model_i(exp_prior, 0, data_l$max_t, data_l$biomarkers)

  model_null$evaluateLogLikelihood <- function(params, jump, datalist) {0}
  list(null = model_null, cop = model)
}

log_dirichlet_pdf <- function(x, alpha) {
  if (!is.vector(x) || !is.numeric(x) || !is.vector(alpha) || !is.numeric(alpha)) {
    stop("Input must be numeric vectors.")
  }
  
  k <- length(alpha)
  
  if (length(x) != k) {
    stop("Length of 'x' must match the length of 'alpha'.")
  }
  
 # if (any(x < 0) || sum(x) != 1) {
  #  stop("x must be a probability vector (summing to 1) with non-negative elements.", sum(x), " ", x < 0)
 # }
  
  if (any(alpha <= 0)) {
    stop("alpha must be a vector of positive values.")
  }
  
  log_beta <- sum(lgamma(alpha)) - lgamma(sum(alpha))
  
  log_pdf <- log_beta + sum((alpha - 1) * log(x))
  
  return(log_pdf)
}

create_model_i <- function(exp_prior, inf_prior, max_t, biomarkers) {
    
  c_value <- 0.33

  birth_prop <- function(N, N_known, k, c) {
     # pk <- log(1 / (N - N_known))#(k - N_known) * log(5) - 5 - lfactorial(k - N_known)#log(1 / (N - N_known)) # log(1 / choose(N - N_known, k - N_known))
     # pk_1 <- log(1 / (N - N_known))#(k - N_known) * log(5) - 5 - lfactorial(k - N_known)#log(1 / (N - N_known)) # log(1 / choose(N - N_known, k - N_known + 1))
      c * (1) #min(1, pk_1 / pk)
  }

  death_prop <- function(N, N_known, k, c) {
      #pk <- log(1 / (N - N_known))#(k - N_known) * log(5) - 5 - lfactorial(k - N_known)#log(1 / (N - N_known)) # log(1 / choose(N - N_known, k - N_known - 1))
      #pk_1 <- log(1 / (N - N_known))#(k - N_known) * log(5) - 5 - lfactorial(k - N_known)# log(1 / (N - N_known)) # log(1 / choose(N - N_known, k - N_known))
      c * (1) # min(1, pk / pk_1)
  }

  model <- list(

    exp_prior = exp_prior, 
    inf_prior = inf_prior,
    max_t = max_t,

    lowerParSupport_fitted = c(map(seq_len(length(biomarkers)), ~c(-max_t[.x], -1) ) %>% unlist),#, rep(0, length(biomarkers)) ),
    upperParSupport_fitted = c(map(seq_len(length(biomarkers)), ~c(max_t[.x], 1) ) %>% unlist),#, rep(10, length(biomarkers))) ,

    namesOfParameters = c(
      map(biomarkers, ~paste0(c("beta_0_", "beta_1_"), .x)) %>% unlist
   #   map(seq_len(length(biomarkers)), ~paste0(c("p_"), .x)) %>% unlist
      ),

    sampleInitPrior = function(datalist) {
     # cat("sampleInitPrior", "\n")

      logs <- c(map(seq_len(length(biomarkers)), ~ c(runif(1, -max_t[.x], max_t[.x]), runif(1, -1, 1)) ) %>% unlist)
     # p_simplex <- c(rexp(seq_len(length(biomarkers)), 1))
     # cat("lol: ", c(logs, p_simplex), "\n")
      c(logs)
    },

    sampleInitJump = function(params, datalist) {   
     # cat("sampleInitJump", "\n")

      start <- sample(which(datalist$obs_y == 0 ), 2)
      data_NIH_exp <- datalist$obs_y
      data_NIH_exp[start] <- 1
      data_NIH_exp <- c(data_NIH_exp, rep(1 / length(biomarkers), length(biomarkers)))
      data_NIH_inf <- c(datalist$obs_y, rep(0, length(biomarkers)))
      jump_new <- matrix(c(data_NIH_exp, data_NIH_inf), nrow = 2, byrow = TRUE)
      jump_new
    },

    
    evaluateLogPrior = function(params, jump, datalist) {
      #cat("evaluateLogPrior", "\n")
      N <- datalist$N_data

      p <- 0
      # Evlauate logistic values
      for (b in 1:length(biomarkers)) {
        beta0 <- params[2 * (b - 1) + 1]
        beta1 <- params[2 * (b - 1) + 2]
        p <- p + dunif(beta0, -max_t[b], max_t[b], log = TRUE)
        p <- p + dunif(beta1, -1, 1, log = TRUE)
      }
      p_simplex <- jump[1, N + seq_len(length(biomarkers))]

      p <- p + log_dirichlet_pdf(p_simplex, rep(0.5, length(biomarkers)))
      # Evaluate simplex values
     # for (b in 1:seq_len(length(biomarkers))) {
     #   p_virus <- params[2 * (length(biomarkers)) + b]
     #   p <- p + dexp(p_virus, 1, log = TRUE)
     # }
     # cat("in:  Lol", p)

      # cat(N, "\n")
      N_exp <- sum(jump[1, 1:N])
      N_inf <- sum(jump[2, 1:N])
      N_known <- datalist$N_known

      p_est <- (N_inf / N) * exp_prior / (1 - exp_prior) 

      p <- p + dnbinom(N_exp - N_known, N - N_known, 1-p_est, log = TRUE)  
      p <- p + log(1 / choose(N - N_known, N_exp - N_known))
      p <- p + log(1 / choose(N_exp - N_known, N_inf - N_known))
      #cat("p: ", p, "\n")

      p
    },

    evaluateLogLikelihood = function(params, jump, datalist) {
      #cat("in:  evaluateLogLikelihood", "\n")

      ll <- 0
      N <- datalist$N_data
      jump_e <- jump[1, 1:N]
      jump_i <- jump[2, 1:N]

      N_exp <- sum(jump[1, 1:N])
      N_inf <- sum(jump[2, 1:N])

      y <- datalist$obs_y
      x <- datalist$obs_x

      N_data <- datalist$N_data
      ll_mean <- 0
      j <- 0
    #  cat("in:  evaluateLogLikelihood", "\n")
      if (length(biomarkers) > 1) {
        p_simplex <- jump[1, N + seq_len(length(biomarkers))]
      } else {
        p_simplex <- 1
      }
      ll_b <- vector(mode = "numeric", length = length(biomarkers))
      p_est <- vector(mode = "numeric", length = length(biomarkers))

       ll <- 0

        for (i in 1:N_data) {
          if (jump_e[i] == 1) {
            for (b in 1:length(biomarkers)) {
                beta0 <- params[2 * (b - 1) + 1]
                beta1 <- params[2 * (b - 1) + 2]
                p_est[b] <- 1.0 / (1.0 + exp(- (beta0 - (beta0 / max_t[b] + beta1) *  as.numeric(x[i, b])) ) ) * p_simplex[b]
            }
            p <- sum(p_est)
            ll_i <- y[i] * log(p) + (1 - y[i]) * log(1 - p)
            ll <- ll + ll_i
            j <- j + 1
          } else {
          }
        }
        if (j == 0) {
          ll <- -1e10
        }
     # cat(ll_b, "\n")
      ll
    },

    sampleBirthProposal = function(params, jump, i_idx, datalist) {
     # cat("in:  sampleBirthProposal")
      N <- datalist$N_data
      jump_e <- jump[1, 1:N]
      jump_i <- jump[2, 1:N]

      N_exp <- sum(jump[1, 1:N])
      N_inf <- sum(jump[2, 1:N])

      id <- sample(1:N, 1)
      while(jump[1, id] == 1) {
        id <- sample(1:N, 1)
      }
      p_id <<- id

      jump[1, id] <- 1

      beta0 <- params[1]
      beta1 <- params[2]
      x <- datalist$obs_x
    #p <- 1.0 / (1.0 + exp(- (beta0 + beta1 *  x[p_id]) ) )
      #p <- 1.0 / (1.0 + exp(- (beta0 - (beta0 / max_t + beta1) * x[id]) ) )

      jump[2, id] <- rbinom(1, 1, inf_prior)

      jump

    },

    sampleDeathProposal = function(params, jump, i_idx, datalist) {
      # cat("in:  sampleDeathProposal")

      N <- datalist$N_data
      jump_e <- jump[1, 1:N]
      jump_i <- jump[2, 1:N]

      N_exp <- sum(jump[1, 1:N])
      N_inf <- sum(jump[2, 1:N])

      obs_y <- datalist$obs_y

      id <- sample(1:N, 1)
      while(jump[1, id] == 0 || obs_y[id] == 1) {
        id <- sample(1:N, 1)
      }
      d_id <<- id

      jump[1, id] <- 0
      jump[2, id] <- 0
      jump

    },

    evaluateBirthProposal = function(params, jump, i_idx, datalist) {
      #cat("evaluateBirthProposal")
      N <- datalist$N_data
      jump_e <- jump[1, 1:N]
      jump_i <- jump[2, 1:N]

      obs_y <- datalist$obs_y

      inf_i <- obs_y[p_id]

      N_exp <- sum(jump[1, 1:N])
      N_inf <- sum(jump[2, 1:N])

      N_known <- datalist$N_known

      beta0 <- params[1]
      beta1 <- params[2]

      x <- datalist$obs_x
      y <- datalist$obs_y

      p_birth <- birth_prop(N, N_known, N_exp - 1, c_value)
      p_death <- death_prop(N, N_known, N_exp, c_value)

      beta0 <- params[1]
      beta1 <- params[2]
      x <- datalist$obs_x

      log(N - N_exp + 1) - log(N_exp - N_known) - log(p_birth) + log(p_death) - (dbinom(inf_i, 1, inf_prior, log = TRUE))
    },

    evaluateDeathProposal = function(params, jump, i_idx, datalist) {
      #      cat("in:  evaluateDeathProposal")
      N <- datalist$N_data
      obs_y <- datalist$obs_y
      inf_i <- obs_y[d_id]

      jump_e <- jump[1, 1:N]
      jump_i <- jump[2, 1:N]

      N_exp <- sum(jump[1, 1:N])
      N_inf <- sum(jump[2, 1:N])

      beta0 <- params[1]
      beta1 <- params[2]

      x <- datalist$obs_x
      y <- datalist$obs_y

      N_known <- datalist$N_known

      p_birth <- birth_prop(N, N_known, N_exp, c_value)
      p_death <- death_prop(N, N_known, N_exp + 1, c_value)

      beta0 <- params[1]
      beta1 <- params[2]
      x <- datalist$obs_x

      log(N_exp - N_known + 1) - log(N - N_exp) + log(p_birth) - log(p_death) + (dbinom(inf_i, 1, inf_prior, log = TRUE))
    },

    sampleJump = function(params, jump, i_idx, datalist) {
      N <- datalist$N_data
      B <- length(biomarkers)
    #  B <- 4
      p_simplex <- jump[1, N + seq_len(B)]

    #  p_simplex <- c(0.999, 0.00033, 0.00033, 0.00033)
      p_idx <- sample(1:length(p_simplex), 1)
   #  p_idx <- 1
     # p_simplex <- c(0, 0, 0, 1)
     # p_idx <- 4
      q <- runif(1, 0, 1)
      if (q < 0.2) {
        p_simplex <- sample(p_simplex, B)
        jump[1, N + seq_len(B)] <- p_simplex
      } else {

        simp_update <- p_simplex[p_idx]
        sd <- runif(1, 0, 0.2)
        p_new <- min(max(simp_update[1] + rnorm(1, 0, sd), 0), 1)

        diff = (simp_update[1] - p_new) / (B - 1)
        p_simplex_new <- p_simplex + diff
        p_simplex_new[p_idx] <- p_new
        if (any(p_simplex_new < 0) || any(p_simplex_new > 1) ) {
          jump[1, N + seq_len(B)] <- p_simplex
        } else{
          p_simplex_new[p_idx] <- p_new
          jump[1, N + seq_len(B)] <- p_simplex_new
        }
      }

      jump
    },

    sampleProposal = function(params, jump, datalist) {
        #     cat("in:  sampleProposal")

      N <- datalist$N_data

      jump_e <- jump[1, 1:N]
      jump_i <- jump[2, 1:N]

      N_exp <- sum(jump[1, 1:N])
      N_inf <- sum(jump[2, 1:N])

      N_known <- datalist$N_known
    
      if (N_exp == N_known) {
        b_k <- birth_prop(N, N_known, N_exp, c_value)
        d_k <- 0
      } else if (N_exp == N) {
        b_k <- 0
        d_k <- death_prop(N, N_known, N_exp, c_value)
      } else {
        b_k <- birth_prop(N, N_known, N_exp, c_value)
        d_k <- death_prop(N, N_known, N_exp, c_value)
      }
      p_k <- 1 - b_k - d_k
      q <- c(d_k, d_k + p_k, d_k + p_k + b_k)

      q
    }
  )
  model
}


run_model <- function(model_exp, data_NIH, data_name) {

    model_exp_names <- names(model_exp)
 #   cat(str(model_exp))
  #    cat(model_exp_names)

    purrr::map2(model_exp[2], model_exp_names[2],
        function(model, names) {
        model_name <- paste0("base_", model$exp_prior)
        dir.create(here::here("outputs", "prior_pred", data_name))

        settings <-  list(
            numberChainRuns = 2,
            numberCores = 2,
            iterations = 80000,
            burninPosterior = 40000,
            thin = 100,
            consoleUpdates = 100,
            onAdaptiveCov = TRUE,
            updatesAdaptiveCov = 100,
            burninAdaptiveCov = 1000,
            updatesAdaptiveTemp = 10,
            covarInitVal = 1e-6, # make very small if struggling to sample to beginning
            covarInitValAdapt = 1e-2, # make very small if struggling to sample to beginning
            covarMaxVal = 1, # decrease if struggling to sample in the middle
            runParallel = TRUE,
            onDebug = FALSE,
            noGibbsSteps = 1,
            numberFittedPar = 2 * length(data_NIH$biomarkers),
            lowerParBounds = model$lowerParSupport_fitted,
            upperParBounds = model$upperParSupport_fitted
        )

        outputs <- rjmc_func_pp(model, data_NIH, settings)
        saveRDS(outputs, here::here("outputs", "prior_pred", data_name, paste0(model_name, "_", names, ".RDS")))
        }
    )
}

plot_outputs <- function(outputs, model, data_X, filename, typename, n_chains, plot_helper) {
 # outputs, data_unif_A, "fixed_I_unif_A", "null"
 #outputs
# model <- model_exp[[2]]
 # n_chains <- 2
 #data_X <- data_NIH_2023_H1

  S <- outputs$mcmc[[1]] %>% dim %>% .[1]

  #paste0(data_name, "/", model_name), "cop", 2, plot_helper

  dir.create(here::here("outputs", "prior_pred", filename))
  N <- data_X$N_data
  max_t <- model$max_t
  logisticoutput <- map_df(1:n_chains, ~as.data.frame(outputs$mcmc[[.x]]) %>% mutate(chain = .x)) 
 
  post_summary <- map_df(1:length(data_X$biomarkers),
    function(b) {
      1:nrow(logisticoutput) %>% 
        map_df(
          ~data.frame(
            biomarker = data_X$biomarkers[b],
            chain = logisticoutput[.x, length(data_X$biomarkers) * 2 + 1],
            t = seq(0, max_t[b], length = 50),
            p = 1.0 / (1.0 + exp(- (logisticoutput[.x, 2*(b-1) + 1] - (logisticoutput[.x, 2*(b-1) + 1] / max_t[b] + logisticoutput[.x, 2*(b-1) + 2]) *  seq(0, max_t[b], length = 50)) ) )
          )
        ) %>% group_by(t, biomarker) %>% mean_qi(p)
    }
  ) %>% mutate(biomarker = factor(biomarker, levels =  names(plot_helper$virus_color)))

  n_chain <- n_chains
  s <<- 1
  summary_post <- map_df(1:n_chain,
    function(c) {
      list_outputs_mode <- list()
      j <- 1
      for (i in 1:length(outputs$jump[[c]])) {
        sample_E <- sum(outputs$jump[[c]][[i]][1, ])
        sample_I <- sum(outputs$jump[[c]][[i]][2, ])

        list_outputs_mode[[j]] <- data.frame(
            E = sample_E,
            I = sample_I,
            sample = s,
            chain = c
        )
        s <<- s + 1
        j <- j + 1
      }
      list_outputs_mode %>% bind_rows()
    }
  ) 

  post_pcr <- data_X$obs_x[data_X$ids_pcr, ]
  post_sero <- data_X$obs_x[data_X$ids_sero, ]

  p1 <- data.frame(
    value = c(data_X$N_pcr, data_X$N_sero, data_X$N_known * model$exp_prior / (1 -  model$exp_prior), data_X$N_pcr, data_X$N_sero, summary_post %>% summarise(mean(E)) %>% as.numeric),
    param = c("I_pcr", "I_sero", "E", "I_pcr", "I_sero", "E"),
    bayes_type = c("Prior", "Prior", "Prior", "Posterior", "Posterior", "Posterior")
  ) %>%
    mutate(name_plot = case_when(
      param == "E" ~"Exposed and not infected",
      param == "I_pcr" ~"PCR confirmed",
      param == "I_sero" ~"Serologically inferred",
      )) %>%
    ggplot() + 
      geom_col(aes(y = bayes_type, x = value, color = name_plot, fill = name_plot), alpha = 0.6) + 
      scale_fill_manual(values = c( "#fdeb8b","#938a59", "#ce5837")) + 
      scale_color_manual(values = c("#fdeb8b", "#938a59", "#ce5837")) + 
      scale_x_continuous(expand = c(0,0)) + 
      scale_y_discrete(expand = c(0,0)) + 
      guides(color = "none") + 
        theme(
    plot.title = element_text(size = 16, hjust = 0.5),  # Title style
    axis.title.x = element_text(size = 14),  # X-axis label style
    axis.title.y = element_text(size = 14),  # Y-axis label style
    axis.text.x = element_text(size = 12),  # X-axis tick label style
    axis.text.y = element_text(size = 12),  # Y-axis tick label style
   # panel.grid.major = element_blank(),  # No major gridlines
    panel.grid.minor = element_blank(),  # No minor gridlines
    panel.background = element_blank(),  # No panel background
    axis.line.x = element_line(color = "black"),  # Color of x-axis line
    axis.line.y = element_blank(),  # Color of y-axis line    legend.position = "none"  # No legend
     plot.margin = unit(c(0, 0, 0, 0), "cm") ,
     aspect.ratio = 0.1,
     legend.position = "top"
  ) +
  labs(x = paste0("Number of individuals (of ", data_X$N_data, ")"), y = "", fill = "") + 
  ggtitle(paste0("Assuming VE is ", model$exp_prior) )

  df_point <- map_df(1:length(data_X$biomarkers),
    function(b) {
      data.frame(
        biomarker =  data_X$biomarkers[b],
        value = c(post_pcr[, b], post_sero[, b]),
        inf_prop = rep(1, length(c(post_pcr[, b], post_sero[, b]))),
        type = c(rep("PCR confirmed", length(post_pcr[,b])), rep("Serologically inferred", length(post_sero[,b])))
      )
    }
  ) %>%
    mutate(biomarker = factor(biomarker, levels =  names(plot_helper$virus_color)))

  df_post <- 
  map_df(1:length(data_X$biomarkers),
  function(b) {
      map_df(1:n_chains, 
        ~data.frame(
          id = 1:data_X$N_data,
          biomarker = data_X$biomarkers[b],
          x = data_X$obs_x[, b],
          chain = .x,
          exp_mean = outputs$jump[[.x]] %>% unlist %>%
                matrix(nrow = S, byrow = TRUE) %>% apply(2, mean) %>% .[seq(1, N * 2, 2)],
          inf_mean = outputs$jump[[.x]] %>% unlist %>%
                matrix(nrow = S, byrow = TRUE) %>% apply(2, mean) %>% .[seq(2, N * 2, 2)]
        ) 
      ) 
    }
  )

  df_post_exp <- df_post %>% group_by(id, biomarker, x) %>% summarise(exp_mean = mean(exp_mean), inf_mean = mean(inf_mean)) %>%
    filter(!id %in% c(data_X$ids_pcr, data_X$ids_sero)) %>% ungroup %>% group_by(x, biomarker) %>%
    summarise(exp_mean = mean(exp_mean), inf_mean = mean(inf_mean), n = n()) %>%
    mutate(biomarker = factor(biomarker, levels =  names(plot_helper$virus_color)))

  p2 <- post_summary %>% 
    ggplot() + 
    geom_hline(yintercept = 0.5, color = "gray40" ) + 
    geom_count(data = df_point, aes(value, y = inf_prop, fill = type), shape = 21, position = position_dodge(0.5) ) + 
    geom_point(data = df_post_exp, aes(x, y = inf_mean, alpha = exp_mean, size = n), fill = "#fdeb8b", shape = 21 ) + 
    geom_ribbon(aes(x = t, ymin = .lower, ymax = .upper, fill = biomarker),  alpha = 0.4) + 
      geom_line(aes(x = t, y = p, color = biomarker), size = 2) + 
     guides(size = "none", fill = "none", alpha = "none", color = "none") + 
      theme_bw() + theme_minimal() + 
    scale_fill_manual(values = 
      c("PCR confirmed" = "#938a59", "Serologically inferred" = "#ce5837",
        plot_helper$virus_color)) +
    scale_color_manual(values = 
     plot_helper$virus_color)  +
        theme_minimal() +  # Minimalistic theme
  scale_x_continuous(breaks = plot_helper$titre_convert,
    labels = names(plot_helper$titre_convert), expand = c(0, 0))  + 
    scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),  # Title style
    axis.title.x = element_text(size = 14),  # X-axis label style
    axis.title.y = element_text(size = 14),  # Y-axis label style
    axis.text.x = element_text(size = 12),  # X-axis tick label style
    axis.text.y = element_text(size = 12),  # Y-axis tick label style
   # panel.grid.major = element_blank(),  # No major gridlines
    panel.grid.minor = element_blank(),  # No minor gridlines
    panel.background = element_blank(),  # No panel background
    axis.line.x = element_line(color = "black"),  # Color of x-axis line
    axis.line.y = element_blank(),  # Color of y-axis line    legend.position = "none"  # No legend
     plot.margin = unit(c(1, 1, 1, 1), "cm") 
  ) +
  labs(y = "Probability of \ninfection given exposure",
    x = "HAI titre at exposure") + 
    facet_wrap(vars(biomarker))

## simplex
    simplex_df <-  map_df(1:n_chains, 
        ~data.frame(
          biomarker = data_X$biomarkers,
          chain = .x,
          post_simplex = outputs$jump[[.x]] %>% unlist %>%
                matrix(nrow = S, byrow = TRUE) %>% apply(2, mean) %>% .[seq(N * 2 + 1, 2*N + length(data_X$biomarkers)*2, 2)]
        ) 
      ) 
    

  p3 <- simplex_df %>% mutate(biomarker = factor(biomarker, levels =  names(plot_helper$virus_color))) %>% 
    ggplot() + 
      geom_col(aes(y = chain, x = post_simplex, color = biomarker, fill = biomarker), alpha = 0.6) + 
      scale_fill_manual(values = plot_helper$virus_color) + 
      scale_color_manual(values = plot_helper$virus_color) + 
      scale_x_continuous(expand = c(0,0)) + 
      scale_y_discrete(expand = c(0,0)) + 
      guides(color = "none") + 
      theme(
    plot.title = element_text(size = 16, hjust = 0.5),  # Title style
    axis.title.x = element_text(size = 14),  # X-axis label style
    axis.title.y = element_text(size = 14),  # Y-axis label style
    axis.text.x = element_text(size = 12),  # X-axis tick label style
    axis.text.y = element_text(size = 12),  # Y-axis tick label style
   # panel.grid.major = element_blank(),  # No major gridlines
    panel.grid.minor = element_blank(),  # No minor gridlines
    panel.background = element_blank(),  # No panel background
    axis.line.x = element_line(color = "black"),  # Color of x-axis line
    axis.line.y = element_blank(),  # Color of y-axis line    legend.position = "none"  # No legend
     plot.margin = unit(c(0, 0, 0, 0), "cm") ,
     aspect.ratio = 0.1,
     legend.position = "bottom"
  ) +
  labs(x = paste0("Contribution of each measured \nbiomarker to likelihood"), y = "", fill = "") 

  p1 / p2 / p3 + plot_layout(heights = c(1, 5, 1))

  ggsave(here::here("outputs", "prior_pred", filename, paste("cop_", typename, ".pdf")), width = 10, height = 10)

  require(patchwork)


  summary_post %>% group_by(E, I) %>% summarise(n = n()/ nrow(.)) %>%
    ggplot() + 
        geom_hline(yintercept =  data_X$N_known, color = "gray") + 
        geom_vline(xintercept =  data_X$N_known, color = "gray") + 
        geom_abline(slope = 1, color = "gray90", linetype = "dashed") + 
        geom_tile(aes(x = I, y = E, fill = n)) + 
        geom_histogram(aes(x = I), bins = 30, fill = "transparent", color = "black") +  # Histogram for x-axis
        geom_histogram(aes(y = E), bins = 30, fill = "transparent", color = "black") +  # Histogram for y-axis
        theme_minimal() +
        theme(
            panel.grid = element_blank()
        ) + 
        labs(x = "Number of infected individuals", y = "Number of exposed individuals")
    ggsave(here::here("outputs", "prior_pred", filename, paste("hist_EI_", typename, ".pdf")))

  p1 <- summary_post %>% group_by(E, chain) %>% summarise(n = n()) %>%
    ggplot() + geom_col(aes(x = E, y = n, fill = chain)) + 
      geom_smooth(data = summary_post %>% group_by(E) %>%summarise(n = n()), aes(x = E, y = n))
  p2 <- p1 + 
    facet_grid(vars(chain))

  p3 <- summary_post %>% group_by(I, chain) %>% summarise(n = n()) %>%
    ggplot() + geom_col(aes(x = I, y = n, fill = chain)) + 
      geom_smooth(data = summary_post %>% group_by(I) %>% summarise(n = n()), aes(x = I, y = n))
  p4 <- p3 + 
    facet_grid(vars(chain))
  (p1 / p2) | (p3 / p4)
  ggsave(here::here("outputs", "prior_pred", filename, paste("E_dist_", typename, ".pdf")))
}





#  logisticoutput_i <- logisticoutput %>% filter(chain == 1)
  
#  full_sample_I <- map_df(1:N,
#    function(j) {
#      COP_i <- (1.0 / (1.0 + exp(- (logisticoutput_i[, 1] - (logisticoutput_i[, 1] / max_t + logisticoutput_i[, 2]) * data_X$obs_x[j]) ) )) 
#      E <- outputs$jump[[1]] %>% unlist %>%
#            matrix(nrow = 200, byrow = TRUE) %>% .[, 2*(j - 1) + 1] 
#      I <- outputs$jump[[1]] %>% unlist %>%
#            matrix(nrow = 200, byrow = TRUE) %>% .[, 2*(j - 1) + 2] 
#      data.frame(
#        s = 1:200,
#        i = j,
#        E = E,
#        I = I,
#        I_inf = E * COP_i
#      )
#    }
#  )

  #df_sum_I <- full_sample_I %>% group_by(s) %>% summarise(I_sum = sum(I), I_inf_sum = sum(I_inf))
  #p0A <- df_sum_I %>% 
  #  ggplot() + geom_histogram(aes(I_sum))+ theme_bw() + labs(x = "Number of infections (bin)")
  #p0B <- df_sum_I %>% 
  #  ggplot() + geom_histogram(aes(I_inf_sum))+ theme_bw() + labs(x = "Number of infections (cop)")

  ##post_E <- outputs$jump[[1]] %>% unlist %>%
  #      matrix(nrow = 200, byrow = TRUE) %>% apply(2, mean) %>% .[seq(1, 2*N, 2)]
  #post_I <- outputs$jump[[1]] %>% unlist %>%
  #      matrix(nrow = 200, byrow = TRUE) %>% apply(2, mean)  %>% .[seq(2, 2*N, 2)]

  #logisticoutput_i <- logisticoutput %>% filter(chain == 1)
  #post_COP <- map(1:N, ~(1.0 / (1.0 + exp(- (logisticoutput_i[, 1] - (logisticoutput_i[, 1] / max_t + logisticoutput_i[, 2]) * data_X$obs_x[.x]) ) )) %>% mean) %>%
  #  unlist

  #df_post_alt <- data.frame(
  #  chain = 1,
  #  x = data_X$obs_x,
  #  post_COP = post_COP,
  #  post_E = post_E,
 # #  post_I = (post_COP * post_E) 
  #)
  
 # p1 <- df_post_alt %>% 
  #  ggplot() + geom_point(aes(x = x, y = post_E, fill = as.character(chain)), shape = 21, alpha = 0.7) + theme_bw() + 
  #  labs(x = "Pre-study titre", y = "Proportion of posterior exposed")
  #p2 <- df_post_alt %>% 
  #  ggplot() + geom_point(aes(x = x, y = post_I, fill = as.character(chain)), shape = 21, alpha = 0.7) + theme_bw() + 
  #  labs(x = "Pre-study titre", y = "Proportion of posterior infected")
  #p1 / p2 / (p0A + p0B)
  #ggsave(here::here("outputs", "prior_pred", filename, paste("post_dist_", typename, ".pdf")))


#}