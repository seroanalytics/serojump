library(testthat)
library(data.table)
library(dplyr)

modelhier_A <- readRDS(here::here("outputs", "fits", "fudan_e3", "base_hier", "model_summary.rds"))


model_summary <- modelhier_A

# Sample test data
test_df <- data.table(
  id = c(1, 1, 2, 2),
  sample = c(1, 2, 1, 2),
  biomarker = c("X", "Y", "X", "Y"),
  time = c(1, 2, 1, 2)
)

test_that("compute_individual_trajectories returns expected output", {
  model_summary <- modelhier_A

  params <- initialize_parameters(model_summary)
  bio_markers <- params$model_outline$infoModel$biomarkers
  sample_ids <- sample(1:params$M, 10)
  param_samples <- prepare_posterior_samples(params$post, params$chain_samples)
  
  exposure_order <- get_exposure_order(params$model_outline, params$N)
  types_traj <- unique(exposure_order$exp_type)

  traj_type <- types_traj[1]
  
  selected_ids <- exposure_order %>% filter(exp_type == traj_type) %>% pull(id) %>% unique()
  selected_ids <- sample(selected_ids, min(20, length(selected_ids)))

  traj_data <- compute_individual_trajectories(
    id = selected_ids,
    bio_markers = bio_markers,
    samples = sample_ids,
    exposure_order = exposure_order,
    model = params$model_outline,
    param_samples = param_samples,
    T_max = params$T_max
  )
    
  
  expect_true(is.data.frame(output))
  expect_true(all(output$id %in% test_ids))
  expect_true(all(output$biomarker %in% test_bio_markers))
  expect_true(all(output$t <= test_T_max))
})


test_that("compute_individual_trajectories returns expected output", {

  model_summary <- modelhier_A

  params <- initialize_parameters(model_summary, 10)
    # Extract the biomarkers, the sampled numbers and the resulting posterior values 
    bio_markers <- params$model_outline$infoModel$biomarkers
    sample_ids <- sample(1:params$M, params$S)
    param_samples <- prepare_posterior_samples(params$post, params$chain_samples)
    
    # Get the exposure order for each indivdual and each chain in the posterior
    df_inferred_exp <- compute_exposure_order(params, bio_markers, sample_ids)

    # Get the original serologicla data
    data_fit_list <- extract_serological_data(params$data_t, bio_markers, params$N)

    # Get the trajectory types
    types_traj <- unique(exposure_order$exp_type)
  traj_type <- types_traj[1]
  
  selected_ids <- df_inferred_exp %>% filter(exp_type == traj_type) %>% pull(id) %>% unique()
  selected_ids <- sample(selected_ids, min(20, length(selected_ids)))

  model <- params$model_outline

  i <- selected_ids[1]
  bio <- bio_markers[1]
  s <- sample_ids[1]

  exp_data <- filter_sorted_data(df_inferred_exp, i, s, bio) %>% mutate(row_id = row_number())
  
  times <- c(exp_data$time, params$T_max)
  time_diff <- diff(times)
  titre_traj <- c()
  titre_anchor <- NULL

  titre_start <- data_fit_list[[i]] %>% filter(row_num == 1, bio == bio) %>% pull(titre)

  for (row_idx in seq_len(nrow(exp_data))) {
    exp_type <- exp_data[row_idx, ]$exp_type
    id_key <- which(model$abkineticsModel %>% map(~.x$exposureType ) %>% unlist() == exp_data[row_idx, ]$exp_type)
    model_func <- model$abkineticsModel[[id_key]]$funcForm
    base_params <- model$abkineticsModel[[id_key]]$pars
    
    group_idx <- extract_hierarchical_param(model, id_key, i)
    param_list <- calculate_param_list(base_params, model$abkineticsModel[[id_key]]$hierFlag, model, id_key)
    model_params <- extract_parameters(model, id_key, param_samples, param_list, model_summary, s, group_idx) %>% as.numeric
    
    new_titre <- compute_trajectory(if (is.null(titre_anchor)) titre_start else titre_anchor, time_diff[row_idx], model_func, model_params)
    titre_traj <- c(titre_traj, new_titre)
    titre_anchor <- last(new_titre)
  }

  # map_dbl(seq_len(time_vector), ~ model_function(anchor_titre, .x, params))
  
  data.table(
    id = i,
    sample = s,
    t = 1:T_max,
    biomarker = bio,
    titre_trajectory = titre_traj
  ) %>% filter(!is.na(titre_trajectory))

 
  }
)
