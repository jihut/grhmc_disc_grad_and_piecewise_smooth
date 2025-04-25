source("piecewise_smooth/implementation_scripts/VARI/VARI_adaptive_step_grhmc_piecewise_smooth_density_transformed_function.R")
source("piecewise_smooth/implementation_scripts/general_scripts/grhmc_piecewise_smooth_density_transformed_function.R")

VARI_grhmc_piecewise_smooth_density_transformed_function <- function(
    model_list, 
    lambda_initial, 
    lambda_lower_limit = 0.001, 
    time_period_adaptive = 5000,
    time_period_generating_samples = 5000,
    n_generated_samples = 5000,
    n_adaptive_samples = 1000,
    diag_s_elements_initial = NULL,
    m_initial = NULL,
    qbar_initial = NULL,
    pbar_initial = NULL,
    random_state = NULL,
    rtol = NULL,
    atol = NULL,
    h_max = 1.0,
    proportion_time_until_adaptive_start = 0.05,
    min_proportion_of_previous_state = 0.5,
    max_proportion_of_previous_state = 2,
    sampling_compute_temporal_averages_of_moments = FALSE, 
    reflection_type = NULL,
    verbose_at_refresh = FALSE,
    last.root.offset.lin.root.finder = 1e-10,
    last.root.offset.non.lin.root.finder = 1e-10,
    precision_real_root_lin_root_finder = 1e-13,
    num_subdiv_non_lin_root_finder = 8L
) {
  
  if(!is.null(random_state)){
    set.seed(random_state)
  }
  
  print("Running adaptive step")
  
  time_period_adaptive <<- time_period_adaptive
  
  adaptive_step_run <- VARI_adaptive_step_grhmc_piecewise_smooth_density_transformed_function(
    model_list = model_list, 
    lambda = lambda_initial,
    lambda_lower_limit = lambda_lower_limit,
    T = time_period_adaptive,
    n_samples = n_adaptive_samples,
    m_initial = m_initial,
    diag_s_elements_initial = diag_s_elements_initial,
    qbar_initial = qbar_initial,
    pbar_initial = pbar_initial,
    rtol = rtol,
    atol = atol,
    h_max = h_max,
    proportion_time_until_adaptive_start = proportion_time_until_adaptive_start,
    min_proportion_of_previous_state = min_proportion_of_previous_state,
    max_proportion_of_previous_state = max_proportion_of_previous_state,
    reflection_type = reflection_type,
    verbose_at_refresh = verbose_at_refresh,
    last.root.offset.lin.root.finder = last.root.offset.lin.root.finder,
    last.root.offset.non.lin.root.finder = last.root.offset.non.lin.root.finder,
    precision_real_root_lin_root_finder = precision_real_root_lin_root_finder,
    num_subdiv_non_lin_root_finder = num_subdiv_non_lin_root_finder
  )
  
  print("Running sample step")
  
  sample_step_run <- grhmc_piecewise_smooth_density_transformed_function(
    model_list = model_list,
    lambda = adaptive_step_run$lambda_adaptive,
    T = time_period_generating_samples,
    n_samples = n_generated_samples,
    diag_s_elements_initial = adaptive_step_run$s_adaptive,
    m_initial = adaptive_step_run$m_adaptive,
    qbar_initial = adaptive_step_run$qbar_end,
    pbar_initial = adaptive_step_run$pbar_end,
    Lambda_initial = 0,
    u_initial = rexp(1),
    rtol = rtol,
    atol = atol,
    h_max = h_max,
    sampling_compute_temporal_averages_of_moments = sampling_compute_temporal_averages_of_moments,
    reflection_type = reflection_type,
    verbose_at_refresh = verbose_at_refresh,    
    last.root.offset.lin.root.finder = last.root.offset.lin.root.finder,
    last.root.offset.non.lin.root.finder = last.root.offset.non.lin.root.finder,
    precision_real_root_lin_root_finder = precision_real_root_lin_root_finder,
    num_subdiv_non_lin_root_finder = num_subdiv_non_lin_root_finder
  )
  
  list(
    adaptive_step_run = adaptive_step_run,
    sample_step_run = sample_step_run
  )
  
}
