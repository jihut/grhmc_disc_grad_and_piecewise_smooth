grhmc_discontinuous_gradient_function_cpp <- function(
    model_list, # this should be a function with argument "settings" defined in the cpp model file
    lambda,
    T = 5000,
    n_samples = 1000,
    q_initial,
    p_initial,
    u_initial, 
    random_state = NULL,
    rtol = NULL,
    atol = NULL,
    h_max = 1.0,
    last_root_offset_lin_root_finder = 1e-10,
    precision_real_root_lin_root_finder = 1.0e-13,
    coeff_tol_real_root_lin_root_finder = 1.0e-13,
    push_eps_tol = 1e-4,
    relevant_indices = NULL,
    return_output_from_ode = FALSE,
    trajectory_info_string = "",
    progress_ratio = 0.1
) {
  
  d <- length(q_initial)
  times <- seq(from = 0, to = T, length.out = n_samples + 1)
  
  if(is.null(relevant_indices)) {
    final_relevant_indices <- 0:(2 * d + 3 - 1) # three additional "variables" in the ODE state of number of momentum events, Lambda and u. 
  } else {
    final_relevant_indices <- c(
      0:2,
      relevant_indices + 3 - 1, # to get out the q's
      relevant_indices + 3 - 1 + d # to get out the p's
    )
  }
  
  settings <- list(
    z0 = c(0, 0, u_initial, q_initial, p_initial),
    times = times,
    h_max = h_max, 
    atol = atol,
    rtol = rtol,
    last_root_offset_lin_root_finder = last_root_offset_lin_root_finder,
    precision_real_root_lin_root_finder = precision_real_root_lin_root_finder,
    coeff_tol = coeff_tol_real_root_lin_root_finder,
    relevant_indices = final_relevant_indices,
    n_parameters = d,
    lambda = lambda,
    push_eps_tol = push_eps_tol,
    trajectory_info = trajectory_info_string,
    progress_ratio = progress_ratio
  )
  
  sim_out <- model_list(
    settings
  )
  
  df_sim_out <- sim_out$samples
  
  df_sim_out  
  if (is.null(relevant_indices)) {
    colnames(df_sim_out) <- c(
      "time",
      "number_of_events",
      "Lambda",
      "u",
      paste0("q", 1:d),
      paste0("p", 1:d)
    )
    
    q_original_samples <- df_sim_out[, 4:(4 + d - 1) + 1][-1, ]
  } else {
    colnames(df_sim_out) <- c(
      "time",
      "number_of_events",
      "Lambda",
      "u",
      paste0("q", relevant_indices),
      paste0("p", relevant_indices)
    )
    
    q_original_samples <- df_sim_out[, 4:(4 + length(relevant_indices) - 1) + 1][-1, ]
  }
  
  if (return_output_from_ode) {
    return(
      list(
        q_original_samples = q_original_samples,
        output_from_ode_solver = df_sim_out,
        n_evals_ode = sim_out$n_evals_ode,
        n_boundary_crossings = sim_out$n_boundary_crossings,
        n_evals = sim_out$n_evals,
        z_final = sim_out$z_final,
        q_final = sim_out$z_final[4:(4 + d - 1)]
      )
    )  
  } else {
    return(
      list(
        q_original_samples = q_original_samples,
        n_evals_ode = sim_out$n_evals_ode,
        n_evals = sim_out$n_evals,
        n_boundary_crossings = sim_out$n_boundary_crossings,
        z_final = sim_out$z_final,
        q_final = sim_out$z_final[4:(4 + d - 1)]
      )
    )
  }
  
  
  
}
