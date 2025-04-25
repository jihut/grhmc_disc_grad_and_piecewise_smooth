source("piecewise_smooth/implementation_scripts/general_scripts/deSolverRoot.R")

grhmc_piecewise_smooth_density_transformed_function <- function(
    model_list,
    lambda,
    T = 5000,
    n_samples = 1000,
    diag_s_elements_initial,
    m_initial,
    qbar_initial,
    pbar_initial,
    Lambda_initial,
    u_initial,
    random_state = NULL,
    rtol = NULL,
    atol = NULL,
    h_max = 1.0,
    sampling_compute_temporal_averages_of_moments = FALSE,
    reflection_type = NULL, 
    verbose_at_refresh = FALSE,
    last.root.offset.lin.root.finder = 1e-10,
    last.root.offset.non.lin.root.finder = 1e-10,
    precision_real_root_lin_root_finder = 1.0e-13,
    num_subdiv_non_lin_root_finder = 8L
) {

  stopifnot(reflection_type == "deterministic" | reflection_type == "randomized_dense" | reflection_type == "randomized_sparse")
  
  if(!is.null(random_state)){
    set.seed(random_state)
  }
  
  n_evals_ode <- 0L # store how many times the ode function below is called in the process
  
  count_n_evals_ode <- function(){
    
    n_evals_ode <<- n_evals_ode + 1 # increment the counter every time the ode function is called
    
  }
  
  d <- length(qbar_initial)
  model_list$target_jump_fun(m_initial + diag_s_elements_initial * qbar_initial)
  sampling_compute_temporal_averages_of_moments <<- sampling_compute_temporal_averages_of_moments
  
  # transposed_model_list_linear_constrained_A <- t(model_list$region_lin_root_list$A) # New addition: Use sparse matrix as a linear equation defining the boundary of change in gradient tends to concern only a certain amount of coordinates
  if (!is.null(model_list$region_lin_root_list)) {
    num_of_lin_constraints <- nrow(model_list$region_lin_root_list$A)
    model_list_linear_constrained_A <- Matrix::Matrix(model_list$region_lin_root_list$A, sparse = TRUE)
    transposed_model_list_linear_constrained_A <- Matrix::t(model_list_linear_constrained_A)
  } else {
    num_of_lin_constraints <- 0
  }
  sparse_diag_s <- Matrix::sparseMatrix(i = 1:d, j = 1:d, x = diag_s_elements_initial)
  
  # define the (first order) ode for state
  # state[1]: number of events up to time point
  # state[2]: Lambda, which resets after each event
  # state[3]: u, simulate from exp(1) after each event
  # state[4:(4 + n_dim - 1)]: qbar
  # state[(4 + n_dim):(4 + 2*n_dim - 1)]: pbar
  # state[(4 + 2*n_dim):(4 + 3*n_dim - 1)]: \int qbar dt
  # state[(4 + 3*n_dim):(4 + 4*n_dim - 1)]: \int qbar^2 dt
  # state[(4 + 4*n_dim):(4 + 5*n_dim - 1)]: \int q dt = \int (m + Sqbar) dt
  # state[(4 + 5*n_dim):(4 + 6*n_dim - 1)]: \int q^2 dt = \int (m + Sqbar)^2 dt
  
  if (sampling_compute_temporal_averages_of_moments) {
    
    ode <- function(t, state, parms){
      
      count_n_evals_ode()
      
      ret <- c(
        
        0, # number of momentum refresh events
        lambda, # integrate lambda to get Lambda
        0, # u
        state[(4 + d):(4 + 2 * d - 1)], # \dot qbar = pbar
        diag_s_elements_initial *
          model_list$log_target_grad(m_initial + diag_s_elements_initial * state[4:(4 + d - 1)]), #\dot pbar = gradient of log transformed density wrt. qbar
        state[4:(4 + d - 1)], # \int qbar dt
        state[4:(4 + d - 1)] ^ 2, # \int qbar^2 dt
        m_initial + diag_s_elements_initial * state[4:(4 + d - 1)], # \int q dt = \int m + Sqbar dt
        (m_initial + diag_s_elements_initial * state[4:(4 + d - 1)]) ^ 2 # \int q^2 dt = \int (m + Sqbar) dt
      )
      
      ret
      
    }  
    
  } else {
    
    ode <- function(t, state, parms){
      
      count_n_evals_ode()
      
      ret <- c(
        
        0, # number of momentum refresh events
        lambda, # integrate lambda to get Lambda
        0, # u
        state[(4 + d):(4 + 2 * d - 1)], # \dot qbar = pbar
        diag_s_elements_initial *
          model_list$log_target_grad(m_initial + diag_s_elements_initial * state[4:(4 + d - 1)]) #\dot pbar = gradient of log transformed density wrt. qbar
      )
      
      ret
      
    }
    
  }
  
  additional_non_lin_root_list <- model_list$additional_non_lin_root_list
  
  if (!is.null(additional_non_lin_root_list)) {
    additional_non_lin_root_fun <- additional_non_lin_root_list$root_fun
    additional_non_lin_event_fun <- additional_non_lin_root_list$event_fun
  }
  
  if (!is.null(additional_non_lin_root_list)) { # in case there are any other predefined non linear root and event functions, e.g. related to some constraints etc. 
    
    final_non_lin_root_fun <- function(t, y, parms) {
      
      additional_non_lin_root_fun(t, y, parms, m_vector = m_initial, s_vector = diag_s_elements_initial, adaptive_yes_or_no = F)
      
    }
    
    final_non_lin_event_fun <- function(t, y, parms) {
      
      additional_non_lin_event_fun(t, y, parms, m_vector = m_initial, s_vector = diag_s_elements_initial, adaptive_yes_or_no = F)
      
    }
    
  } else {
    
    final_non_lin_root_fun <- NULL
    final_non_lin_event_fun <- NULL
    
  }
  
  additional_lin_root_list <- model_list$additional_lin_root_list # in case there are any other predefined linear root and event functions, e.g. related to some constraints etc. 
  if (!is.null(additional_lin_root_list)) {
    additional_lin_root_A <- Matrix::Matrix(additional_lin_root_list$A, sparse = TRUE)
    transposed_additional_lin_root_A <- Matrix::t(additional_lin_root_A)
    num_of_additional_lin_eqs <- nrow(additional_lin_root_A)    
  }
  
  region_lin_root_fun <- function(t, y, parms) {
    
    if (is.null(additional_lin_root_list)) { # if no additional linear root functions are given
      
      if (num_of_lin_constraints == 0) {
        
        lin_constraints_root_fun <- NULL
        
      } else if (num_of_lin_constraints == 1) { # no need for matrix if only one single linear root function related to jump in density
        
        lin_constraints_root_fun <- (diag_s_elements_initial * model_list$region_lin_root_list$A) %*% y[4:(4 + d - 1)] + model_list$region_lin_root_list$A %*% m_initial + model_list$region_lin_root_list$B
        
      } else {
        
        # lin_constraints_root_fun <- t(diag(diag_s_elements_initial) %*% transposed_model_list_linear_constrained_A) %*% y[4:(4 + d - 1)] + model_list$region_lin_root_list$A %*% m_initial + model_list$region_lin_root_list$B
        lin_constraints_root_fun <- as.numeric(Matrix::t(sparse_diag_s %*% transposed_model_list_linear_constrained_A) %*% y[4:(4 + d - 1)]) + as.numeric(model_list_linear_constrained_A %*% m_initial) + model_list$region_lin_root_list$B
        
      }
      
      c(
        y[2] - y[3], # momentum refresh 
        # diag(diag_s_elements_initial) %*% model_list$region_lin_root_list$A %*% y[4:(4 + d - 1)] + model_list$region_lin_root_list$A %*% m_initial + model_list$region_lin_root_list$B 
        lin_constraints_root_fun # root function to detect crossing boundary between two densities
      )
      
    } else { # if additional linear root functions are given, similar as above, but also an extra set of additional linear root functions
      
      if (num_of_lin_constraints == 0) {
        
        lin_constraints_root_fun <- NULL
        
      } else if (num_of_lin_constraints == 1) {
        
        lin_constraints_root_fun <- (diag_s_elements_initial * model_list$region_lin_root_list$A) %*% y[4:(4 + d - 1)] + model_list$region_lin_root_list$A %*% m_initial + model_list$region_lin_root_list$B
        
      } else {
        
        # lin_constraints_root_fun <- t(diag(diag_s_elements_initial) %*% transposed_model_list_linear_constrained_A) %*% y[4:(4 + d - 1)] + model_list$region_lin_root_list$A %*% m_initial + model_list$region_lin_root_list$B
        lin_constraints_root_fun <- as.numeric(Matrix::t(sparse_diag_s %*% transposed_model_list_linear_constrained_A) %*% y[4:(4 + d - 1)]) + as.numeric(model_list_linear_constrained_A %*% m_initial) + model_list$region_lin_root_list$B
        
      }
      
      if (num_of_additional_lin_eqs == 1) {
        
        additional_lin_root_fun <- (diag_s_elements_initial * additional_lin_root_list$A) %*% y[4:(4 + d - 1)] + additional_lin_root_list$A %*% m_initial + additional_lin_root_list$B
        
      } else {
        
        # lin_constraints_root_fun <- t(diag(diag_s_elements_initial) %*% transposed_model_list_linear_constrained_A) %*% y[4:(4 + d - 1)] + model_list$region_lin_root_list$A %*% m_initial + model_list$region_lin_root_list$B
        additional_lin_root_fun <- as.numeric(Matrix::t(sparse_diag_s %*% transposed_additional_lin_root_A) %*% y[4:(4 + d - 1)]) + as.numeric(additional_lin_root_A %*% m_initial) + additional_lin_root_list$B
        
      }
      
      c(
        y[2] - y[3],
        # diag(diag_s_elements_initial) %*% model_list$region_lin_root_list$A %*% y[4:(4 + d - 1)] + model_list$region_lin_root_list$A %*% m_initial + model_list$region_lin_root_list$B 
        lin_constraints_root_fun,
        additional_lin_root_fun
      )
      
    }
    
  }
  
  region_lin_event_fun <- function(t, y, parms, lin.root.func) {
    
    yroot <- lin.root.func(t, y, parms)
    
    whichroot <- which(min(abs(yroot)) == abs(yroot))
    
    if (whichroot == 1) { # if root due to first root function --> momentum refresh
      
      if (verbose_at_refresh) {
        # print("#######################")
        # print(paste0("refresh t: ", t))
        print(paste0("t: ", t))
      }
      
      # new.state <- y
      # new.state[1] <- new.state[1] + 1
      # new.state[2] <- 0
      # new.state[3] <- rexp(1)
      # new.state[(4 + d):(4 + 2 * d - 1)] <- rnorm(d)
      
      if (sampling_compute_temporal_averages_of_moments) {
        new.state <- c(
          y[1] + 1,
          0,
          rexp(1),
          y[4:(4 + d - 1)],
          rnorm(d),
          y[(4 + 2 * d):(4 + 6 * d - 1)]
        )  
      } else {
        new.state <- c(
          y[1] + 1,
          0,
          rexp(1),
          y[4:(4 + d - 1)],
          rnorm(d)
        )
      }
      
      
    } else if (whichroot > 1 & whichroot <= (num_of_lin_constraints + 1)) { # due to crossing boundary
      # print("#######################")
      # print("Boundary crossing")
      # print(paste0("Current region: ", region_id))
      # print(paste0("t: ", t))
      
      lin_constraint_index <- whichroot - 1
      
      # qbar <- y[4:(4 + d - 1)]
      # pbar <- y[(4 + d):(4 + 2 * d - 1)] 
      # # position slightly after boundary is hit (to be improved in code...): 
      # qp <- qbar + 0.0001 * pbar 
      # model_list$grad_jump_fun(m_initial + diag_s_elements_initial * qp)
      # # q and p unchanged, only id is updated
      # # new.state <- y
      # # new.state[4 + 7 * d] <- new_region_id
      # new.state <- c(y[1:3], qbar, pbar, y[(4 + 2 * d):(4 + 7 * d - 1)])
      
      old_qbar <- y[4:(4 + d - 1)]
      old_pbar <- y[(4 + d):(4 + 2 * d - 1)]
      
      # print(paste0("Current qbar: ", old_qbar))
      # print(paste0("Current pbar: ", old_pbar))
      
      old_q <- m_initial + diag_s_elements_initial * old_qbar
      old_q[abs(old_q) <= 1e-13] <- 0
      # print(paste0("Current q: ", old_q))
      
      
      new_qbar_eps <- old_qbar + old_pbar * 1e-5
      old_qbar_eps <- old_qbar - old_pbar * 1e-5
      
      old_potential_energy <- -model_list$log_target_fun(old_q) # potential energy just before hitting the boundary
      # print(paste0("old_potential_energy: ", old_potential_energy))
      
      model_list$target_jump_fun(m_initial + diag_s_elements_initial * new_qbar_eps) # make sure that one is in the new region
      
      new_potential_energy <- -model_list$log_target_fun(old_q) # log det(S) cancel out? - potential energy right after switching region
      # print(paste0("new_potential_energy: ", new_potential_energy))
      
      delta_U <- new_potential_energy - old_potential_energy # difference between the potential energy in the two regions
      
      # print(paste0("delta_U: ", delta_U))
      
      if (num_of_lin_constraints != 1) {
        normal_vec <- diag_s_elements_initial * model_list$region_lin_root_list$A[whichroot - 1, ]
      } else {
        normal_vec <- diag_s_elements_initial * model_list$region_lin_root_list$A
      }
      # print(paste0("normal vec: ", normal_vec))
      
      # Next: Find the projection of pbar onto the normal vec
      old_pbar_perpendicular <- sum(old_pbar * normal_vec) / sum(normal_vec ^ 2) * normal_vec
      # print(paste0("pbar_perpendicular: ", old_pbar_perpendicular))
      old_pbar_parallel <- old_pbar - old_pbar_perpendicular
      # print(paste0("pbar_parallel: ", old_pbar_parallel))
      
      norm_squared_old_pbar_perpendicular <- sum(old_pbar_perpendicular ^ 2)
      # print(paste0("norm_squared_old_pbar_perpendicular: ", norm_squared_old_pbar_perpendicular))
      
      if (norm_squared_old_pbar_perpendicular > 2 * delta_U) { # if the momentum squared along the normal direction is larger than two times the potential energy --> transition to new region
        
        # print("Cross boundary")
        
        new_pbar_perpendicular <- sqrt(norm_squared_old_pbar_perpendicular - 2 * delta_U) * old_pbar_perpendicular / sqrt(norm_squared_old_pbar_perpendicular) # refract the momentum along the normal direction after moving to a different region
        
        new_pbar <- old_pbar_parallel + new_pbar_perpendicular
        # print(paste0("new pbar:", new_pbar))
        
        if (sampling_compute_temporal_averages_of_moments) {
          new.state <- c(
            y[1:3],
            old_qbar,
            # new_qbar_eps,
            new_pbar,
            y[(4 + 2 * d):(4 + 6 * d - 1)]
          )          
        } else {
          new.state <- c(
            y[1:3],
            old_qbar,
            # new_qbar_eps,
            new_pbar
          )
        }

        
        # print(paste0("New region: ", region_id))
        
      } else { # if not, then reflect the component of the momentum along the normal direction, either by deterministic or randomized reflection
        # print("No cross boundary")
        
        if (reflection_type == "deterministic") {
          
          new_pbar_perpendicular <- -old_pbar_perpendicular
          
          new_pbar <- old_pbar_parallel + new_pbar_perpendicular 
          
        } else if (reflection_type == "randomized_dense") {
          
          z <- rnorm(d)
          
          new_pbar <- z - sum((old_pbar + z) * normal_vec) / sum(normal_vec ^ 2) * normal_vec
          
        } else if (reflection_type == "randomized_sparse") {
          
          non_zero_normal_vec_comps <- which(normal_vec != 0)
          
          new_pbar <- old_pbar
          
          z <- rnorm(length(non_zero_normal_vec_comps))
          
          new_pbar[non_zero_normal_vec_comps] <- z - sum((old_pbar[non_zero_normal_vec_comps] + z) * normal_vec[non_zero_normal_vec_comps]) / sum(normal_vec[non_zero_normal_vec_comps] ^ 2) * normal_vec[non_zero_normal_vec_comps]
          
        }
        
        new_qbar <- old_qbar + new_pbar * 1e-5
        # new_qbar <- old_qbar_eps + new_pbar * 1e-10
        # print(paste0("new pbar:", new_pbar))
        
        model_list$target_jump_fun(old_qbar_eps * diag_s_elements_initial + m_initial) # ensure that the region just before hit is relevant again
        # print(paste0("New region: ", region_id))
        
        if (sampling_compute_temporal_averages_of_moments) {
          
          new.state <- c(
            y[1:3],
            old_qbar,
            # new_qbar,
            new_pbar,
            y[(4 + 2 * d):(4 + 6 * d - 1)]
          )
          
        } else {
          
          new.state <- c(
            y[1:3],
            old_qbar,
            # new_qbar,
            new_pbar
          )
          
        }
        
      }
      
    } else {
      
      new.state <- additional_lin_root_list$event_fun(t, y, parms, m_vector = m_initial, s_vector = diag_s_elements_initial, adaptive_yes_or_no = F)
      
    }
    
    return(new.state)
    
  }
  
  if (sampling_compute_temporal_averages_of_moments) {
    
    y0 <- c(
      0,
      Lambda_initial,
      u_initial,
      qbar_initial,
      pbar_initial,
      rep(0, d),
      rep(0, d),
      rep(0, d),
      rep(0, d)
    )  
    
  } else {
    
    y0 <- c(
      0,
      Lambda_initial,
      u_initial,
      qbar_initial,
      pbar_initial
    )
    
  }
  
  if(is.null(rtol)){
    rtol <- 1e-4 # default in deSolve::lsodar
  } else {
    rtol <- rtol
  }
  
  if(is.null(atol)){
    atol <- 1e-4 # default in deSolve::lsodar
  } else {
    atol <- atol
  }
  
  if(is.null(h_max)) {
    h_max <- 0.5
  } else {
    h_max <- h_max
  }
  
  sim_out <- deSolverRoot(
    y = y0, 
    times = seq(from = 0, to = T, length.out = n_samples + 1),
    func = ode,
    root.func = final_non_lin_root_fun,
    event.func = final_non_lin_event_fun,
    lin.root.func = region_lin_root_fun,
    lin.root.event.func = region_lin_event_fun,
    h.max = h_max,
    last.root.offset.lin.root.finder = last.root.offset.lin.root.finder,
    last.root.offset.non.lin.root.finder = last.root.offset.non.lin.root.finder,
    precision_real_root_lin_root_finder = precision_real_root_lin_root_finder,
    num_subdiv_non_lin_root_finder = num_subdiv_non_lin_root_finder
  )
  
  df_sim_out <- sim_out$samples
  
  if (sampling_compute_temporal_averages_of_moments) {
    
    colnames(df_sim_out) <- c(
      "time",
      "number_of_events",
      "Lambda",
      "u",
      paste0("qbar", 1:d),
      paste0("pbar", 1:d),
      paste0("int_qbar", 1:d),
      paste0("int_qbar", 1:d, "_squared"),
      paste0("int_q", 1:d),
      paste0("int_q", 1:d, "_squared")
    )
    
  } else {
    
    colnames(df_sim_out) <- c(
      "time",
      "number_of_events",
      "Lambda",
      "u",
      paste0("qbar", 1:d),
      paste0("pbar", 1:d)
    )
    
  }
  
  q_original_samples <- t(m_initial + diag(diag_s_elements_initial, nrow = d) %*% t(df_sim_out[, 4:(4 + d - 1) + 1]))[-1, ]
  
  return(
    list(
      q_original_samples = q_original_samples,
      output_from_ode_solver = df_sim_out,
      n_evals_ode = n_evals_ode,
      s_elements = diag_s_elements_initial,
      m_elements = m_initial,
      lambda = lambda,
      random_state = random_state,
      rtol = rtol, 
      atol = atol
    )
  )
  
}
