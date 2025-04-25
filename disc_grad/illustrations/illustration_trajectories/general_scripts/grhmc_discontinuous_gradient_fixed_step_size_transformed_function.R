source("disc_grad/illustrations/illustration_trajectories/general_scripts/deSolverRoot_fixed_step_size.R")

grhmc_discontinuous_gradient_fixed_step_size_transformed_function <- function(
    model_list,
    T = 5000,
    n_samples = 1000,
    diag_s_elements_initial,
    m_initial,
    qbar_initial,
    pbar_initial,
    h,
    last.root.offset.lin.root.finder = 1.0e-8,
    last.root.offset.non.lin.root.finder = 1.0e-8,
    precision_real_root_lin_root_finder = 1.0e-13,
    num_subdiv_non_lin_root_finder = 8L
) {
  
  n_evals_ode <- 0L # store how many times the ode function below is called in the process
  
  count_n_evals_ode <- function(){
    
    n_evals_ode <<- n_evals_ode + 1 # increment the counter every time the ode function is called
    
  }
  
  d <- length(qbar_initial)
  model_list$grad_jump_fun(m_initial + diag_s_elements_initial * qbar_initial)
  
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
  # state[1:2]: qbar
  # state[3:4]: pbar
  
  ode <- function(t, state, parms){
    
    count_n_evals_ode()
    
    ret <- c(
      
      
      state[(d + 1):(2 * d)], # \dot qbar = pbar
      diag_s_elements_initial *
        model_list$log_target_grad(m_initial + diag_s_elements_initial * state[1:d]) #\dot pbar = gradient of log transformed density wrt. qbar
      
    )
    
    ret
    
  }
  
  additional_non_lin_root_list <- model_list$additional_non_lin_root_list
  
  if (!is.null(additional_non_lin_root_list)) {
    additional_non_lin_root_fun <- additional_non_lin_root_list$root_fun
    additional_non_lin_event_fun <- additional_non_lin_root_list$event_fun
  }
  
  if (!is.null(additional_non_lin_root_list)) {
    
    final_non_lin_root_fun <- function(t, y, parms) {
      
      additional_non_lin_root_fun(t, y, parms, m_vector = m_initial, s_vector = diag_s_elements_initial, adaptive_yes_or_no = F)
      
    }
    
    final_non_lin_event_fun <- function(t, y, parms) {
      # print("non lin event fun")
      additional_non_lin_event_fun(t, y, parms, m_vector = m_initial, s_vector = diag_s_elements_initial, adaptive_yes_or_no = F)
      
    }
    
  } else {
    
    final_non_lin_root_fun <- NULL
    final_non_lin_event_fun <- NULL
    
  }
  
  additional_lin_root_list <- model_list$additional_lin_root_list
  if (!is.null(additional_lin_root_list)) {
    additional_lin_root_A <- Matrix::Matrix(additional_lin_root_list$A, sparse = TRUE)
    transposed_additional_lin_root_A <- Matrix::t(additional_lin_root_list$A)
    num_of_additional_lin_eqs <- nrow(additional_lin_root_A)    
  }
  
  if (num_of_lin_constraints > 0) {
    region_lin_root_fun <- function(t, y, parms) {
      
      if (is.null(additional_lin_root_list)) {
        
        if (num_of_lin_constraints == 0) {
          
          lin_constraints_root_fun <- NULL
          
        } else if (num_of_lin_constraints == 1) {
          
          lin_constraints_root_fun <- c(diag_s_elements_initial * model_list$region_lin_root_list$A) %*% y[1:d] + model_list$region_lin_root_list$A %*% m_initial + model_list$region_lin_root_list$B
          
        } else {
          
          # lin_constraints_root_fun <- t(diag(diag_s_elements_initial) %*% transposed_model_list_linear_constrained_A) %*% y[1:d] + model_list$region_lin_root_list$A %*% m_initial + model_list$region_lin_root_list$B
          lin_constraints_root_fun <- as.numeric(Matrix::t(sparse_diag_s %*% transposed_model_list_linear_constrained_A) %*% y[1:d]) + as.numeric(model_list_linear_constrained_A %*% m_initial) + model_list$region_lin_root_list$B
          
        }
        
        c(
          lin_constraints_root_fun
        )
        
      } else {
        
        if (num_of_lin_constraints == 0) {
          
          lin_constraints_root_fun <- NULL
          
        } else if (num_of_lin_constraints == 1) {
          
          lin_constraints_root_fun <- c(diag_s_elements_initial * transposed_model_list_linear_constrained_A) %*% y[1:d] + model_list$region_lin_root_list$A %*% m_initial + model_list$region_lin_root_list$B
          
        } else {
          
          lin_constraints_root_fun <- as.numeric(Matrix::t(sparse_diag_s %*% transposed_model_list_linear_constrained_A) %*% y[1:d]) + as.numeric(model_list_linear_constrained_A %*% m_initial) + model_list$region_lin_root_list$B
          
        }
        
        if (num_of_additional_lin_eqs == 1) {
          
          additional_lin_root_fun <- c(diag_s_elements_initial * transposed_additional_lin_root_A) %*% y[1:d] + additional_lin_root_list$A %*% m_initial + additional_lin_root_list$B
          
        } else {
          
          additional_lin_root_fun <- as.numeric(Matrix::t(sparse_diag_s %*% transposed_additional_lin_root_A) %*% y[1:d]) + as.numeric(additional_lin_root_A %*% m_initial) + additional_lin_root_list$B
          
        }
        
        c(
          lin_constraints_root_fun,
          additional_lin_root_fun
        )
        
      }
      
    }
    
    region_lin_event_fun <- function(t, y, parms, lin.root.func) {
      # print("lin event fun")
      yroot <- lin.root.func(t, y, parms)
      
      whichroot <- which(min(abs(yroot)) == abs(yroot))
      
      if (whichroot <= (num_of_lin_constraints)) {
        
        qbar <- y[1:d]
        pbar <- y[(d + 1):(2 * d)] 
        # position slightly after boundary is hit (to be improved in code...): 
        qp <- qbar + 0.0001 * pbar 
        model_list$grad_jump_fun(m_initial + diag_s_elements_initial * qp)
        # q and p unchanged, only id is updated
        # new.state <- y
        # new.state[4 + 6 * d] <- new_region_id
        new.state <- c(qbar, pbar)
        
      } else {
        
        new.state <- additional_lin_root_list$event_fun(t, y, parms, m_vector = m_initial, s_vector = diag_s_elements_initial, adaptive_yes_or_no = F)
        
      }
      
      return(new.state)
      
    }    
  } else {
    region_lin_root_fun <- NULL
    region_lin_event_fun <- NULL
  }
  

  
  
  y0 <- c(
    qbar_initial,
    pbar_initial
  )
  
  sim_out <- deSolverRoot(
    y = y0, 
    times = seq(from = 0, to = T, length.out = n_samples + 1),
    func = ode,
    root.func = final_non_lin_root_fun,
    event.func = final_non_lin_event_fun,
    lin.root.func = region_lin_root_fun,
    lin.root.event.func = region_lin_event_fun,
    last.root.offset.lin.root.finder = last.root.offset.lin.root.finder,
    last.root.offset.non.lin.root.finder = last.root.offset.non.lin.root.finder,
    precision_real_root_lin_root_finder = precision_real_root_lin_root_finder,
    num_subdiv_non_lin_root_finder = num_subdiv_non_lin_root_finder,
    h = h
  )
  
  df_sim_out <- sim_out$samples
  colnames(df_sim_out) <- c(
    "time",
    paste0("qbar", 1:d),
    paste0("pbar", 1:d)
  )
  
  q_original_samples <- t(m_initial + diag(diag_s_elements_initial, nrow = d) %*% t(df_sim_out[, 1:d + 1]))[-1, ]
  
  return(
    list(
      q_original_samples = q_original_samples,
      output_from_ode_solver = df_sim_out,
      n_evals_ode = n_evals_ode,
      s_elements = diag_s_elements_initial,
      m_elements = m_initial,
      q_final = df_sim_out[nrow(df_sim_out), 1:d + 1],
      p_final = df_sim_out[nrow(df_sim_out), (d + 2):(2 * d + 1)],
      sim_output = sim_out
    )
  )
  
}
