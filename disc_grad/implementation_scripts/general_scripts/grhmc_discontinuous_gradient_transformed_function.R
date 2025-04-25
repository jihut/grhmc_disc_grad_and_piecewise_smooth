source("disc_grad/implementation_scripts/general_scripts/deSolverRoot.R")

grhmc_discontinuous_gradient_transformed_function <- function(
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
    last.root.offset.lin.root.finder = 1.0e-8,
    last.root.offset.non.lin.root.finder = 1.0e-8,
    precision_real_root_lin_root_finder = 1.0e-13,
    num_subdiv_non_lin_root_finder = 8L,
    verbose_at_refresh = FALSE
) {
  
  if(!is.null(random_state)){
    set.seed(random_state)
  }
  
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
  # state[1]: number of events up to time point
  # state[2]: Lambda, which resets after each event
  # state[3]: u, simulate from exp(1) after each event
  # state[4:(4 + n_dim - 1)]: qbar
  # state[(4 + n_dim):(4 + 2*n_dim - 1)]: pbar
  # state[(4 + 2*n_dim):(4 + 3*n_dim - 1)]: \int qbar dt
  # state[(4 + 3*n_dim):(4 + 4*n_dim - 1)]: \int qbar^2 dt
  # state[(4 + 4*n_dim):(4 + 5*n_dim - 1)]: \int q dt = \int (m + Sqbar) dt
  # state[(4 + 5*n_dim):(4 + 6*n_dim - 1)]: \int q^2 dt = \int (m + Sqbar)^2 dt
  # state[4 + 6*n_dim]: region id
  
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
    transposed_additional_lin_root_A <- Matrix::t(additional_lin_root_list$A)
    num_of_additional_lin_eqs <- nrow(additional_lin_root_A)    
  }
  
  region_lin_root_fun <- function(t, y, parms) {
    
    if (is.null(additional_lin_root_list)) { # if no additional linear root functions are given
      
      if (num_of_lin_constraints == 0) {
        
        lin_constraints_root_fun <- NULL
        
      } else if (num_of_lin_constraints == 1) { # no need for matrix if only one single linear root function related to discontinuous gradient
        
        lin_constraints_root_fun <- c(diag_s_elements_initial * model_list$region_lin_root_list$A) %*% y[4:(4 + d - 1)] + model_list$region_lin_root_list$A %*% m_initial + model_list$region_lin_root_list$B
        
      } else {
        
        # lin_constraints_root_fun <- t(diag(diag_s_elements_initial) %*% transposed_model_list_linear_constrained_A) %*% y[4:(4 + d - 1)] + model_list$region_lin_root_list$A %*% m_initial + model_list$region_lin_root_list$B
        lin_constraints_root_fun <- as.numeric(Matrix::t(sparse_diag_s %*% transposed_model_list_linear_constrained_A) %*% y[4:(4 + d - 1)]) + as.numeric(model_list_linear_constrained_A %*% m_initial) + model_list$region_lin_root_list$B
        
      }
      
      c(
        y[2] - y[3], # momentum refresh 
        # diag(diag_s_elements_initial) %*% model_list$region_lin_root_list$A %*% y[4:(4 + d - 1)] + model_list$region_lin_root_list$A %*% m_initial + model_list$region_lin_root_list$B 
        lin_constraints_root_fun # root function to detect crossing boundary of discontinuous gradient
      )
      
    } else { # if additional linear root functions are given, similar as above, but also an extra set of additional linear root functions
      
      if (num_of_lin_constraints == 0) {
        
        lin_constraints_root_fun <- NULL
        
      } else if (num_of_lin_constraints == 1) {
        
        lin_constraints_root_fun <- c(diag_s_elements_initial * transposed_model_list_linear_constrained_A) %*% y[4:(4 + d - 1)] + model_list$region_lin_root_list$A %*% m_initial + model_list$region_lin_root_list$B
        
      } else {
        
        # lin_constraints_root_fun <- t(diag(diag_s_elements_initial) %*% transposed_model_list_linear_constrained_A) %*% y[4:(4 + d - 1)] + model_list$region_lin_root_list$A %*% m_initial + model_list$region_lin_root_list$B
        lin_constraints_root_fun <- as.numeric(Matrix::t(sparse_diag_s %*% transposed_model_list_linear_constrained_A) %*% y[4:(4 + d - 1)]) + as.numeric(model_list_linear_constrained_A %*% m_initial) + model_list$region_lin_root_list$B
        
      }
      
      if (num_of_additional_lin_eqs == 1) {
        
        additional_lin_root_fun <- c(diag_s_elements_initial * transposed_additional_lin_root_A) %*% y[4:(4 + d - 1)] + additional_lin_root_list$A %*% m_initial + additional_lin_root_list$B
        
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
        
        print(paste0("t: ", t))
        
      }
      
      # new.state <- y
      # new.state[1] <- new.state[1] + 1
      # new.state[2] <- 0
      # new.state[3] <- rexp(1)
      # new.state[(4 + d):(4 + 2 * d - 1)] <- rnorm(d)
      
      new.state <- c(
        y[1] + 1,
        0,
        rexp(1),
        y[4:(4 + d - 1)],
        rnorm(d),
        y[(4 + 2 * d):(4 + 6 * d - 1)]
      )
      
    } else if (whichroot > 1 & whichroot <= (num_of_lin_constraints + 1)) { # due to crossing boundary of discontinuous gradient
      
      qbar <- y[4:(4 + d - 1)]
      pbar <- y[(4 + d):(4 + 2 * d - 1)] 
      # position slightly after boundary is hit (to be improved in code...): 
      qp <- qbar + 0.0001 * pbar 
      model_list$grad_jump_fun(m_initial + diag_s_elements_initial * qp)
      # q and p unchanged, only id is updated
      # new.state <- y
      # new.state[4 + 6 * d] <- new_region_id
      new.state <- c(y[1:3], qbar, pbar, y[(4 + 2 * d):(4 + 6 * d - 1)])
      
    } else { # due to the additional linear root functions that are given
      
      new.state <- additional_lin_root_list$event_fun(t, y, parms, m_vector = m_initial, s_vector = diag_s_elements_initial, adaptive_yes_or_no = F)
      
    }
    
    return(new.state)
    
  }
  
  
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
    last.root.offset.lin.root.finder = last.root.offset.lin.root.finder,
    last.root.offset.non.lin.root.finder = last.root.offset.non.lin.root.finder,
    precision_real_root_lin_root_finder = precision_real_root_lin_root_finder,
    num_subdiv_non_lin_root_finder = num_subdiv_non_lin_root_finder,
    h.max = h_max
  )
  
  df_sim_out <- sim_out$samples
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
