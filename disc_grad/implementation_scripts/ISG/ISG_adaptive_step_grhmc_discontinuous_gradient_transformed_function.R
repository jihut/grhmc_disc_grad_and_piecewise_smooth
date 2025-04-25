source("disc_grad/implementation_scripts/general_scripts/deSolverRoot.R")

ISG_adaptive_step_grhmc_discontinuous_gradient_transformed_function <- function( # tune both lambda (event rate), m and diagonal S at the same time
  model_list,
  lambda,
  lambda_lower_limit = 0.001,
  T = 5000,
  n_samples = 1000,
  m_initial = NULL,
  diag_s_elements_initial = NULL,
  qbar_initial = NULL,
  pbar_initial = NULL,
  random_state = NULL,
  rtol = NULL,
  atol = NULL,
  h_max = 1.0,
  proportion_time_until_adaptive_start = 0.05,
  min_proportion_of_previous_state = 0.5,
  max_proportion_of_previous_state = 2,
  last.root.offset.lin.root.finder = 1.0e-8,
  last.root.offset.non.lin.root.finder = 1.0e-8,
  precision_real_root_lin_root_finder = 1.0e-13,
  num_subdiv_non_lin_root_finder = 8L,
  verbose_at_refresh = FALSE
) {
  
  if(!is.null(random_state)){
    set.seed(random_state)
  }
  
  # n_dim_error <- FALSE
  
  # tryCatch(
  #   {
  #     d <- length(model_list$log_target_grad(0, 1))
  #   },
  #   error = function(e){
  #     print(e)
  #     n_dim_error <<- TRUE
  #   }
  # )
  # 
  # if(n_dim_error == TRUE){
  #   if(is.null(n_parameters)){
  #     stop("Please specify n_parameters")
  #   } else {
  #     d <- n_parameters
  #   }
  # }
  
  d <- model_list$n_parameters
  
  if(is.null(m_initial)){
    m_adaptive <- rep(0, d) # initialize m 
  } else {
    m_adaptive <- m_initial
  }
  old_m_adaptive <- m_adaptive
  
  if(is.null(diag_s_elements_initial)){
    s_adaptive <- rep(1, d) # initialize m 
  } else {
    s_adaptive <- diag_s_elements_initial
  }
  old_s_adaptive <- s_adaptive
  
  if(is.null(qbar_initial)){
    # qbar_initial = rep(0, d)
    qbar_initial <- rnorm(d) # If initial values of qbar are not specified --> simulate randomly from a multivariate normal distribution with zero mean vector and diagonal covariance matrix
  } else {
    qbar_initial <- qbar_initial
  }
  
  if(is.null(pbar_initial)){
    pbar_initial <- rnorm(d)
  } else {
    pbar_initial <- pbar_initial
  }
  u_initial <- rexp(1)
  
  q_initial <- m_adaptive + s_adaptive * qbar_initial
  model_list$grad_jump_fun(q_initial)
  
  # transposed_model_list_linear_constrained_A <- t(model_list$region_lin_root_list$A) # New addition: Use sparse matrix as a linear equation defining the boundary of change in gradient tends to concern only a certain amount of coordinates
  if (!is.null(model_list$region_lin_root_list)) {
    num_of_lin_constraints <- nrow(model_list$region_lin_root_list$A)
    model_list_linear_constrained_A <- Matrix::Matrix(model_list$region_lin_root_list$A, sparse = TRUE)
    transposed_model_list_linear_constrained_A <- Matrix::t(model_list_linear_constrained_A)
  } else {
    num_of_lin_constraints <- 0
  }
  
  sparse_diag_s <- Matrix::sparseMatrix(i = 1:d, j = 1:d, x = s_adaptive)
  
  # Quantities for adaptive tuning of lambda 
  
  lambda_adaptive <<- lambda
  lambda_lower_limit <<- lambda_lower_limit
  nut_time_found_indicator <<- FALSE
  current_qbar <<- qbar_initial
  t_at_current_qbar <<- 0
  sum_all_nut_times <<- 0
  n_uncensored_nut_times <<- 0
  proportion_time_until_adaptive_start <<- proportion_time_until_adaptive_start
  # reflection_indicator <<- FALSE
  
  # define the (first order) ode for state
  # state[1]: number of events up to time point
  # state[2]: Lambda, which resets after each event
  # state[3]: u, simulate from exp(1) after each event
  # state[4:(4 + d - 1)]: qbar
  # state[(4 + d):(4 + 2*d - 1)]: pbar
  # state[(4 + 2*d):(4 + 3*d - 1)]: \int qbar dt
  # state[(4 + 3*d):(4 + 4*d - 1)]: \int qbar^2 dt
  # state[(4 + 4*d):(4 + 5*d - 1)]: \int q dt = \int (m + Sqbar) dt
  # state[(4 + 5*d):(4 + 6*d - 1)]: \int q^2 dt = \int (m + Sqbar)^2 dt
  # state[(4 + 6*d):(4 + 7*d - 1)]: \int gradient of original log target density ^2 dt
  # state[4 + 7*d]: region id
  
  ode <- function(t, y, parms){
    
    # In ISG: we need to use the log target gradient in the tuning of mass matrix as well (integrated squared gradients).
    # More efficient to do this than defining evaluation of log_target_grad twice directly in the set of differential equations.
    
    eval_log_target_grad <- model_list$log_target_grad(m_adaptive + s_adaptive * y[4:(4 + d - 1)])
    
    ret <- c(
      
      0, # number of momentum refresh events
      lambda_adaptive, # integrate lambda to get Lambda
      0, # u
      y[(4 + d):(4 + 2 * d - 1)], # \dot qbar = pbar
      s_adaptive *
        eval_log_target_grad, #\dot pbar = gradient of log transformed density wrt. qbar
      y[4:(4 + d - 1)], # \int qbar dt
      y[4:(4 + d - 1)] ^ 2, # \int qbar^2 dt
      m_adaptive + s_adaptive * y[4:(4 + d - 1)], # \int q dt = \int m + Sqbar dt
      (m_adaptive + s_adaptive * y[4:(4 + d - 1)]) ^ 2, # \int q^2 dt = \int (m + Sqbar) dt
      eval_log_target_grad ^ 2
    )
    
    ret
    
  }
  
  additional_lin_root_list <- model_list$additional_lin_root_list
  if (!is.null(additional_lin_root_list)) {
    additional_lin_root_A <- Matrix::Matrix(additional_lin_root_list$A, sparse = TRUE)
    transposed_additional_lin_root_A <- Matrix::t(additional_lin_root_list$A)
    num_of_additional_lin_eqs <- nrow(additional_lin_root_A)    
  }
  
  
  additional_non_lin_root_list <- model_list$additional_non_lin_root_list
  
  if (!is.null(additional_lin_root_list)) {
    additional_non_lin_root_fun <- model_list$additional_non_lin_root_list$root_fun
    additional_non_lin_event_fun <- model_list$additional_non_lin_root_list$event_fun
  }
  
  
  
  
  # Define non-linear root and event functions
  
  final_non_lin_root_fun <- function(t, y, parms) {
    
    if (is.null(additional_non_lin_root_list)) { # if no other non linear root functions are given 
      
      if (nut_time_found_indicator) { # nut time has been found during a period between two momentum updates, stop finding nut time
        
        1
        
      } else {
        
        # sum((y[4:(4 + d - 1)] - current_qbar) * y[(4 + d):(4 + 2 * d - 1)]) / # nut condition
        #   sqrt(sum((y[4:(4 + d - 1)] - current_qbar) ^ 2) * sum(y[(4 + d):(4 + 2 * d - 1)] ^ 2) + 1.0e-8) # scaling factor due to numerical purposes
        sum((y[4:(4 + d - 1)] - current_qbar) * y[(4 + d):(4 + 2 * d - 1)])
      }
      
    } else {
      
      if (nut_time_found_indicator) {
        
        c(
          1,
          additional_non_lin_root_fun(t, y, parms, m_vector = m_adaptive, s_vector = s_adaptive, adaptive_yes_or_no = T)
        )
        
      } else {
        
        c(
          # sum((y[4:(4 + d - 1)] - current_qbar) * y[(4 + d):(4 + 2 * d - 1)]) / # nut condition
          #   sqrt(sum((y[4:(4 + d - 1)] - current_qbar) ^ 2) * sum(y[(4 + d):(4 + 2 * d - 1)] ^ 2) + 1.0e-8), # scaling factor due to numerical purposes
          sum((y[4:(4 + d - 1)] - current_qbar) * y[(4 + d):(4 + 2 * d - 1)]),
          additional_non_lin_root_fun(t, y, parms, m_vector = m_adaptive, s_vector = s_adaptive, adaptive_yes_or_no = T)
        )
        
      }
      
    }
    
  }
  
  final_non_lin_event_fun <- function(t, y, parms) {
    
    if (t > 0) {
      
      if (is.null(additional_non_lin_root_list)) {
        
        if (nut_time_found_indicator) {
          
          new.state <- y
          
        } else { # if nut time has not been found yet
          
          sum_all_nut_times <<- sum_all_nut_times + (t - t_at_current_qbar)
          n_uncensored_nut_times <<- n_uncensored_nut_times + 1
          nut_time_found_indicator <<- TRUE
          
          new.state <- y
          
        }
        
      } else {
        
        yroot <- final_non_lin_root_fun(t, y, parms)
        whichroot <- which(min(abs(yroot)) == abs(yroot))
        
        if (whichroot == 1) {
          
          if (nut_time_found_indicator) {
            
            new.state <- y
            
          } else { # if nut time has not been found yet
            
            sum_all_nut_times <<- sum_all_nut_times + (t - t_at_current_qbar)
            n_uncensored_nut_times <<- n_uncensored_nut_times + 1 # nut event has happened
            nut_time_found_indicator <<- TRUE
            
            new.state <- y
            
          }
          
        } else { # root due to the additional non linear root functions
          
          new.state <- additional_non_lin_event_fun(t, y, parms, m_vector = m_adaptive, s_vector = s_adaptive, adaptive_yes_or_no = T)
          
        }
        
      }
      
    } else {
      
      new.state <- y
      
    }
    
    new.state
    
  }
  
  region_lin_root_fun <- function(t, y, parms) {
    
    if (is.null(additional_lin_root_list)) { # if no additional linear root functions are given
      
      if (num_of_lin_constraints == 0) {
        
        lin_constraints_root_fun <- NULL
        
      } else if (num_of_lin_constraints == 1) { # no need for matrix if only one single linear root function related to discontinuous gradient
        
        lin_constraints_root_fun <- c(s_adaptive * model_list$region_lin_root_list$A) %*% y[4:(4 + d - 1)] + model_list$region_lin_root_list$A %*% m_adaptive + model_list$region_lin_root_list$B
        
      } else {
        
        # lin_constraints_root_fun <- t(diag(s_adaptive) %*% transposed_model_list_linear_constrained_A) %*% y[4:(4 + d - 1)] + model_list$region_lin_root_list$A %*% m_adaptive + model_list$region_lin_root_list$B
        lin_constraints_root_fun <- as.numeric(Matrix::t(sparse_diag_s %*% transposed_model_list_linear_constrained_A) %*% y[4:(4 + d - 1)]) + as.numeric(model_list_linear_constrained_A %*% m_adaptive) + model_list$region_lin_root_list$B
        
      }
      
      c(
        y[2] - y[3], # momentum refresh
        # diag(s_adaptive) %*% model_list$region_lin_root_list$A %*% y[4:(4 + d - 1)] + model_list$region_lin_root_list$A %*% m_adaptive + model_list$region_lin_root_list$B 
        lin_constraints_root_fun # root function to detect crossing boundary of discontinuous gradient
      )
      
    } else { # if additional linear root functions are given, similar as above, but also an extra set of additional linear root functions
      
      if (num_of_lin_constraints == 0) {
        
        lin_constraints_root_fun <- NULL
        
      } else if (num_of_lin_constraints == 1) {
        
        lin_constraints_root_fun <- c(s_adaptive * transposed_model_list_linear_constrained_A) %*% y[4:(4 + d - 1)] + model_list$region_lin_root_list$A %*% m_adaptive + model_list$region_lin_root_list$B
        
      } else {
        
        # lin_constraints_root_fun <- t(diag(s_adaptive) %*% transposed_model_list_linear_constrained_A) %*% y[4:(4 + d - 1)] + model_list$region_lin_root_list$A %*% m_adaptive + model_list$region_lin_root_list$B
        lin_constraints_root_fun <- as.numeric(Matrix::t(sparse_diag_s %*% transposed_model_list_linear_constrained_A) %*% y[4:(4 + d - 1)]) + as.numeric(model_list_linear_constrained_A %*% m_adaptive) + model_list$region_lin_root_list$B
        
      }
      
      if (num_of_additional_lin_eqs == 1) {
        
        additional_lin_root_fun <- c(s_adaptive * transposed_additional_lin_root_A) %*% y[4:(4 + d - 1)] + additional_lin_root_list$A %*% m_adaptive + additional_lin_root_list$B
        
      } else {
        
        # lin_constraints_root_fun <- t(diag(s_adaptive) %*% transposed_model_list_linear_constrained_A) %*% y[4:(4 + d - 1)] + model_list$region_lin_root_list$A %*% m_adaptive + model_list$region_lin_root_list$B
        additional_lin_root_fun <- as.numeric(Matrix::t(sparse_diag_s %*% transposed_additional_lin_root_A) %*% y[4:(4 + d - 1)]) + as.numeric(additional_lin_root_A %*% m_adaptive) + additional_lin_root_list$B
        
      }
      
      c(
        y[2] - y[3],
        # diag(s_adaptive) %*% model_list$region_lin_root_list$A %*% y[4:(4 + d - 1)] + model_list$region_lin_root_list$A %*% m_adaptive + model_list$region_lin_root_list$B 
        lin_constraints_root_fun,
        additional_lin_root_fun
      )
      
    }
    
  }
  
  region_lin_event_fun <- function(t, y, parms, lin.root.func) {
    
    yroot <- lin.root.func(t, y, parms)
    
    whichroot <- which(min(abs(yroot)) == abs(yroot))
    
    if (t > 0) {
      
      if (whichroot == 1) { # if root due to first root function --> momentum refresh
        
        if (verbose_at_refresh) {
          
          print(paste0("t: ", t))
          
        }
        
        if (t > proportion_time_until_adaptive_start * T) { # adaptive period --> update s, m and lambda
          
          old_m_adaptive <<- m_adaptive
          old_s_adaptive <<- s_adaptive
          
          m_adaptive <<- y[(4 + 4 * d):(4 + 5 * d - 1)] / t
          proposed_s_adaptive <-
            sqrt(1 / (y[(4 + 6 * d):(4 + 7 * d - 1)] / t))
          
          s_adaptive <<- pmin(max_proportion_of_previous_state * s_adaptive, pmax(min_proportion_of_previous_state * s_adaptive, proposed_s_adaptive))
          sparse_diag_s <<- Matrix::sparseMatrix(i = 1:d, j = 1:d, x = s_adaptive)
          
          if (!nut_time_found_indicator) { # for this part, the nut time is censored as one has not observed such an event
            
            sum_all_nut_times <<- sum_all_nut_times + (t - t_at_current_qbar)
            
          } else {
            
            nut_time_found_indicator <<- FALSE # nut event has happened before momentum refresh --> refresh to enable root related to nut again
            
          }
          
          # reflection_indicator <<- FALSE
          
          current_qbar <<- (1 / s_adaptive) * old_s_adaptive * y[4:(4 + d - 1)] + (1 / s_adaptive) * (old_m_adaptive - m_adaptive)
          t_at_current_qbar <<- t
          lambda_adaptive <<- n_uncensored_nut_times / sum_all_nut_times
          
          if (lambda_adaptive < lambda_lower_limit) {
            lambda_adaptive <<- lambda_lower_limit
          }
          
          new.state <- c(
            y[1] + 1,
            0,
            rexp(1),
            current_qbar, # qbar needs to be changed in order for q to stay still after an event where m and S have both changed
            rnorm(d),
            y[(4 + 2 * d):(4 + 7 * d - 1)]
          )
          
          
        } else { # safeguard period at the beginning, no adaptive yet
          
          if (!nut_time_found_indicator) {
            
            sum_all_nut_times <<- sum_all_nut_times + (t - t_at_current_qbar)
            
          } else {
            
            nut_time_found_indicator <<- FALSE
            
          }
          
          t_at_current_qbar <<- t
          current_qbar <<- y[4:(4 + d - 1)]
          
          # reflection_indicator <<- FALSE
          
          new.state <- c(
            y[1] + 1,
            0,
            rexp(1),
            y[4:(4 + d - 1)],
            rnorm(d),
            y[(4 + 2 * d):(4 + 7 * d - 1)]
          )
          
        }
        
      } else if (whichroot > 1 & whichroot <= (num_of_lin_constraints + 1)) { # due to crossing boundary of discontinuous gradient
        
        qbar <- y[4:(4 + d - 1)]
        pbar <- y[(4 + d):(4 + 2 * d - 1)] 
        # position slightly after boundary is hit (to be improved in code...): 
        qp <- qbar + 0.0001 * pbar 
        model_list$grad_jump_fun(m_adaptive + s_adaptive * qp)
        # q and p unchanged, only id is updated
        # new.state <- y
        # new.state[4 + 7 * d] <- new_region_id
        new.state <- c(y[1:3], qbar, pbar, y[(4 + 2 * d):(4 + 7 * d - 1)])
        
      } else { # due to the additional linear root functions that are given
        
        new.state <- additional_lin_root_list$event_fun(t, y, parms, m_vector = m_adaptive, s_vector = s_adaptive, adaptive_yes_or_no = T)
        
      }
      
      
    } else {
      
      new.state <- y
      
    }
    
    new.state
    
  }
  
  y0 <- c(
    0,
    0,
    u_initial,
    qbar_initial,
    pbar_initial,
    rep(0, d),
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
    paste0("int_q", 1:d, "_squared"),
    paste0("int_original_grad_log_target_squared_", 1:d)
  )
  
  q_original_samples <- t(m_adaptive + diag(s_adaptive, nrow = d) %*% t(df_sim_out[, 4:(4 + d - 1) + 1]))[-1, ]
  
  return(
    list(
      output_from_ode_solver = df_sim_out, 
      n_evals_ode = sim_out$n.evals, 
      m_adaptive = as.numeric(m_adaptive), 
      s_adaptive = as.numeric(s_adaptive),
      lambda_adaptive = lambda_adaptive,
      sum_all_nut_times = sum_all_nut_times,
      n_uncensored_nut_times = n_uncensored_nut_times,
      qbar_end = as.numeric(df_sim_out[nrow(df_sim_out), 4:(4 + d - 1) + 1]), # plus one since lsodar adds a time column as well
      pbar_end = as.numeric(df_sim_out[nrow(df_sim_out), (4 + d):(4 + 2 * d - 1) + 1]),
      u_end = df_sim_out[nrow(df_sim_out), 3 + 1], # value of u at the time the adaptive process ends
      Lambda_end = df_sim_out[nrow(df_sim_out), 2 + 1] # value of Lambda at the time the adaptive process ends
    )
  )
  
}
