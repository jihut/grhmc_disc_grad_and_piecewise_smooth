source("sticky/implementation_scripts/general_scripts/deSolverRoot.R")

grhmc_sticky_transformed_function <- function(
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
  unfreeze_indices <- rep(1, d) 
  length_sticky_indices <- length(model_list$sticky_indices)
  
  # define the (first order) ode for state
  # state[1]: number of events up to time point
  # state[2]: Lambda, which resets after each event
  # state[3]: u, simulate from exp(1) after each event
  # state[4:(4 + d - 1)]: qbar
  # state[(4 + d):(4 + 2*d - 1)]: pbar
  # state[(4 + 2 * d):(4 + 3 * d - 2)]: Given that a component is freeze, how much time before unfreezes?
  
  ode <- function(t, state, parms){
    
    count_n_evals_ode()
    ret <- c(
      
      0, # number of momentum refresh events
      lambda, # integrate lambda to get Lambda
      0, # u
      state[(4 + d):(4 + 2 * d - 1)] * unfreeze_indices, # \dot qbar = pbar
      diag_s_elements_initial *
        model_list$log_target_grad(m_initial + diag_s_elements_initial * state[4:(4 + d - 1)]) * 
        unfreeze_indices, #\dot pbar = gradient of log transformed density wrt. qbar
      rep(-1, length_sticky_indices) # countdown how much time before unfreeze
    )
    
    ret
    
  }
  
  region_lin_root_fun <- function(t, y, parms) {
    c(
      y[2] - y[3],
      diag_s_elements_initial[model_list$sticky_indices] * y[4:(4 + d - 1)][model_list$sticky_indices] + 
        m_initial[model_list$sticky_indices] + (1 - unfreeze_indices[model_list$sticky_indices]) * 1e+10, # boundary hitting
      y[(4 + 2 * d):(4 +  2 * d + length_sticky_indices - 1)] # freeze time elapsed done, now need to unfreeze
    )
  }
  
  region_lin_event_fun <- function(t, y, parms, region_lin_root_fun) {
    
    yroot <- region_lin_root_fun(t, y, parms)
    # print(t)
    # print(yroot)
    whichroot <- which(min(abs(yroot)) == abs(yroot))
    
    if (whichroot == 1) {
      
      if (verbose_at_refresh) {
        
        print(paste0("t: ", t))
        
      }
      
      # new.state <- y
      # new.state[1] <- new.state[1] + 1
      # new.state[2] <- 0
      # new.state[3] <- rexp(1)
      # new.state[(4 + d):(4 + 2 * d - 1)] <- rnorm(d)
      
      pbar <- y[(4 + d):(4 + 2 * d - 1)]
      if (sum(unfreeze_indices) > 0) {
        pbar[unfreeze_indices == 1] <- rnorm(sum(unfreeze_indices), 0, 1) 
      }
      # HERE: ONLY UPDATE MOMENTUM OF UNFREEZED COORDINATES
      
      new.state <- c(
        y[1] + 1,
        0,
        rexp(1),
        y[4:(4 + d - 1)],
        pbar,
        y[(4 + 2 * d):(4 + 2 * d + length_sticky_indices - 1)]
      )
      
    } else if (whichroot > 1 & whichroot <= length_sticky_indices + 1) {
      
      qbar <- y[4:(4 + d - 1)]
      pbar <- y[(4 + d):(4 + 2 * d - 1)]
      time_unfreeze <- y[(4 + 2 * d):(4 + 2 * d + length_sticky_indices - 1)]
      # print("####")
      # print(t)
      
      relevant_index <- whichroot - 1
      index_component_freeze <- model_list$sticky_indices[whichroot - 1]
      # print(index_component_freeze)
      unfreeze_indices[index_component_freeze] <<- 0
      kappa <- (model_list$weight_slab[relevant_index]) / (1 - model_list$weight_slab[relevant_index]) * 
        dnorm(0, sd = model_list$sigma_slab[relevant_index])
      # print(kappa)
      exp_rate <- kappa * abs(pbar[index_component_freeze] * diag_s_elements_initial[index_component_freeze])
      # print(exp_rate)
      sim_unfreeze_time <- rexp(1, rate = exp_rate)
      time_unfreeze[relevant_index] <- sim_unfreeze_time
      qbar[index_component_freeze] <- 0
      new.state <- c(
        y[1:3],
        qbar,
        pbar, 
        time_unfreeze
      )
      
    } else if (whichroot > length_sticky_indices + 1 & whichroot <= (2 * length_sticky_indices + 1)) {
      
      relevant_index <- whichroot - (length_sticky_indices + 1)
      index_component_unfreeze <- model_list$sticky_indices[relevant_index]
      unfreeze_indices[index_component_unfreeze] <<- 1
      qbar <- y[4:(4 + d - 1)]
      qbar[index_component_unfreeze] <- y[(4 + d):(4 + 2 * d - 1)][index_component_unfreeze] * 1e-10
      
      new.state <- c(
        y[1:3],
        qbar, 
        y[(4 + d):(4 + 2 * d - 1)],
        y[(4 + 2 * d):(4 + 2 * d + length_sticky_indices - 1)]
      )
      
    } else {
      
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
    rep(-1, length_sticky_indices)
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
    root.func = NULL,
    event.func = NULL,
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
    paste0("time_until_unfreeze", model_list$sticky_indices)
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
