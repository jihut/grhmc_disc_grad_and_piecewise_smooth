rm(list = ls())

library(doParallel)
library(dplyr)
library(ggplot2)
source("disc_grad/implementation_scripts/ISG/full_ISG_grhmc_discontinuous_gradient_transformed_function.R")

# For a range of integrator tolerance (absolute and relative)
# Run 100 independent trajectories and look at the samples generated all together. 
# Compare marginal of q_2 to verify that everything is correct. 

true_cdf_q2 <- function(x) {
  
  0.5 * (pnorm(x) + pnorm(sqrt(2) / 2 * x) ^ 2)
  
}

true_pdf_q2 <- function(q2) {
  
  0.5 * (
    dnorm(q2) + 
      sqrt(2) * pnorm(sqrt(2) * q2 / 2) * dnorm(sqrt(2) * q2 / 2)
  )
  
}

n_dim_vec <- 2 ^ {1:11}
n_dim_vec <- sort(n_dim_vec, decreasing = T)
n_cores <- 10
n_trajectories <- 100
trajectory_length <- 10000
n_samples <- 10000
init_cluster <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(init_cluster)

for (k in 1:length(n_dim_vec)) {
  print(k)
  n_dim <- n_dim_vec[k] # choose another value here if wanted
  
  doRNG::registerDoRNG(seed = 42)
  
  store_chains <- foreach::foreach(i = 1:n_trajectories) %dopar% {
    dir.create("disc_grad/example_scripts/toy_model/grhmc/log", showWarnings = F)
    sink(
      paste0(
        "disc_grad/example_scripts/toy_model/grhmc/log/log_nr",
        i,
        ".txt"
      )
    )
    n_dim <- n_dim
    
    max_model <- list(
      
      n_parameters = n_dim,
      
      grad_jump_fun = function(q) {
        
        grad_indices <<- c(
          as.integer(q[1] > 0)
        )
        
      },
      
      additional_non_lin_root_list = NULL, # usually: list(root_fun = ..., event_fun = ...)
      
      additional_lin_root_list = NULL, # usually: list(A = ..., B = ..., event_fun = ...)
      
      region_lin_root_list = list(
        A = matrix(c(1, rep(0, n_dim - 1)), nrow = 1), 
        B = 0
      ), # also need to specify the break point here again
      
      sim_q0 = function() {
        q <- numeric(n_dim)
        q[1] <- rnorm(1)
        q[2:n_dim] <- rnorm(n_dim - 1, mean = max(0, q[1]), sd = 1)
        q
      }
      
    )
    
    max_model$log_target_grad <- function(q) { # define the log target grad
      
      c(
        -q[1] + sum((q[2:n_dim] - q[1]) * grad_indices),
        -q[2:n_dim] + q[1] * grad_indices
      )
      
    }
    
    model_list <- max_model
    qbar_initial <- model_list$sim_q0()
    pbar_initial <- rnorm(n_dim)
    u_initial <- rexp(1)
    print("warmup")
    test_run_warmup <- grhmc_discontinuous_gradient_transformed_function(
      model_list = model_list,
      lambda = 0.2,
      T = trajectory_length,
      n_samples = n_samples / 20,
      diag_s_elements_initial = rep(1, n_dim),
      m_initial = rep(0, n_dim),
      qbar_initial = qbar_initial,
      pbar_initial = pbar_initial,
      Lambda_initial = 0,
      u_initial = u_initial,
      random_state = NULL,
      rtol = 1e-4,
      atol = 1e-4,
      verbose_at_refresh = T,
      return_output_from_ode = F,
      relevant_indices = 1:10 # only keep the first ten coordinates as safety, the remaining should be the same as q2 anyways
    )
    
    new_pbar_initial <- rnorm(n_dim)
    new_u_initial <- rexp(1)
    
    print("sampling")
    test_run_sampling <- grhmc_discontinuous_gradient_transformed_function(
      model_list = model_list,
      lambda = 0.2,
      T = trajectory_length,
      n_samples = n_samples,
      diag_s_elements_initial = rep(1, n_dim),
      m_initial = rep(0, n_dim),
      qbar_initial = test_run_warmup$qbar_final,
      pbar_initial = new_pbar_initial,
      Lambda_initial = 0,
      u_initial = new_u_initial,
      random_state = NULL,
      rtol = 1e-4,
      atol = 1e-4,
      verbose_at_refresh = T,
      return_output_from_ode = F,
      relevant_indices = 1:10
    )
    sink()
    test_run_sampling
    
  }
  
  saveRDS(store_chains, paste0("disc_grad/example_scripts/toy_model/grhmc/max_model_d_", n_dim, "_grhmc.RDS"))
  rm(store_chains)
  gc()
}

parallel::stopCluster(init_cluster)
