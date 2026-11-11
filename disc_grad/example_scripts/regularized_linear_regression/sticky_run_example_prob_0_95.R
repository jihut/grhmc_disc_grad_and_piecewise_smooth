# Run the standard spike and slab variable selection strategy
# For comparison with GRHMC setup: For each variable, the prior probability of beta_i being zero is set to 0.95 with the variance of the slab being 1. 

rm(list = ls())

library(MASS)
# library(brms)
library(ggplot2)
library(dplyr)
library(rstan)
library(foreach)
library(gridExtra)

source("sticky/implementation_scripts/general_scripts/grhmc_sticky_transformed_function.R")

n_cores <- 10
n_iterations <- 10
n_samples_per_iteration <- 1000000
time_period <- 50000

head(Boston)

x <- as.matrix(Boston[, -ncol(Boston)]) # 13 covariates
y <- Boston$medv
n_rows <- length(y)
length_beta <- ncol(x)

scale_x <- scale(x)
scale_y <- scale(y)
mean_y <- mean(y)

Boston_scale <- data.frame(scale_x, medv = scale_y)

init_cluster <- parallel::makeCluster(n_cores)

doParallel::registerDoParallel(init_cluster)

doRNG::registerDoRNG(seed = 42)

store_runs <- foreach::foreach(l = 1:n_iterations) %dopar% {
  
  og_regmod <- lm(medv ~ ., data = Boston_scale)
  
  sink(paste0("disc_grad/example_scripts/regularized_linear_regression/sticky_log_folder/log_nr", l, ".txt"))
  
  print("Start run")
  
  # Also fit using scaled x and y
  
  length_beta <- ncol(x)
  
  # beta_j = max(0, beta_j+) - max(0, beta_j-) with beta_j- and beta_j+ having N(mu, rho^2) as prior
  
  # Also define gamma = log (sigma ^ 2) and set a standard normal distribution prior on gamma as well
  # Default prior for sigma in stan_glm is Exp(k) where k is the rate parameter, this can be used here as well
  k <- 1 
  
  grhmc_model_list <- list(
    
    n_parameters = length_beta + 1,
    
    sticky_indices = 1:length_beta,
    
    sim_q0 = function() rnorm(length_beta + 1)
    
  )
  prior_prob_zero <- 0.95
  prior_var_slab <- 1
  grhmc_model_list$weight_slab <- rep(1 - prior_prob_zero, (grhmc_model_list$n_parameters) - 1) # 95% prior probability of beta_i = 0
  grhmc_model_list$sigma_slab <- rep(sqrt(prior_var_slab), grhmc_model_list$n_parameters - 1)
  
  grhmc_model_list$log_target_grad <- function(q) { # new version, vectorized and therefore slightly faster
    
    beta <- q[1:length_beta]
    y_minus_linear_predictor <- as.numeric(scale_y - (scale_x %*% beta))
    exp_gamma <- exp(q[grhmc_model_list$n_parameters])
    log_target_grad_value_vec <- numeric(grhmc_model_list$n_parameters)
    
    prior_beta_term <- -beta / grhmc_model_list$sigma_slab ^ 2
    
    y_minus_linear_predictor_times_scale_x <- y_minus_linear_predictor %*% scale_x
    
    # If grad_indices = 0 for i = 1:length_beta --> beta_i+ < 0 so only the prior part in this region, no term related to observations and covariates
    log_target_grad_value_vec[1:length_beta] <- prior_beta_term + y_minus_linear_predictor_times_scale_x / exp_gamma
    log_target_grad_value_vec[grhmc_model_list$n_parameters] <- - (n_rows / 2) - q[grhmc_model_list$n_parameters] + (sum(y_minus_linear_predictor ^ 2) / (2 * exp_gamma)) # if prior of gamma = log(sigma ^ 2) is N(0, 1)
    
    log_target_grad_value_vec
    
  }
  
  qbar_initial <- c(og_regmod$coefficients[-1], sum((og_regmod$residuals)^2) / og_regmod$df.residual)
  pbar_initial <- rnorm(grhmc_model_list$n_parameters)
  
  print("warmup")
  boston_scale_fit_grhmc_warmup <- grhmc_sticky_transformed_function(
    model_list = grhmc_model_list,
    lambda = 2.5,
    T = time_period,
    n_samples = n_samples_per_iteration,
    diag_s_elements_initial = rep(1, grhmc_model_list$n_parameters),
    m_initial = rep(0, grhmc_model_list$n_parameters),
    qbar_initial = qbar_initial,
    pbar_initial = pbar_initial,
    Lambda_initial = 0,
    u_initial = rexp(1),
    random_state = NULL,
    rtol = 1e-4,
    atol = 1e-4,
    h_max = 1.0,
    last.root.offset.lin.root.finder = 1.0e-8,
    last.root.offset.non.lin.root.finder = 1.0e-8,
    precision_real_root_lin_root_finder = 1.0e-13,
    num_subdiv_non_lin_root_finder = 8L,
    verbose_at_refresh = T
  ) 
  
  new_qbar_initial <- tail(boston_scale_fit_grhmc_warmup$q_original_samples, 1)
  new_qbar_initial[new_qbar_initial == 0] <- rnorm(sum(new_qbar_initial == 0), 0, 0.01)
  new_pbar_initial <- rnorm(grhmc_model_list$n_parameters)
  
  print("sampling")
  boston_scale_fit_grhmc_sampling <- grhmc_sticky_transformed_function(
    model_list = grhmc_model_list,
    lambda = 2.5,
    T = time_period,
    n_samples = n_samples_per_iteration,
    diag_s_elements_initial = rep(1, grhmc_model_list$n_parameters),
    m_initial = rep(0, grhmc_model_list$n_parameters),
    qbar_initial = new_qbar_initial,
    pbar_initial = new_pbar_initial,
    Lambda_initial = 0,
    u_initial = rexp(1),
    random_state = NULL,
    rtol = 1e-4,
    atol = 1e-4,
    h_max = 1.0,
    last.root.offset.lin.root.finder = 1.0e-8,
    last.root.offset.non.lin.root.finder = 1.0e-8,
    precision_real_root_lin_root_finder = 1.0e-13,
    num_subdiv_non_lin_root_finder = 8L,
    verbose_at_refresh = T
  ) 
  
  print("End run")
  
  sink()
  
  list(
    q_original_samples = boston_scale_fit_grhmc_sampling$q_original_samples,
    n_evals_ode = boston_scale_fit_grhmc_sampling$n_evals_ode
  )
  
}

saveRDS(store_runs, "disc_grad/example_scripts/regularized_linear_regression/sticky_run_example_prob_0_95.RDS")

parallel::stopCluster(init_cluster)
