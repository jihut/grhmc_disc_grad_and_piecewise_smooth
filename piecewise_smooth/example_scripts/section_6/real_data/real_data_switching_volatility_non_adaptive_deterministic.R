rm(list = ls())
library(doParallel)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
source("piecewise_smooth/implementation_scripts/general_scripts/grhmc_piecewise_smooth_density_transformed_function.R")

data <- read.table("piecewise_smooth/example_scripts/section_6/real_data/data.dat")
data

y_vec <- data$V1

num_time_points <- length(y_vec)

init_cluster <- parallel::makeCluster(5)
doParallel::registerDoParallel(init_cluster)

# Split up into two parts to avoid memory issues

# First: 5 independent trajectories where z_initial is from a normal distribution with mean -1

store_matrix_nr1 <- foreach::foreach(l = 1:5, .combine = "rbind") %dopar% {
  # sink(paste0("piecewise_smooth/example_scripts/section_6/real_data/log/log_nr", l, ".txt"))
  y_vec <- data$V1
  
  num_time_points <- length(y_vec)
  
  prior_param_gamma_low <- 1 # exponential prior on sigma here
  prior_param_gamma_high <- 0.5
  
  prior_shape1_rho <- 2 # prior parameters when rho prior is transformed beta
  prior_shape2_rho <- 2
  
  # q[1:num_time_points] <- z(1):z(num_time_points)
  # q[num_time_points + 1] <- rho_star (with rho = tanh(rho_star))
  # q[num_time_points + 2] <- gamma_low (with sigma_low = exp(0.5 * gamma_low))
  # q[num_time_points + 3] <- gamma_high (with sigma_high = exp(0.5 * gamma_high))
  
  n_parameters <- num_time_points + 3
  
  A <- diag(num_time_points)
  A <- cbind(A, matrix(0, nrow = num_time_points, ncol = 3))
  dim(A)
  
  B <- rep(0, num_time_points)
  
  grhmc_model_list <- list(
    
    n_parameters = n_parameters,
    
    target_jump_fun = function(q) {
      
      region_id <<- sapply(1:num_time_points, function(i) as.integer(q[i] > 0))
      
    },
    
    region_lin_root_list = list(
      A = A, 
      B = B
    )
    
  )
  
  additional_A <- matrix(c(rep(0, num_time_points + 1), -1, 1), nrow = 1) # impose gamma_high > gamma_low
  additional_B <- 0
  
  additional_event_fun <- function(t, y, parms, m_vector, s_vector, adaptive_yes_or_no) { # reflection based on the constraint above
    
    d <- grhmc_model_list$n_parameters
    
    old_qbar <- y[4:(4 + d - 1)]
    old_pbar <- y[(4 + d):(4 + 2 * d - 1)]
    
    normal_vec <- s_vector * additional_A
    
    # Next: Find the projection of pbar onto the normal vec
    old_pbar_perpendicular <- sum(old_pbar * normal_vec) / sum(normal_vec ^ 2) * normal_vec
    # print(paste0("pbar_perpendicular: ", old_pbar_perpendicular))
    old_pbar_parallel <- old_pbar - old_pbar_perpendicular
    
    new_pbar_perpendicular <- -old_pbar_perpendicular # no randomized reflection for now
    
    new_pbar <- old_pbar_parallel + new_pbar_perpendicular
    
    if (adaptive_yes_or_no) {
      
      if (!nut_time_found_indicator) {
        
        sum_all_nut_times <<- sum_all_nut_times + (t - t_at_current_qbar)
        
      } else {
        
        nut_time_found_indicator <<- FALSE
        
      }
      
      t_at_current_qbar <<- t
      current_qbar <<- old_qbar
      
      if (t > proportion_time_until_adaptive_start * T) {
        
        lambda_adaptive <<- n_uncensored_nut_times / sum_all_nut_times
        
        if (lambda_adaptive < lambda_lower_limit) {
          lambda_adaptive <<- lambda_lower_limit
        }
        
      }
      
      if (sampling_compute_temporal_averages_of_moments) {
        
        new.state <- c(
          y[1:3],
          old_qbar,
          # new_qbar,
          new_pbar,
          y[(4 + 2 * d):(4 + 7 * d - 1)]
        ) 
        
      } else {
        
        new.state <- c(
          y[1:3],
          old_qbar,
          # new_qbar,
          new_pbar,
          y[(4 + 2 * d):(4 + 3 * d - 1)]
        )
        
      }
      
    } else {
      
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
    
    print("reflection as gamma_low = gamma_high")
    
    new.state
    
  }
  
  grhmc_model_list$additional_lin_root_list <- list(A = additional_A, B = additional_B, event_fun = additional_event_fun)
  
  grhmc_model_list$log_target_grad <- function(q) {
    
    time_indices_2_to_t_minus_1 <- 2:(num_time_points - 1)
    time_indices_3_to_t <- 3:(num_time_points)
    time_indices_1_to_t_minus_2 <- 1:(num_time_points - 2)
    time_indices_2_to_t <- c(time_indices_2_to_t_minus_1, num_time_points)
    time_indices_1_to_t_minus_1 <- c(time_indices_1_to_t_minus_2, num_time_points - 1)
    
    
    rho_star <- q[num_time_points + 1]
    rho <- tanh(rho_star)
    cosh_rho_star <- cosh(rho_star)
    sinh_rho_star <- sinh(rho_star)
    
    gamma_low <- q[num_time_points + 2]
    sigma_low <- exp(gamma_low / 2)
    
    gamma_high <- q[num_time_points + 3]
    sigma_high <- exp(gamma_high / 2)
    
    sigma <- sigma_high * region_id + sigma_low * (1 - region_id)
    
    z_t <- q[1:num_time_points]
    
    log_target_grad_value_vec <- numeric(grhmc_model_list$n_parameters)
    
    log_target_grad_value_vec[1] <-
      (-rho * y_vec[1] / (sigma[1]) + 
         rho * y_vec[2] / (sigma[2]) + 
         2 * z_t[1] - z_t[2]) / (rho ^ 2 - 1)
    
    log_target_grad_value_vec[time_indices_2_to_t_minus_1] <-
      (-rho * y_vec[time_indices_2_to_t_minus_1] / (sigma[time_indices_2_to_t_minus_1]) + 
         rho * y_vec[time_indices_3_to_t] / (sigma[time_indices_3_to_t]) -
         z_t[time_indices_1_to_t_minus_2] + 2 * z_t[time_indices_2_to_t_minus_1] - z_t[time_indices_3_to_t]) / (rho ^ 2 - 1)
    
    
    log_target_grad_value_vec[num_time_points] <- 
      (-rho * y_vec[num_time_points] / (sigma[num_time_points]) -
         z_t[num_time_points - 1] + z_t[num_time_points]) / (rho ^ 2 - 1)
    
    log_target_grad_value_vec[num_time_points + 1] <- 
      num_time_points * rho + 
      (-z_t[1] / sigma[1] + rho * y_vec[1] / sigma[1] ^ 2) * (y_vec[1] + rho * sigma[1] * (-z_t[1])) / (rho ^ 2 - 1) +
      sum(
        ((-z_t[time_indices_2_to_t] + z_t[time_indices_1_to_t_minus_1]) / sigma[time_indices_2_to_t] + rho * y_vec[time_indices_2_to_t] / sigma[time_indices_2_to_t] ^ 2) * (y_vec[time_indices_2_to_t] + rho * sigma[time_indices_2_to_t] * (-z_t[time_indices_2_to_t] + z_t[time_indices_1_to_t_minus_1])) / (rho ^ 2 - 1)
      ) + 
      (-prior_shape1_rho - prior_shape2_rho) * tanh(rho_star) + prior_shape1_rho - prior_shape2_rho # prior term if prior of rho is transformed beta
    
    log_target_grad_value_vec[num_time_points + 2] <- 
      (-1 / 2 + (cosh_rho_star ^ 2 * y_vec[1] ^ 2) / (2 * sigma_low ^ 2) - cosh_rho_star * sinh_rho_star * y_vec[1] * z_t[1] / (2 * sigma_low)) * (1 - region_id[1]) + 
      sum(
        (-1 / 2 + (cosh_rho_star ^ 2 * y_vec[time_indices_2_to_t] ^ 2) / (2 * sigma_low ^ 2) - cosh_rho_star * sinh_rho_star * y_vec[time_indices_2_to_t] * (z_t[time_indices_2_to_t] - z_t[time_indices_1_to_t_minus_1]) / (2 * sigma_low)) * (1 - region_id[time_indices_2_to_t])
      ) + 
      -prior_param_gamma_low / 2 * exp(0.5 * gamma_low) + 1 / 2 # prior term if prior of sigma_low is exponential
    
    log_target_grad_value_vec[num_time_points + 3] <- 
      (-1 / 2 + (cosh_rho_star ^ 2 * y_vec[1] ^ 2) / (2 * sigma_high ^ 2) - cosh_rho_star * sinh_rho_star * y_vec[1] * z_t[1] / (2 * sigma_high)) * (region_id[1]) + 
      sum(
        (-1 / 2 + (cosh_rho_star ^ 2 * y_vec[time_indices_2_to_t] ^ 2) / (2 * sigma_high ^ 2) - cosh_rho_star * sinh_rho_star * y_vec[time_indices_2_to_t] * (z_t[time_indices_2_to_t] - z_t[time_indices_1_to_t_minus_1]) / (2 * sigma_high)) * (region_id[time_indices_2_to_t])
      ) + 
      -prior_param_gamma_high / 2 * exp(0.5 * gamma_high) + 1 / 2 # prior term if prior of sigma_high is exponential
    
    
    log_target_grad_value_vec
    
  }
  
  grhmc_model_list$log_target_fun <- function(q) {
    time_indices_2_to_t <- 2:num_time_points
    time_indices_1_to_t_minus_1 <- time_indices_2_to_t - 1
    rho_star <- q[num_time_points + 1]
    rho <- tanh(rho_star)
    cosh_rho_star <- cosh(rho_star)
    sinh_rho_star <- sinh(rho_star)
    
    gamma_low <- q[num_time_points + 2]
    sigma_low <- exp(gamma_low / 2)
    
    gamma_high <- q[num_time_points + 3]
    sigma_high <- exp(gamma_high / 2)
    
    sigma <- sigma_high * region_id + sigma_low * (1 - region_id)
    
    z_t <- q[1:num_time_points]
    
    log_target_value <- dnorm(z_t[1], log = T) + dnorm(y_vec[1], mean = rho * sigma[1] * z_t[1], sd = sqrt(1 - rho ^ 2) * sigma[1], log = T) + 
      sum(
        dnorm(z_t[time_indices_2_to_t], mean = z_t[time_indices_1_to_t_minus_1], log = T) + 
          dnorm(y_vec[time_indices_2_to_t], mean = rho * sigma[time_indices_2_to_t] * (z_t[time_indices_2_to_t] - z_t[time_indices_1_to_t_minus_1]), sd = sqrt(1 - rho ^ 2) * sigma[time_indices_2_to_t], log = T)
      )
    
    # log_prior_rho_star <- log(2) + 2 * rho_star - 2 * log(exp(2 * rho_star) + 1) # if prior of rho is uniform
    log_prior_rho_star <- log(gamma(prior_shape1_rho + prior_shape2_rho)) - log(gamma(prior_shape1_rho)) - log(gamma(prior_shape2_rho)) +
      (prior_shape1_rho - 1) * log((tanh(rho_star) + 1) / 2) + (prior_shape2_rho - 1) * log(1 - (tanh(rho_star) + 1) / 2) - log(2) +
      log(1 - tanh(rho_star) ^ 2) # if prior of rho is transformed beta
    
    log_prior_gamma_low <- log(prior_param_gamma_low) - log(2) - prior_param_gamma_low * exp(1 / 2 * gamma_low) + 1 / 2 * gamma_low # if prior of sigma_low is exponential
    # log_prior_gamma_low <- dnorm(gamma_low, mean = prior_mu_gamma_low, sd = prior_sigma_gamma_low, log = T) # if prior of sigma_low^2 is normal
    # log_prior_gamma_low <- -log(gamma(prior_shape_gamma_low)) - prior_shape_gamma_low * log(prior_scale_gamma_low) + prior_shape_gamma_low * gamma_low / 2 -
    #   exp(gamma_low / 2) / prior_scale_gamma_low - log(2) # if prior of sigma_low is gamma
    
    log_prior_gamma_high <- log(prior_param_gamma_high) - log(2) - prior_param_gamma_high * exp(1 / 2 * gamma_high) + 1 / 2 * gamma_high # if prior of sigma_high is exponential
    # log_prior_gamma_high <- dnorm(gamma_high, mean = prior_mu_gamma_high, sd = prior_sigma_gamma_high, log = T) # if prior of sigma_high^2 is normal
    # log_prior_gamma_high <- -log(gamma(prior_shape_gamma_high)) - prior_shape_gamma_high * log(prior_scale_gamma_high) + prior_shape_gamma_high * gamma_high / 2 -
    #   exp(gamma_high / 2) / prior_scale_gamma_high - log(2) # if prior of sigma_high is gamma
    
    log_target_value <- log_target_value + log_prior_rho_star + log_prior_gamma_low + log_prior_gamma_high
    
    log_target_value
    
  }
  
  set.seed(l)
  
  # rho_initial <- runif(1, -1, 1)
  rho_initial <- 2 * rbeta(1, shape1 = prior_shape1_rho, shape2 = prior_shape2_rho) - 1
  rho_star_initial <- 1 / 2 * log((1 + rho_initial) / (1 - rho_initial))
  
  sigma_low_initial <- rexp(1, rate = prior_param_gamma_low) # if exponential prior for sigma
  # sigma_low_initial <- rgamma(1, shape = prior_shape_gamma_low, scale = prior_scale_gamma_low) # if gamma prior for sigma
  gamma_low_initial <- log(sigma_low_initial ^ 2)
  # gamma_low_initial <- rnorm(1, mean = prior_mu_gamma_low, sd = prior_sigma_gamma_low) # if log normal prior for sigma^2
  
  sigma_high_initial <- rexp(1, rate = prior_param_gamma_high) # if exponential prior for sigma
  # sigma_high_initial <- rgamma(1, shape = prior_shape_gamma_high, scale = prior_scale_gamma_high) # if gamma prior for sigma
  gamma_high_initial <- log(sigma_high_initial ^ 2)
  # gamma_high_initial <- rnorm(1, mean = prior_mu_gamma_high, sd = prior_sigma_gamma_high) # if log normal prior for sigma^2
  
  while (gamma_high_initial < gamma_low_initial) {
    
    sigma_low_initial <- rexp(1, rate = prior_param_gamma_low) # if exponential prior for sigma
    # sigma_low_initial <- rgamma(1, shape = prior_shape_gamma_low, scale = prior_scale_gamma_low) # if gamma prior for sigma
    gamma_low_initial <- log(sigma_low_initial ^ 2)
    # gamma_low_initial <- rnorm(1, mean = prior_mu_gamma_low, sd = prior_sigma_gamma_low) # if log normal prior for sigma^2
    
    sigma_high_initial <- rexp(1, rate = prior_param_gamma_high) # if exponential prior for sigma
    # sigma_high_initial <- rgamma(1, shape = prior_shape_gamma_high, scale = prior_scale_gamma_high) # if gamma prior for sigma
    gamma_high_initial <- log(sigma_high_initial ^ 2)
    # gamma_high_initial <- rnorm(1, mean = prior_mu_gamma_high, sd = prior_sigma_gamma_high) # if log normal prior for sigma^2
    
  }
  
  # z_initial <- numeric(num_time_points)
  # z_initial[1] <- rnorm(1)
  # for (i in 2:num_time_points) {
  #   z_initial[i] <- rnorm(1, mean = z_initial[i - 1], sd = 1)
  # }
  
  if (l <= 5) {
    z_initial <- rnorm(num_time_points, mean = -1) # an assumption that the majority of the time points are in the low volatility state
  } else {
    z_initial <- rnorm(num_time_points, mean = 1) # an assumption that the majority of the time points are in the low volatility state
  }
  
  qbar_initial <- c(z_initial, rho_star_initial, gamma_low_initial, gamma_high_initial)
  pbar_initial <- rnorm(grhmc_model_list$n_parameters)
  
  region_id <<- grhmc_model_list$target_jump_fun(qbar_initial)
  
  system.time(
    test_run <- grhmc_piecewise_smooth_density_transformed_function(
      model_list = grhmc_model_list,
      lambda = 0.2,
      T = 100000,
      n_samples = 100000,
      diag_s_elements_initial = rep(1, grhmc_model_list$n_parameters),
      m_initial = rep(0, grhmc_model_list$n_parameters),
      qbar_initial = qbar_initial,
      pbar_initial = pbar_initial,
      Lambda_initial = 0,
      u_initial = rexp(1),
      random_state = NULL,
      rtol = NULL,
      atol = NULL,
      h_max = 1.0,
      sampling_compute_temporal_averages_of_moments = FALSE,
      reflection_type = "deterministic", 
      verbose_at_refresh = T,
      last.root.offset.lin.root.finder = 1e-10,
      last.root.offset.non.lin.root.finder = 1e-10,
      precision_real_root_lin_root_finder = 1.0e-13,
      num_subdiv_non_lin_root_finder = 8L
    )
  )
  print("Finished.")
  cbind(test_run$q_original_samples, rep(l, nrow(test_run$q_original_samples)))
  
}
parallel::stopCluster(init_cluster)
# saveRDS(store_matrix_nr1, "piecewise_smooth/example_scripts/section_6/real_data/real_data_switching_volatility_non_adaptive_deterministic_samples_nr1.RDS")
store_matrix_nr1 <- readRDS("piecewise_smooth/example_scripts/section_6/real_data/real_data_switching_volatility_non_adaptive_deterministic_samples_nr1.RDS")

# Second: 5 independent trajectories where z_initial is from a normal distribution with mean 1

init_cluster <- parallel::makeCluster(5)
doParallel::registerDoParallel(init_cluster)
store_matrix_nr2 <- foreach::foreach(l = 6:10, .combine = "rbind") %dopar% {
  sink(paste0("piecewise_smooth/example_scripts/section_6/real_data/log/log_nr", l, ".txt"))
  y_vec <- data$V1
  
  num_time_points <- length(y_vec)
  
  prior_param_gamma_low <- 1 # exponential prior on sigma here
  prior_param_gamma_high <- 0.5
  
  prior_shape1_rho <- 2 # prior parameters when rho prior is transformed beta
  prior_shape2_rho <- 2
  
  # q[1:num_time_points] <- z(1):z(num_time_points)
  # q[num_time_points + 1] <- rho_star (with rho = tanh(rho_star))
  # q[num_time_points + 2] <- gamma_low (with sigma_low = exp(0.5 * gamma_low))
  # q[num_time_points + 3] <- gamma_high (with sigma_high = exp(0.5 * gamma_high))
  
  n_parameters <- num_time_points + 3
  
  A <- diag(num_time_points)
  A <- cbind(A, matrix(0, nrow = num_time_points, ncol = 3))
  dim(A)
  
  B <- rep(0, num_time_points)
  
  grhmc_model_list <- list(
    
    n_parameters = n_parameters,
    
    target_jump_fun = function(q) {
      
      region_id <<- sapply(1:num_time_points, function(i) as.integer(q[i] > 0))
      
    },
    
    region_lin_root_list = list(
      A = A, 
      B = B
    )
    
  )
  
  additional_A <- matrix(c(rep(0, num_time_points + 1), -1, 1), nrow = 1) # impose gamma_high > gamma_low
  additional_B <- 0
  
  additional_event_fun <- function(t, y, parms, m_vector, s_vector, adaptive_yes_or_no) { # reflection based on the constraint above
    
    d <- grhmc_model_list$n_parameters
    
    old_qbar <- y[4:(4 + d - 1)]
    old_pbar <- y[(4 + d):(4 + 2 * d - 1)]
    
    normal_vec <- s_vector * additional_A
    
    # Next: Find the projection of pbar onto the normal vec
    old_pbar_perpendicular <- sum(old_pbar * normal_vec) / sum(normal_vec ^ 2) * normal_vec
    # print(paste0("pbar_perpendicular: ", old_pbar_perpendicular))
    old_pbar_parallel <- old_pbar - old_pbar_perpendicular
    
    new_pbar_perpendicular <- -old_pbar_perpendicular # no randomized reflection for now
    
    new_pbar <- old_pbar_parallel + new_pbar_perpendicular
    
    if (adaptive_yes_or_no) {
      
      if (!nut_time_found_indicator) {
        
        sum_all_nut_times <<- sum_all_nut_times + (t - t_at_current_qbar)
        
      } else {
        
        nut_time_found_indicator <<- FALSE
        
      }
      
      t_at_current_qbar <<- t
      current_qbar <<- old_qbar
      
      if (t > proportion_time_until_adaptive_start * T) {
        
        lambda_adaptive <<- n_uncensored_nut_times / sum_all_nut_times
        
        if (lambda_adaptive < lambda_lower_limit) {
          lambda_adaptive <<- lambda_lower_limit
        }
        
      }
      
      if (sampling_compute_temporal_averages_of_moments) {
        
        new.state <- c(
          y[1:3],
          old_qbar,
          # new_qbar,
          new_pbar,
          y[(4 + 2 * d):(4 + 7 * d - 1)]
        ) 
        
      } else {
        
        new.state <- c(
          y[1:3],
          old_qbar,
          # new_qbar,
          new_pbar,
          y[(4 + 2 * d):(4 + 3 * d - 1)]
        )
        
      }
      
    } else {
      
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
    
    print("reflection as gamma_low = gamma_high")
    
    new.state
    
  }
  
  grhmc_model_list$additional_lin_root_list <- list(A = additional_A, B = additional_B, event_fun = additional_event_fun)
  
  grhmc_model_list$log_target_grad <- function(q) {
    
    time_indices_2_to_t_minus_1 <- 2:(num_time_points - 1)
    time_indices_3_to_t <- 3:(num_time_points)
    time_indices_1_to_t_minus_2 <- 1:(num_time_points - 2)
    time_indices_2_to_t <- c(time_indices_2_to_t_minus_1, num_time_points)
    time_indices_1_to_t_minus_1 <- c(time_indices_1_to_t_minus_2, num_time_points - 1)
    
    
    rho_star <- q[num_time_points + 1]
    rho <- tanh(rho_star)
    cosh_rho_star <- cosh(rho_star)
    sinh_rho_star <- sinh(rho_star)
    
    gamma_low <- q[num_time_points + 2]
    sigma_low <- exp(gamma_low / 2)
    
    gamma_high <- q[num_time_points + 3]
    sigma_high <- exp(gamma_high / 2)
    
    sigma <- sigma_high * region_id + sigma_low * (1 - region_id)
    
    z_t <- q[1:num_time_points]
    
    log_target_grad_value_vec <- numeric(grhmc_model_list$n_parameters)
    
    log_target_grad_value_vec[1] <-
      (-rho * y_vec[1] / (sigma[1]) + 
         rho * y_vec[2] / (sigma[2]) + 
         2 * z_t[1] - z_t[2]) / (rho ^ 2 - 1)
    
    log_target_grad_value_vec[time_indices_2_to_t_minus_1] <-
      (-rho * y_vec[time_indices_2_to_t_minus_1] / (sigma[time_indices_2_to_t_minus_1]) + 
         rho * y_vec[time_indices_3_to_t] / (sigma[time_indices_3_to_t]) -
         z_t[time_indices_1_to_t_minus_2] + 2 * z_t[time_indices_2_to_t_minus_1] - z_t[time_indices_3_to_t]) / (rho ^ 2 - 1)
    
    
    log_target_grad_value_vec[num_time_points] <- 
      (-rho * y_vec[num_time_points] / (sigma[num_time_points]) -
         z_t[num_time_points - 1] + z_t[num_time_points]) / (rho ^ 2 - 1)
    
    log_target_grad_value_vec[num_time_points + 1] <- 
      num_time_points * rho + 
      (-z_t[1] / sigma[1] + rho * y_vec[1] / sigma[1] ^ 2) * (y_vec[1] + rho * sigma[1] * (-z_t[1])) / (rho ^ 2 - 1) +
      sum(
        ((-z_t[time_indices_2_to_t] + z_t[time_indices_1_to_t_minus_1]) / sigma[time_indices_2_to_t] + rho * y_vec[time_indices_2_to_t] / sigma[time_indices_2_to_t] ^ 2) * (y_vec[time_indices_2_to_t] + rho * sigma[time_indices_2_to_t] * (-z_t[time_indices_2_to_t] + z_t[time_indices_1_to_t_minus_1])) / (rho ^ 2 - 1)
      ) + 
      (-prior_shape1_rho - prior_shape2_rho) * tanh(rho_star) + prior_shape1_rho - prior_shape2_rho # prior term if prior of rho is transformed beta
    
    log_target_grad_value_vec[num_time_points + 2] <- 
      (-1 / 2 + (cosh_rho_star ^ 2 * y_vec[1] ^ 2) / (2 * sigma_low ^ 2) - cosh_rho_star * sinh_rho_star * y_vec[1] * z_t[1] / (2 * sigma_low)) * (1 - region_id[1]) + 
      sum(
        (-1 / 2 + (cosh_rho_star ^ 2 * y_vec[time_indices_2_to_t] ^ 2) / (2 * sigma_low ^ 2) - cosh_rho_star * sinh_rho_star * y_vec[time_indices_2_to_t] * (z_t[time_indices_2_to_t] - z_t[time_indices_1_to_t_minus_1]) / (2 * sigma_low)) * (1 - region_id[time_indices_2_to_t])
      ) + 
      -prior_param_gamma_low / 2 * exp(0.5 * gamma_low) + 1 / 2 # prior term if prior of sigma_low is exponential
    
    log_target_grad_value_vec[num_time_points + 3] <- 
      (-1 / 2 + (cosh_rho_star ^ 2 * y_vec[1] ^ 2) / (2 * sigma_high ^ 2) - cosh_rho_star * sinh_rho_star * y_vec[1] * z_t[1] / (2 * sigma_high)) * (region_id[1]) + 
      sum(
        (-1 / 2 + (cosh_rho_star ^ 2 * y_vec[time_indices_2_to_t] ^ 2) / (2 * sigma_high ^ 2) - cosh_rho_star * sinh_rho_star * y_vec[time_indices_2_to_t] * (z_t[time_indices_2_to_t] - z_t[time_indices_1_to_t_minus_1]) / (2 * sigma_high)) * (region_id[time_indices_2_to_t])
      ) + 
      -prior_param_gamma_high / 2 * exp(0.5 * gamma_high) + 1 / 2 # prior term if prior of sigma_high is exponential
    
    
    log_target_grad_value_vec
    
  }
  
  grhmc_model_list$log_target_fun <- function(q) {
    time_indices_2_to_t <- 2:num_time_points
    time_indices_1_to_t_minus_1 <- time_indices_2_to_t - 1
    rho_star <- q[num_time_points + 1]
    rho <- tanh(rho_star)
    cosh_rho_star <- cosh(rho_star)
    sinh_rho_star <- sinh(rho_star)
    
    gamma_low <- q[num_time_points + 2]
    sigma_low <- exp(gamma_low / 2)
    
    gamma_high <- q[num_time_points + 3]
    sigma_high <- exp(gamma_high / 2)
    
    sigma <- sigma_high * region_id + sigma_low * (1 - region_id)
    
    z_t <- q[1:num_time_points]
    
    log_target_value <- dnorm(z_t[1], log = T) + dnorm(y_vec[1], mean = rho * sigma[1] * z_t[1], sd = sqrt(1 - rho ^ 2) * sigma[1], log = T) + 
      sum(
        dnorm(z_t[time_indices_2_to_t], mean = z_t[time_indices_1_to_t_minus_1], log = T) + 
          dnorm(y_vec[time_indices_2_to_t], mean = rho * sigma[time_indices_2_to_t] * (z_t[time_indices_2_to_t] - z_t[time_indices_1_to_t_minus_1]), sd = sqrt(1 - rho ^ 2) * sigma[time_indices_2_to_t], log = T)
      )
    
    # log_prior_rho_star <- log(2) + 2 * rho_star - 2 * log(exp(2 * rho_star) + 1) # if prior of rho is uniform
    log_prior_rho_star <- log(gamma(prior_shape1_rho + prior_shape2_rho)) - log(gamma(prior_shape1_rho)) - log(gamma(prior_shape2_rho)) +
      (prior_shape1_rho - 1) * log((tanh(rho_star) + 1) / 2) + (prior_shape2_rho - 1) * log(1 - (tanh(rho_star) + 1) / 2) - log(2) +
      log(1 - tanh(rho_star) ^ 2) # if prior of rho is transformed beta
    
    log_prior_gamma_low <- log(prior_param_gamma_low) - log(2) - prior_param_gamma_low * exp(1 / 2 * gamma_low) + 1 / 2 * gamma_low # if prior of sigma_low is exponential
    # log_prior_gamma_low <- dnorm(gamma_low, mean = prior_mu_gamma_low, sd = prior_sigma_gamma_low, log = T) # if prior of sigma_low^2 is normal
    # log_prior_gamma_low <- -log(gamma(prior_shape_gamma_low)) - prior_shape_gamma_low * log(prior_scale_gamma_low) + prior_shape_gamma_low * gamma_low / 2 -
    #   exp(gamma_low / 2) / prior_scale_gamma_low - log(2) # if prior of sigma_low is gamma
    
    log_prior_gamma_high <- log(prior_param_gamma_high) - log(2) - prior_param_gamma_high * exp(1 / 2 * gamma_high) + 1 / 2 * gamma_high # if prior of sigma_high is exponential
    # log_prior_gamma_high <- dnorm(gamma_high, mean = prior_mu_gamma_high, sd = prior_sigma_gamma_high, log = T) # if prior of sigma_high^2 is normal
    # log_prior_gamma_high <- -log(gamma(prior_shape_gamma_high)) - prior_shape_gamma_high * log(prior_scale_gamma_high) + prior_shape_gamma_high * gamma_high / 2 -
    #   exp(gamma_high / 2) / prior_scale_gamma_high - log(2) # if prior of sigma_high is gamma
    
    log_target_value <- log_target_value + log_prior_rho_star + log_prior_gamma_low + log_prior_gamma_high
    
    log_target_value
    
  }
  
  set.seed(l)
  
  # rho_initial <- runif(1, -1, 1)
  rho_initial <- 2 * rbeta(1, shape1 = prior_shape1_rho, shape2 = prior_shape2_rho) - 1
  rho_star_initial <- 1 / 2 * log((1 + rho_initial) / (1 - rho_initial))
  
  sigma_low_initial <- rexp(1, rate = prior_param_gamma_low) # if exponential prior for sigma
  # sigma_low_initial <- rgamma(1, shape = prior_shape_gamma_low, scale = prior_scale_gamma_low) # if gamma prior for sigma
  gamma_low_initial <- log(sigma_low_initial ^ 2)
  # gamma_low_initial <- rnorm(1, mean = prior_mu_gamma_low, sd = prior_sigma_gamma_low) # if log normal prior for sigma^2
  
  sigma_high_initial <- rexp(1, rate = prior_param_gamma_high) # if exponential prior for sigma
  # sigma_high_initial <- rgamma(1, shape = prior_shape_gamma_high, scale = prior_scale_gamma_high) # if gamma prior for sigma
  gamma_high_initial <- log(sigma_high_initial ^ 2)
  # gamma_high_initial <- rnorm(1, mean = prior_mu_gamma_high, sd = prior_sigma_gamma_high) # if log normal prior for sigma^2
  
  while (gamma_high_initial < gamma_low_initial) {
    
    sigma_low_initial <- rexp(1, rate = prior_param_gamma_low) # if exponential prior for sigma
    # sigma_low_initial <- rgamma(1, shape = prior_shape_gamma_low, scale = prior_scale_gamma_low) # if gamma prior for sigma
    gamma_low_initial <- log(sigma_low_initial ^ 2)
    # gamma_low_initial <- rnorm(1, mean = prior_mu_gamma_low, sd = prior_sigma_gamma_low) # if log normal prior for sigma^2
    
    sigma_high_initial <- rexp(1, rate = prior_param_gamma_high) # if exponential prior for sigma
    # sigma_high_initial <- rgamma(1, shape = prior_shape_gamma_high, scale = prior_scale_gamma_high) # if gamma prior for sigma
    gamma_high_initial <- log(sigma_high_initial ^ 2)
    # gamma_high_initial <- rnorm(1, mean = prior_mu_gamma_high, sd = prior_sigma_gamma_high) # if log normal prior for sigma^2
    
  }
  
  # z_initial <- numeric(num_time_points)
  # z_initial[1] <- rnorm(1)
  # for (i in 2:num_time_points) {
  #   z_initial[i] <- rnorm(1, mean = z_initial[i - 1], sd = 1)
  # }
  
  if (l <= 5) {
    z_initial <- rnorm(num_time_points, mean = -1) # an assumption that the majority of the time points are in the low volatility state
  } else {
    z_initial <- rnorm(num_time_points, mean = 1) # an assumption that the majority of the time points are in the low volatility state
  }
  
  qbar_initial <- c(z_initial, rho_star_initial, gamma_low_initial, gamma_high_initial)
  pbar_initial <- rnorm(grhmc_model_list$n_parameters)
  
  region_id <<- grhmc_model_list$target_jump_fun(qbar_initial)
  
  system.time(
    test_run <- grhmc_piecewise_smooth_density_transformed_function(
      model_list = grhmc_model_list,
      lambda = 0.2,
      T = 100000,
      n_samples = 100000,
      diag_s_elements_initial = rep(1, grhmc_model_list$n_parameters),
      m_initial = rep(0, grhmc_model_list$n_parameters),
      qbar_initial = qbar_initial,
      pbar_initial = pbar_initial,
      Lambda_initial = 0,
      u_initial = rexp(1),
      random_state = NULL,
      rtol = NULL,
      atol = NULL,
      h_max = 1.0,
      sampling_compute_temporal_averages_of_moments = FALSE,
      reflection_type = "deterministic", 
      verbose_at_refresh = T,
      last.root.offset.lin.root.finder = 1e-10,
      last.root.offset.non.lin.root.finder = 1e-10,
      precision_real_root_lin_root_finder = 1.0e-13,
      num_subdiv_non_lin_root_finder = 8L
    )
  )
  print("Finished.")
  cbind(test_run$q_original_samples, rep(l, nrow(test_run$q_original_samples)))
  
}
parallel::stopCluster(init_cluster)
# saveRDS(store_matrix_nr2, "piecewise_smooth/example_scripts/section_6/real_data/real_data_switching_volatility_non_adaptive_deterministic_samples_nr2.RDS")
store_matrix_nr2 <- readRDS("piecewise_smooth/example_scripts/section_6/real_data/real_data_switching_volatility_non_adaptive_deterministic_samples_nr2.RDS")

# store_matrix <- rbind(store_matrix_nr1, store_matrix_nr2)
# saveRDS(store_matrix, "piecewise_smooth/example_scripts/section_6/real_data/real_data_switching_volatility_non_adaptive_deterministic_samples.RDS")
store_matrix <- readRDS("piecewise_smooth/example_scripts/section_6/real_data/real_data_switching_volatility_non_adaptive_deterministic_samples.RDS")

store_matrix <- cbind(store_matrix, rep(1:100000, 10))
store_matrix <- store_matrix[store_matrix[, ncol(store_matrix)] > 50000, ] # remove first half as burn in samples
store_matrix

# rho

rho_samples <- tanh(store_matrix[, num_time_points + 1])
plot(rho_samples, col = store_matrix[, num_time_points + 4])
hist(rho_samples, probability = T)
mean(rho_samples)
median(rho_samples)
sd(rho_samples)
quantile(rho_samples, probs = c(0.025, 0.975))

# sigma_low

sigma_low_samples <- exp(0.5 * store_matrix[, num_time_points + 2])
plot(sigma_low_samples, col = store_matrix[, num_time_points + 4])
hist(sigma_low_samples, probability = T)
mean(sigma_low_samples)
median(sigma_low_samples)
sd(sigma_low_samples)
quantile(sigma_low_samples, probs = c(0.025, 0.975))

# sigma_high

sigma_high_samples <- exp(0.5 * store_matrix[, num_time_points + 3])
plot(sigma_high_samples, col = store_matrix[, num_time_points + 4])
hist(sigma_high_samples, probability = T)
mean(sigma_high_samples)
median(sigma_high_samples)
sd(sigma_high_samples)
quantile(sigma_high_samples, probs = c(0.025, 0.975))

# Z_t > 0

df_ggplot <- data.frame(
  y = y_vec,
  posterior_mean_z = colMeans(store_matrix[, 1:num_time_points]),
  posterior_ratio_z_larger_than_zero = sapply(1:num_time_points, function(i) mean(store_matrix[, i] > 0)),
  posterior_quantile_z = t(apply(store_matrix[, 1:num_time_points], 2, function(x) quantile(x, probs = c(0.025, 0.975)))),
  t = 1:num_time_points
)

colnames(df_ggplot)[4:5] <- c("posterior_quantile_z_2_5", "posterior_quantile_z_97_5")

plot_nr1 <- df_ggplot %>%
  ggplot(aes(x = t, y = y_vec)) + 
  geom_line(col = "blue") + 
  labs(y = expression(Y[t]))
plot_nr1

plot_nr2 <- df_ggplot %>%
  ggplot(aes(x = t, y = posterior_mean_z)) +
  geom_line(col = "red") + 
  labs(y = expression(Z[t]))
plot_nr2

plot_nr2 <- df_ggplot %>%
  ggplot() + 
  geom_line(aes(x = t, y = posterior_mean_z), col = "red") + 
  geom_ribbon(aes(x = t, ymin = posterior_quantile_z_2_5, ymax = posterior_quantile_z_97_5), alpha = 0.15, fill = "red", linetype = 0) + 
  labs(y = expression(Z[t]))
plot_nr2

plot_nr3 <- df_ggplot %>% 
  ggplot(aes(x = t, y = posterior_ratio_z_larger_than_zero)) + 
  geom_line(col = "green") + 
  ylim(c(0, 1)) + 
  labs(y = expression("Posterior " * P(Z[t] > 0)))
plot_nr3

# grid.arrange(plot_nr1, plot_nr2, plot_nr3, ncol = 3)
grid.arrange(plot_nr1, plot_nr2, plot_nr3, nrow = 3)
