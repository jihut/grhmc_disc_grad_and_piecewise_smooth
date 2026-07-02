rm(list = ls())
library(doParallel)
library(dplyr)
library(tidyr)
library(gridExtra)

# Reflective HMC implementation tailored for this setting

single_iteration_reflective_hmc <- function(
    model_list, 
    h, 
    L, 
    q_initial, 
    p_initial,
    num_time_points,
    grad_q_initial
) {
  hamiltonian_initial <- -model_list$log_target_fun(q_initial) + sum(p_initial ^ 2) / 2
  p <- p_initial
  q <- q_initial
  grad_current_q <- grad_q_initial
  current_sign_vec <- sign(q) # used to ensure change in region when checking potential energy difference
  for (i in 1:L) {
    p <- p + h * grad_current_q / 2
    t_0 <- 0
    length_time_point_collision <- num_time_points
    while(length_time_point_collision != 0) {
      # print("t_0")
      # print(t_0)
      # store_q <<- q
      # store_p <<- p
      current_q <- q
      q <- q + (h - t_0) * p
      check_collision <- current_q[1:num_time_points] * q[1:num_time_points]
      which_time_point_collision <- which(check_collision < 0)
      # print("which_time_point_collision")
      # print(which_time_point_collision)
      length_time_point_collision <- length(which_time_point_collision)
      if (length_time_point_collision != 0) {
        t_hit_vec <- -current_q[which_time_point_collision] / p[which_time_point_collision]
        t_hit <- min(t_hit_vec)
        # print("t_hit")
        # print(t_hit)
        index_hit <- which_time_point_collision[which.min(t_hit_vec)]
        # print("index_hit")
        # print(index_hit)
        t_0 <- t_0 + t_hit
        q <- current_q + t_hit * p
        q_hit <<- q
        q[index_hit] <- 0 # ensure that one does not get "double" collision as the boundary in each dimension is simply at zero
        normal_vec <- numeric(model_list$n_parameters)
        normal_vec[index_hit] <- 1
        p_perpendicular <- sum(p * normal_vec) / sum(normal_vec ^ 2) * normal_vec
        p_parallel <- p - p_perpendicular
        norm_squared_p_perpendicular <- sum(p_perpendicular ^ 2)
        old_potential_energy <- -model_list$log_target_fun(q)
        new_sign_vec <- current_sign_vec
        new_sign_vec[index_hit] <- -1 * new_sign_vec[index_hit]
        current_region_id <- region_id
        # print(current_sign_vec)
        # print(new_sign_vec)
        model_list$target_jump_fun(new_sign_vec)
        new_region_id <- region_id
        # print(current_region_id)
        # print(new_region_id)
        new_potential_energy <- -model_list$log_target_fun(q)
        delta_U <- new_potential_energy - old_potential_energy 
        if (delta_U == 0) {
          if (mean(new_region_id == current_region_id) == 1) {
            stop("delta_U == 0, not good!")
          }
        }
        # print("delta_U")
        # print(delta_U)
        # print("norm_squared_p_perpendicular")
        # print(norm_squared_p_perpendicular)
        
        if (norm_squared_p_perpendicular > 2 * delta_U) {
          # print("refract")
          p_perpendicular <- sqrt(norm_squared_p_perpendicular - 2 * delta_U) * p_perpendicular / 
            sqrt(norm_squared_p_perpendicular) # refract the momentum along the normal direction after moving to a different region
          
          p <- p_parallel + p_perpendicular
          
          current_sign_vec <- new_sign_vec
          
        } else {
          # print("reflection")
          p_perpendicular <- -p_perpendicular
          
          p <- p_parallel + p_perpendicular
          
          model_list$target_jump_fun(current_sign_vec) # need to update back to the old region
          
        }
        
      } else { # no boundary collision event --> full half step update of momentum at the end
        grad_current_q <- model_list$log_target_grad(q)
        p <- p + h * grad_current_q / 2
      }
      
    }
    
  }
  
  hamiltonian_final <- -model_list$log_target_fun(q) + sum(p ^ 2) / 2
  
  if (runif(1) < exp(-(hamiltonian_final - hamiltonian_initial))) {
    final_q <- q
    final_p <- -p
    accept_indicator <- 1
    final_grad_q <- grad_current_q
  } else {
    final_q <- q_initial
    final_p <- p_initial
    accept_indicator <- 0
    final_grad_q <- grad_q_initial
    model_list$target_jump_fun(q_initial)
  }
  list(
    final_q = final_q, 
    final_p = final_p, 
    accept_indicator = accept_indicator, 
    hamiltonian_initial = hamiltonian_initial,
    hamiltonian_final = hamiltonian_final,
    final_grad_q = final_grad_q
  )
  
}

#####

data <- read.table("piecewise_smooth/example_scripts/switching_volatility/real_data/data.dat")
y_vec <- data$V1

num_time_points <- length(y_vec)

init_cluster <- parallel::makeCluster(10)
doParallel::registerDoParallel(init_cluster)

final_run <- foreach::foreach(l = 1:10) %dopar% {
  sink(paste0("piecewise_smooth/example_scripts/switching_volatility/real_data/log/log_nr", l, ".txt"))
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
      
    }
    
  )
  n_evals_ode <- 0
  grhmc_model_list$log_target_grad <- function(q) {
    n_evals_ode <<- 1 + n_evals_ode
    time_indices_2_to_t_minus_1 <- 2:(num_time_points - 1)
    time_indices_3_to_t <- 3:(num_time_points)
    time_indices_1_to_t_minus_2 <- 1:(num_time_points - 2)
    time_indices_2_to_t <- c(time_indices_2_to_t_minus_1, num_time_points)
    time_indices_1_to_t_minus_1 <- c(time_indices_1_to_t_minus_2, num_time_points - 1)
    
    
    rho_star <- q[num_time_points + 1]
    rho <- tanh(rho_star)
    cosh_rho_star <- cosh(rho_star)
    sinh_rho_star <- sinh(rho_star)
    
    gamma_star_low <- q[num_time_points + 2]
    exp_gamma_star_low <- exp(gamma_star_low) 
    gamma_low <- gamma_star_low
    sigma_low <- exp(gamma_low / 2)
    
    gamma_star_high <- q[num_time_points + 3]
    exp_gamma_star_high <- exp(gamma_star_high) 
    gamma_high <- log(exp(gamma_star_low) + exp(gamma_star_high))
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
    
    # log_target_grad_value_vec[num_time_points + 2] <- 
    #   (-1 / 2 + (cosh_rho_star ^ 2 * y_vec[1] ^ 2) / (2 * sigma_low ^ 2) - cosh_rho_star * sinh_rho_star * y_vec[1] * z_t[1] / (2 * sigma_low)) * (1 - region_id[1]) + 
    #   sum(
    #     (-1 / 2 + (cosh_rho_star ^ 2 * y_vec[time_indices_2_to_t] ^ 2) / (2 * sigma_low ^ 2) - cosh_rho_star * sinh_rho_star * y_vec[time_indices_2_to_t] * (z_t[time_indices_2_to_t] - z_t[time_indices_1_to_t_minus_1]) / (2 * sigma_low)) * (1 - region_id[time_indices_2_to_t])
    #   ) + 
    #   -prior_param_gamma_low / 2 * exp(0.5 * gamma_low) + 1 / 2 # prior term if prior of sigma_low is exponential
    
    # log_target_grad_value_vec[num_time_points + 3] <- 
    #   (-1 / 2 + (cosh_rho_star ^ 2 * y_vec[1] ^ 2) / (2 * sigma_high ^ 2) - cosh_rho_star * sinh_rho_star * y_vec[1] * z_t[1] / (2 * sigma_high)) * (region_id[1]) + 
    #   sum(
    #     (-1 / 2 + (cosh_rho_star ^ 2 * y_vec[time_indices_2_to_t] ^ 2) / (2 * sigma_high ^ 2) - cosh_rho_star * sinh_rho_star * y_vec[time_indices_2_to_t] * (z_t[time_indices_2_to_t] - z_t[time_indices_1_to_t_minus_1]) / (2 * sigma_high)) * (region_id[time_indices_2_to_t])
    #   ) + 
    #   -prior_param_gamma_high / 2 * exp(0.5 * gamma_high) + 1 / 2 # prior term if prior of sigma_high is exponential
    
    log_target_grad_value_gamma_low <-
      (-1 / 2 + (cosh_rho_star ^ 2 * y_vec[1] ^ 2) / (2 * sigma_low ^ 2) - cosh_rho_star * sinh_rho_star * y_vec[1] * z_t[1] / (2 * sigma_low)) * (1 - region_id[1]) +
      sum(
        (-1 / 2 + (cosh_rho_star ^ 2 * y_vec[time_indices_2_to_t] ^ 2) / (2 * sigma_low ^ 2) - cosh_rho_star * sinh_rho_star * y_vec[time_indices_2_to_t] * (z_t[time_indices_2_to_t] - z_t[time_indices_1_to_t_minus_1]) / (2 * sigma_low)) * (1 - region_id[time_indices_2_to_t])
      ) +
      -prior_param_gamma_low / 2 * exp(0.5 * gamma_low) + 1 / 2 # prior term if prior of sigma_low is exponential
    
    log_target_grad_value_gamma_high <-
      (-1 / 2 + (cosh_rho_star ^ 2 * y_vec[1] ^ 2) / (2 * sigma_high ^ 2) - cosh_rho_star * sinh_rho_star * y_vec[1] * z_t[1] / (2 * sigma_high)) * (region_id[1]) +
      sum(
        (-1 / 2 + (cosh_rho_star ^ 2 * y_vec[time_indices_2_to_t] ^ 2) / (2 * sigma_high ^ 2) - cosh_rho_star * sinh_rho_star * y_vec[time_indices_2_to_t] * (z_t[time_indices_2_to_t] - z_t[time_indices_1_to_t_minus_1]) / (2 * sigma_high)) * (region_id[time_indices_2_to_t])
      ) +
      -prior_param_gamma_high / 2 * exp(0.5 * gamma_high) + 1 / 2 # prior term if prior of sigma_high is exponential
    
    log_target_grad_value_vec[num_time_points + 2] <- 
      log_target_grad_value_gamma_low * 1 + 
      log_target_grad_value_gamma_high * exp_gamma_star_low / (exp_gamma_star_low + exp_gamma_star_high) # chain rule, using grad of gamma_high and gamma_low to get grad of gamma_star_high etc. via Jacobian
    
    log_target_grad_value_vec[num_time_points + 3] <- 
      log_target_grad_value_gamma_high * exp_gamma_star_high / (exp_gamma_star_low + exp_gamma_star_high)
    
    log_target_grad_value_vec
    
  }
  
  grhmc_model_list$log_target_fun <- function(q) {
    time_indices_2_to_t <- 2:num_time_points
    time_indices_1_to_t_minus_1 <- time_indices_2_to_t - 1
    rho_star <- q[num_time_points + 1]
    rho <- tanh(rho_star)
    cosh_rho_star <- cosh(rho_star)
    sinh_rho_star <- sinh(rho_star)
    
    # gamma_low <- q[num_time_points + 2]
    # sigma_low <- exp(gamma_low / 2)
    # 
    # gamma_high <- q[num_time_points + 3]
    # sigma_high <- exp(gamma_high / 2)
    
    gamma_star_low <- q[num_time_points + 2]
    exp_gamma_star_low <- exp(gamma_star_low) 
    gamma_low <- gamma_star_low
    sigma_low <- exp(gamma_low / 2)
    
    gamma_star_high <- q[num_time_points + 3]
    exp_gamma_star_high <- exp(gamma_star_high) 
    gamma_high <- log(exp(gamma_star_low) + exp(gamma_star_high))
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
    
    log_prior_gamma_star_jacobian <- log(
      exp_gamma_star_high / (exp_gamma_star_high + exp_gamma_star_low) # jacobian determinant of transformation (gamma_star_low, gamma_star_high) --> (gamma_low, gamma_high)
    )
    
    log_target_value <- log_target_value + log_prior_rho_star + log_prior_gamma_low + log_prior_gamma_high + log_prior_gamma_star_jacobian
    
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
  gamma_star_low_initial <- gamma_low_initial
  
  sigma_high_initial <- rexp(1, rate = prior_param_gamma_high) # if exponential prior for sigma
  # sigma_high_initial <- rgamma(1, shape = prior_shape_gamma_high, scale = prior_scale_gamma_high) # if gamma prior for sigma
  gamma_high_initial <- log(sigma_high_initial ^ 2)
  # gamma_high_initial <- rnorm(1, mean = prior_mu_gamma_high, sd = prior_sigma_gamma_high) # if log normal prior for sigma^2
  gamma_star_high_initial <- log(exp(gamma_high_initial) - exp(gamma_star_low_initial))
  
  while (gamma_high_initial < gamma_low_initial) {
    
    sigma_low_initial <- rexp(1, rate = prior_param_gamma_low) # if exponential prior for sigma
    # sigma_low_initial <- rgamma(1, shape = prior_shape_gamma_low, scale = prior_scale_gamma_low) # if gamma prior for sigma
    gamma_low_initial <- log(sigma_low_initial ^ 2)
    # gamma_low_initial <- rnorm(1, mean = prior_mu_gamma_low, sd = prior_sigma_gamma_low) # if log normal prior for sigma^2
    gamma_star_low_initial <- gamma_low_initial
    
    sigma_high_initial <- rexp(1, rate = prior_param_gamma_high) # if exponential prior for sigma
    # sigma_high_initial <- rgamma(1, shape = prior_shape_gamma_high, scale = prior_scale_gamma_high) # if gamma prior for sigma
    gamma_high_initial <- log(sigma_high_initial ^ 2)
    # gamma_high_initial <- rnorm(1, mean = prior_mu_gamma_high, sd = prior_sigma_gamma_high) # if log normal prior for sigma^2
    gamma_star_high_initial <- log(exp(gamma_high_initial) - exp(gamma_star_low_initial))
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
  
  q_initial <- c(z_initial, rho_star_initial, gamma_star_low_initial, gamma_star_high_initial)
  p_initial <- rnorm(grhmc_model_list$n_parameters)
  
  region_id <<- grhmc_model_list$target_jump_fun(q_initial)
  
  n_warmup_nr1_iterations <- 5000
  # store_warmup_nr1_samples <- matrix(nrow = n_warmup_nr1_iterations, ncol = grhmc_model_list$n_parameters)
  # store_warmup_nr1_accept_indicator <- numeric(n_warmup_nr1_iterations)
  # store_warmup_nr1_energy_error <- numeric(n_warmup_iterations)
  grhmc_model_list$target_jump_fun(q_initial)
  grad_q_initial <- grhmc_model_list$log_target_grad(q_initial)
  
  path_length <- 0.25
  h_warmup <- 0.005
  L_warmup <- ceiling(path_length / h_warmup)
  print("warmup nr1")
  for (i in 1:n_warmup_nr1_iterations) {
    print(i)
    single_iteration <- single_iteration_reflective_hmc(
      model_list = grhmc_model_list,
      h = h_warmup, 
      L = L_warmup,
      q_initial = q_initial,
      p_initial = rnorm(grhmc_model_list$n_parameters),
      num_time_points = length(y_vec),
      grad_q_initial = grad_q_initial
    )
    q_initial <- single_iteration$final_q
    grad_q_initial <- single_iteration$final_grad_q
    # store_warmup_nr1_samples[i, ] <- q_initial
    # store_warmup_nr1_accept_indicator[i] <- single_iteration$accept_indicator
    # store_warmup_nr1_energy_error[i] <- single_iteration$hamiltonian_final - single_iteration$hamiltonian_initial 
  }
  
  n_warmup_nr2_iterations <- 5000
  # store_warmup_nr2_samples <- matrix(nrow = n_warmup_nr2_iterations, ncol = grhmc_model_list$n_parameters)
  # store_warmup_nr2_accept_indicator <- numeric(n_warmup_nr2_iterations)
  # store_warmup_nr2_energy_error <- numeric(n_warmup_nr2_iterations)
  
  path_length <- 1
  h_warmup <- 0.01
  L_warmup <- ceiling(path_length / h_warmup)
  print("warmup nr2")
  for (i in 1:n_warmup_nr2_iterations) {
    print(i)
    single_iteration <- single_iteration_reflective_hmc(
      model_list = grhmc_model_list,
      h = h_warmup, 
      L = L_warmup,
      q_initial = q_initial,
      p_initial = rnorm(grhmc_model_list$n_parameters),
      num_time_points = length(y_vec),
      grad_q_initial = grad_q_initial
    )
    q_initial <- single_iteration$final_q
    grad_q_initial <- single_iteration$final_grad_q
    # store_warmup_nr2_samples[i, ] <- q_initial
    # store_warmup_nr2_accept_indicator[i] <- single_iteration$accept_indicator
    # store_warmup_nr2_energy_error[i] <- single_iteration$hamiltonian_final - single_iteration$hamiltonian_initial
  }
  
  n_warmup_nr3_iterations <- 990000
  # store_warmup_nr2_samples <- matrix(nrow = n_warmup_nr2_iterations, ncol = grhmc_model_list$n_parameters)
  # store_warmup_nr2_accept_indicator <- numeric(n_warmup_nr2_iterations)
  # store_warmup_nr2_energy_error <- numeric(n_warmup_nr2_iterations)
  
  path_length <- 5
  h_warmup <- 0.025
  L_warmup <- ceiling(path_length / h_warmup)
  print("warmup nr3")
  for (i in 1:n_warmup_nr3_iterations) {
    print(i)
    single_iteration <- single_iteration_reflective_hmc(
      model_list = grhmc_model_list,
      h = h_warmup, 
      L = L_warmup,
      q_initial = q_initial,
      p_initial = rnorm(grhmc_model_list$n_parameters),
      num_time_points = length(y_vec),
      grad_q_initial = grad_q_initial
    )
    q_initial <- single_iteration$final_q
    grad_q_initial <- single_iteration$final_grad_q
    # store_warmup_nr2_samples[i, ] <- q_initial
    # store_warmup_nr2_accept_indicator[i] <- single_iteration$accept_indicator
    # store_warmup_nr2_energy_error[i] <- single_iteration$hamiltonian_final - single_iteration$hamiltonian_initial
  }
  
  n_evals_ode <- 0 # reset gradient eval counter
  n_sampling_iterations <- 1000000
  store_sampling_samples <- matrix(nrow = n_sampling_iterations, ncol = grhmc_model_list$n_parameters)
  store_sampling_accept_indicator <- numeric(n_sampling_iterations)
  store_sampling_energy_error <- numeric(n_sampling_iterations)
  
  path_length <- 5
  h_sampling <- 0.03
  L_sampling <- ceiling(path_length / h_sampling)
  print("sampling")
  for (i in 1:n_sampling_iterations) {
    print(i)
    single_iteration <- single_iteration_reflective_hmc(
      model_list = grhmc_model_list,
      h = h_sampling, 
      L = L_sampling,
      q_initial = q_initial,
      p_initial = rnorm(grhmc_model_list$n_parameters),
      num_time_points = length(y_vec),
      grad_q_initial = grad_q_initial
    )
    q_initial <- single_iteration$final_q
    grad_q_initial <- single_iteration$final_grad_q
    store_sampling_samples[i, ] <- q_initial
    store_sampling_accept_indicator[i] <- single_iteration$accept_indicator
    store_sampling_energy_error[i] <- single_iteration$hamiltonian_final - single_iteration$hamiltonian_initial 
  }
  
  list(
    samples = store_sampling_samples,
    accept_indicator = store_sampling_accept_indicator,
    energy_error = store_sampling_energy_error,
    n_evals_ode = n_evals_ode
  )
  
}

saveRDS(final_run, "piecewise_smooth/example_scripts/switching_volatility/real_data/real_data_switching_volatility_with_constraint_reflective_hmc.RDS")

parallel::stopCluster(init_cluster)

final_result <- # concatenate same-named elements together across the chains
  lapply(names(final_run[[1]]), function(element_name) {
    elements <- lapply(final_run, '[[', element_name)
    if (is.matrix(elements[[1]])) {
      do.call(rbind, elements)
    } else {
      unlist(elements)
    }
  })

names(final_result) <- names(final_run[[1]])
final_result

final_samples_3d_array <- abind::abind(lapply(final_run, '[[', "samples"), along = 3)
final_samples_3d_array <- aperm(final_samples_3d_array, perm = c(1, 3, 2))
dim(final_samples_3d_array)

rstan_monitor_summary <- rstan::monitor(final_samples_3d_array, warmup = 0)
rstan_monitor_summary$n_eff
rstan_monitor_summary$n_eff * 1000000 / sum(final_result$n_evals_ode)

store_matrix <- cbind(final_result$samples, chain_index = sort(rep(1:10, 20000)))

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

# sigma_high_samples <- exp(0.5 * store_matrix[, num_time_points + 3])
sigma_high_samples <- exp(0.5 * log(exp(store_matrix[, num_time_points + 2]) + exp(store_matrix[, num_time_points + 3])))
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
