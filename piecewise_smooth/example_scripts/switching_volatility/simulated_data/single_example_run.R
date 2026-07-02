rm(list = ls())
source("piecewise_smooth/implementation_scripts/general_scripts/grhmc_piecewise_smooth_density_transformed_function.R")

set.seed(42)

sigma_y_high <- 5
sigma_y_low <- 2
rho <- -0.75

z1 <- rnorm(1)

relevant_sigma_y <- ifelse(z1 > 0, sigma_y_high, sigma_y_low)
y1 <- rnorm(1, mean = rho * relevant_sigma_y * z1, sd = sqrt((1 - rho ^ 2) * relevant_sigma_y))

num_time_points <- 100
z_vec <- numeric(num_time_points)
z_vec[1] <- z1
y_vec <- numeric(num_time_points)
y_vec[1] <- y1

for (t in 2:num_time_points) {
  
  z_vec[t] <- z_vec[t - 1] + rnorm(1)
  relevant_sigma_y <- ifelse(z_vec[t] > 0, sigma_y_high, sigma_y_low)
  y_vec[t] <- rnorm(1, mean = rho * relevant_sigma_y * (z_vec[t] - z_vec[t - 1]), sd = sqrt((1 - rho ^ 2) * relevant_sigma_y))
  
}

plot(y_vec, type = "l")
lines(z_vec, col = "red")

###########################

prior_shape_gamma_low <- 5 # prior parameters when sigma prior is gamma
prior_scale_gamma_low <- 1

prior_shape_gamma_high <- 5 # prior parameters when sigma prior is gamma
prior_scale_gamma_high <- 2

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
    1 / 2 * (prior_shape_gamma_low - exp(gamma_low / 2) / prior_scale_gamma_low) # prior term if prior of sigma_low is gamma
    
  log_target_grad_value_vec[num_time_points + 3] <- 
    (-1 / 2 + (cosh_rho_star ^ 2 * y_vec[1] ^ 2) / (2 * sigma_high ^ 2) - cosh_rho_star * sinh_rho_star * y_vec[1] * z_t[1] / (2 * sigma_high)) * (region_id[1]) + 
    sum(
      (-1 / 2 + (cosh_rho_star ^ 2 * y_vec[time_indices_2_to_t] ^ 2) / (2 * sigma_high ^ 2) - cosh_rho_star * sinh_rho_star * y_vec[time_indices_2_to_t] * (z_t[time_indices_2_to_t] - z_t[time_indices_1_to_t_minus_1]) / (2 * sigma_high)) * (region_id[time_indices_2_to_t])
    ) + 
    # -prior_param_gamma_high / 2 * exp(0.5 * gamma_high) + 1 / 2 # prior term if prior of sigma_high is exponential
    1 / 2 * (prior_shape_gamma_high - exp(gamma_high / 2) / prior_scale_gamma_high) # prior term if prior of sigma_low is gamma
  
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
  
  log_prior_rho_star <- log(gamma(prior_shape1_rho + prior_shape2_rho)) - log(gamma(prior_shape1_rho)) - log(gamma(prior_shape2_rho)) +
    (prior_shape1_rho - 1) * log((tanh(rho_star) + 1) / 2) + (prior_shape2_rho - 1) * log(1 - (tanh(rho_star) + 1) / 2) - log(2) +
    log(1 - tanh(rho_star) ^ 2) # if prior of rho is transformed beta
  
  # log_prior_gamma_low <- log(prior_param_gamma_low) - log(2) - prior_param_gamma_low * exp(1 / 2 * gamma_low) + 1 / 2 * gamma_low # if prior of sigma_low is exponential
  log_prior_gamma_low <- -log(gamma(prior_shape_gamma_low)) - prior_shape_gamma_low * log(prior_scale_gamma_low) + prior_shape_gamma_low * gamma_low / 2 -
    exp(gamma_low / 2) / prior_scale_gamma_low - log(2) # if prior of sigma_low is gamma
  
  # log_prior_gamma_high <- log(prior_param_gamma_high) - log(2) - prior_param_gamma_high * exp(1 / 2 * gamma_high) + 1 / 2 * gamma_high # if prior of sigma_high is exponential
  log_prior_gamma_high <- -log(gamma(prior_shape_gamma_high)) - prior_shape_gamma_high * log(prior_scale_gamma_high) + prior_shape_gamma_high * gamma_high / 2 -
    exp(gamma_high / 2) / prior_scale_gamma_high - log(2) # if prior of sigma_high is gamma
  
  log_target_value <- log_target_value + log_prior_rho_star + log_prior_gamma_low + log_prior_gamma_high
  
  log_target_value
  
}

set.seed(420)

# rho_initial <- runif(1, -1, 1)
rho_initial <- 2 * rbeta(1, shape1 = prior_shape1_rho, shape2 = prior_shape2_rho) - 1
rho_star_initial <- 1 / 2 * log((1 + rho_initial) / (1 - rho_initial))

# sigma_low_initial <- rexp(1, rate = prior_param_gamma_low) # if exponential prior for sigma
sigma_low_initial <- rgamma(1, shape = prior_shape_gamma_low, scale = prior_scale_gamma_low) # if gamma prior for sigma
gamma_low_initial <- log(sigma_low_initial ^ 2)

# sigma_high_initial <- rexp(1, rate = prior_param_gamma_high) # if exponential prior for sigma
sigma_high_initial <- rgamma(1, shape = prior_shape_gamma_high, scale = prior_scale_gamma_high) # if gamma prior for sigma
gamma_high_initial <- log(sigma_high_initial ^ 2)

while (gamma_high_initial < gamma_low_initial) {
  
  # sigma_low_initial <- rexp(1, rate = prior_param_gamma_low) # if exponential prior for sigma
  sigma_low_initial <- rgamma(1, shape = prior_shape_gamma_low, scale = prior_scale_gamma_low) # if gamma prior for sigma
  gamma_low_initial <- log(sigma_low_initial ^ 2)

  # sigma_high_initial <- rexp(1, rate = prior_param_gamma_high) # if exponential prior for sigma
  sigma_high_initial <- rgamma(1, shape = prior_shape_gamma_high, scale = prior_scale_gamma_high) # if gamma prior for sigma
  gamma_high_initial <- log(sigma_high_initial ^ 2)

}

# z_initial <- numeric(num_time_points)
# z_initial[1] <- rnorm(1)
# for (i in 2:num_time_points) {
#   z_initial[i] <- rnorm(1, mean = z_initial[i - 1], sd = 1)
# }

z_initial <- rnorm(num_time_points, mean = -1)

qbar_initial <- c(z_initial, rho_star_initial, gamma_low_initial, gamma_high_initial)
pbar_initial <- rnorm(grhmc_model_list$n_parameters)

region_id <<- grhmc_model_list$target_jump_fun(qbar_initial)

system.time(
  test_run <- grhmc_piecewise_smooth_density_transformed_function(
    model_list = grhmc_model_list,
    lambda = 0.2,
    T = 10000,
    n_samples = 10000,
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

# rho_star and rho

plot(test_run$q_original_samples[, num_time_points + 1])
rho_samples <- tanh(test_run$q_original_samples[, num_time_points + 1])
plot(rho_samples)
hist(rho_samples)
mean(rho_samples)
median(rho_samples)
mean(rho_samples[-c(1:5000)])
median(rho_samples[-c(1:5000)])

# gamma_low and sigma_low

plot(test_run$q_original_samples[, num_time_points + 2])
sigma_low_samples <- exp(0.5 * test_run$q_original_samples[, num_time_points + 2])
plot(sigma_low_samples)
hist(sigma_low_samples)
mean(sigma_low_samples)
median(sigma_low_samples)
mean(sigma_low_samples[-c(1:5000)])
median(sigma_low_samples[-c(1:5000)])

# gamma_low and sigma_low

plot(test_run$q_original_samples[, num_time_points + 3])
sigma_high_samples <- exp(0.5 * test_run$q_original_samples[, num_time_points + 3])
plot(sigma_high_samples)
hist(sigma_high_samples)
mean(sigma_high_samples)
median(sigma_high_samples)
mean(sigma_high_samples[-c(1:5000)])
median(sigma_high_samples[-c(1:5000)])
