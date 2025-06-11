rm(list = ls())

library(foreach)
source("disc_grad/implementation_scripts/ISG/full_ISG_grhmc_discontinuous_gradient_transformed_function.R")

set.seed(42)

n_cores <- 10
n_iterations <- 10
n_samples_per_iteration <- 100000
sampling_time_period_per_iteration <- 100000
n_adaptive_samples_per_iteration <- 10000
adaptive_time_period_per_iteration <- 100000

alpha <- 0
delta_1 <- 0.5 
delta_2 <- -0.5
beta_1 <- c(1, 0)
beta_2 <- c(-0.1, 1)
w1 <- w2 <- 1
sigma <- 0.1
gamma <- log(sigma ^ 2)
n <- 100
n_test <- 100
x <- matrix(rnorm(n * 2), ncol = 2)

prior_mu <- rep(0, 9)
prior_sigma <- rep(1, 9)

y <- rnorm(
  n = n,
  mean = alpha + w1 * pmax(0, delta_1 + x %*% beta_1) + w2 * pmax(0, delta_2 + x %*% beta_2),
  sd = sigma
)

x_test <- matrix(rnorm(n * 2), ncol = 2)
y_test <- rnorm(
  n = n_test,
  mean = alpha + w1 * pmax(0, delta_1 + x_test %*% beta_1) + w2 * pmax(0, delta_2 + x_test %*% beta_2),
  sd = sigma
)

og_neuron_1 <- pmax(0, delta_1 + x %*% beta_1)
og_neuron_2 <- pmax(0, delta_2 + x %*% beta_2)

lm(y ~ og_neuron_1 + og_neuron_2)


init_cluster <- parallel::makeCluster(n_cores)
parallel::clusterExport(
  init_cluster, 
  varlist = c(
    "x",
    "y",
    "n",
    "n_samples_per_iteration",
    "sampling_time_period_per_iteration"
  ))
doParallel::registerDoParallel(init_cluster)

samples_matrix <- foreach::foreach(l = 1:10, .combine = "rbind") %dopar% {
  
  # q[1] <- gamma = log(sigma^2)
  # q[2] <- alpha
  # q[3:4] <- w* = log(w) in order to get w >= 0
  # q[5:6] <- delta_j
  # q[7:8] <- beta_1
  # q[9:10] <- beta_2
  
  A <- matrix(nrow = 2 * n, ncol = 10)
  
  for (i in 1:n) { # related to first neuron 
    
    A[i, ] <- c(rep(0, 4), 1, 0, x[i,  ], rep(0, 2)) # essentially q[5] + x %*% q[7:8], i.e. input to the first neuron
    
  }
  
  for (i in 1:n) { # related to second neuron
    
    A[n + i, ] <- c(rep(0, 4), 0, 1, rep(0, 2), x[i, ])
    
  }
  
  A
  
  B <- rep(0, 2 * n)
  
  model_list <- list(
    
    n_parameters = 10,
    
    grad_jump_fun = function(q) {
      
      neuron_indices_list <<- list(
        as.integer(q[5] + x %*% q[7:8] > 0), # observations that will contribute to the first neuron
        as.integer(q[6] + x %*% q[9:10] > 0) # observations that will contribute to the second neuron
      )
      
    },
    
    region_non_lin_root_list = list(root_fun = NULL, event_fun = NULL),
    
    region_lin_root_list = list(A = A, B = B),
    
    sim_q0 = function() rep(1, 10)
    
  )
  
  model_list$log_target_grad <- function(q) {
    # q[1] <- gamma = log(sigma^2)
    # q[2] <- alpha
    # q[3:4] <- w* = log(w) in order to get w >= 0
    # q[5:6] <- delta_j
    # q[7:8] <- beta_1
    # q[9:10] <- beta_2
    
    neuron_1 <- (q[5] + x %*% q[7:8]) * neuron_indices_list[[1]] # check which observation is active for the given neuron
    neuron_2 <- (q[6] + x %*% q[9:10]) * neuron_indices_list[[2]]
    exp_gamma <- exp(q[1])
    w <- exp(q[3:4])
    mu <- q[2] + cbind(neuron_1, neuron_2) %*% w # output of the network
    
    diff_squared_full <- (y - mu) ^ 2
    diff_full <- y - mu
    diff_with_x1_full <- (y - mu) * x[, 1]
    diff_with_x2_full <- (y - mu) * x[, 2]
    diff_with_neuron_1 <- (y - mu) * neuron_1
    diff_with_neuron_2 <- (y - mu) * neuron_2
    
    sum_diff_squared_full <- sum(diff_squared_full)
    sum_diff_full <- sum(diff_full)
    sum_diff_with_x1_full <- sum(diff_with_x1_full)
    sum_diff_with_x2_full <- sum(diff_with_x2_full)
    sum_diff_neuron_1 <- sum(diff_with_neuron_1)
    sum_diff_neuron_2 <- sum(diff_with_neuron_2)
    
    c(
      - n / 2 + (1 / (2 * exp_gamma)) * sum_diff_squared_full - sqrt(exp_gamma) / 2 + 1 / 2, # grad wrt gamma
      - (q[2] - prior_mu[1]) / (prior_sigma[1] ^ 2) + (1 / exp_gamma) * sum_diff_full, # grad wrt alpha
      - (q[3:4] - prior_mu[2:3]) / (prior_sigma[2:3] ^ 2) + c(sum_diff_neuron_1, sum_diff_neuron_2) / exp_gamma * w[1:2], # grad wrt w
      - (q[5:6] - prior_mu[4:5]) / (prior_sigma[4:5] ^ 2)  + c(sum(diff_full * neuron_indices_list[[1]]), sum(diff_full * neuron_indices_list[[2]])) / exp_gamma * w, # grad wrt delta_j
      - (q[7:8] - prior_mu[6:7]) / (prior_sigma[6:7] ^ 2) + c(sum(diff_with_x1_full * neuron_indices_list[[1]]), sum(diff_with_x2_full * neuron_indices_list[[1]])) / exp_gamma * w[1], # grad wrt beta_1
      - (q[9:10] - prior_mu[8:9]) / (prior_sigma[8:9] ^ 2)  + c(sum(diff_with_x1_full * neuron_indices_list[[2]]), sum(diff_with_x2_full * neuron_indices_list[[2]])) / exp_gamma * w[2] # grad wrt beta_2
      
    )
    
  }
  
  set.seed(l)
  qbar_initial <- rnorm(model_list$n_parameters)
  pbar_initial <- rnorm(model_list$n_parameters)
  u_initial <- rexp(1)
  
  neuron_indices_list <<- model_list$grad_jump_fun(qbar_initial)
  
  model_list$log_target_grad(qbar_initial)
  
  system.time(
    neural_network_run <- ISG_grhmc_discontinuous_gradient_transformed_function(
      model_list = model_list,
      lambda_initial = 0.2,
      time_period_adaptive = adaptive_time_period_per_iteration,
      time_period_generating_samples = sampling_time_period_per_iteration,
      n_generated_samples = n_samples_per_iteration,
      n_adaptive_samples = n_adaptive_samples_per_iteration,
      # diag_s_elements_initial = rep(1, model_list$n_parameters),
      diag_s_elements_initial = rep(1, model_list$n_parameters),
      m_initial = rep(0, model_list$n_parameters),
      qbar_initial = qbar_initial,
      pbar_initial = pbar_initial, 
      random_state = NULL, 
      rtol = NULL,
      atol = NULL,
      verbose_at_refresh = TRUE
    )
  )
  
  cbind(neural_network_run$sample_step_run$q_original_samples, rep(l, n_samples_per_iteration))
  
}

parallel::stopCluster(init_cluster)

# saveRDS(samples_matrix, "disc_grad/example_scripts/section_6/bayesian_neural_network/bayesian_neural_network_with_constrained_outer_weights_store_samples.RDS")
samples_matrix <- readRDS("disc_grad/example_scripts/section_6/bayesian_neural_network/bayesian_neural_network_with_constrained_outer_weights_store_samples.RDS")

posterior_samples_sigma <- exp(samples_matrix[, 1] * 0.5)
# plot(posterior_samples_sigma, col = samples_matrix[, 11], ylab = expression(sigma), main = expression(paste("Posterior samples of ", sigma)))
plot(posterior_samples_sigma[seq(from = 50, to = n_samples_per_iteration * n_iterations, by = 50)], col = samples_matrix[seq(from = 50, to = n_samples_per_iteration * n_iterations, by = 50), 11], ylab = expression(sigma), main = expression(paste("Posterior samples of ", sigma)), xlab = "Sample number", pch = 19, cex = 0.1)

mean(posterior_samples_sigma)
sd(posterior_samples_sigma)
quantile(posterior_samples_sigma, probs = c(0.025, 0.975))
