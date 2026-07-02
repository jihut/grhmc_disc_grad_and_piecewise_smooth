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

init_cluster <- parallel::makeCluster(10)
doParallel::registerDoParallel(init_cluster)

integrator_tol <- c(1e-3, 1e-4, 1e-5, 1e-6)
n_trajectories <- 100

list_chains_tol_scenario <- list()

for (k in 1:length(integrator_tol)) {
  doRNG::registerDoRNG(seed = 42)
  
  store_chains <- foreach::foreach(i = 1:n_trajectories) %dopar% {
    
    max_model <- list(
      
      n_parameters = 2,
      
      grad_jump_fun = function(q) {
        
        grad_indices <<- c(
          as.integer(q[1] > 0)
        )
        
      },
      
      additional_non_lin_root_list = NULL, # usually: list(root_fun = ..., event_fun = ...)
      
      additional_lin_root_list = NULL, # usually: list(A = ..., B = ..., event_fun = ...)
      
      region_lin_root_list = list(A = matrix(c(1, 0), nrow = 1), B = 0), # also need to specify the break point here again
      
      sim_q0 = function() rnorm(2)
      
    )
    
    max_model$log_target_grad <- function(q) { # define the log target grad
      
      c(
        -q[1] + (q[2] - q[1]) * grad_indices,
        -q[2] + q[1] * grad_indices
      )
      
    }
    
    # model_list <- max_model
    # 
    # qbar_initial <- model_list$sim_q0()
    # pbar_initial <- rnorm(2)
    # 
    # test_run <- ISG_grhmc_discontinuous_gradient_transformed_function(
    #   model_list = model_list,
    #   lambda_initial = 0.2,
    #   time_period_adaptive = 10000,
    #   time_period_generating_samples = 10000,
    #   n_generated_samples = 10000,
    #   n_adaptive_samples = 10000,
    #   diag_s_elements_initial = rep(1, 2),
    #   m_initial = rep(0, 2),
    #   qbar_initial = qbar_initial,
    #   pbar_initial = pbar_initial, 
    #   random_state = NULL, 
    #   rtol = integrator_tol[k],
    #   atol = integrator_tol[k],
    #   verbose_at_refresh = T
    # )
    # 
    # test_run$sample_step_run
    
    model_list <- max_model
    qbar_initial <- model_list$sim_q0()
    pbar_initial <- rnorm(2)
    u_initial <- rexp(1)
    
    test_run_warmup <- grhmc_discontinuous_gradient_transformed_function(
      model_list = model_list,
      lambda = 0.2,
      T = 10000,
      n_samples = 10000,
      diag_s_elements_initial = rep(1, 2),
      m_initial = rep(0, 2),
      qbar_initial = qbar_initial,
      pbar_initial = pbar_initial,
      Lambda_initial = 0,
      u_initial = u_initial,
      random_state = NULL,
      rtol = integrator_tol[k],
      atol = integrator_tol[k],
      verbose_at_refresh = T
    )
    
    new_pbar_initial <- rnorm(2)
    new_u_initial <- rexp(1)
    
    test_run_sampling <- grhmc_discontinuous_gradient_transformed_function(
      model_list = model_list,
      lambda = 0.2,
      T = 10000,
      n_samples = 100000,
      diag_s_elements_initial = rep(1, 2),
      m_initial = rep(0, 2),
      qbar_initial = test_run_warmup$q_original_samples[nrow(test_run_warmup$q_original_samples), ],
      pbar_initial = new_pbar_initial,
      Lambda_initial = 0,
      u_initial = new_u_initial,
      random_state = NULL,
      rtol = integrator_tol[k],
      atol = integrator_tol[k],
      verbose_at_refresh = T
    )
    
    test_run_sampling
    
  }
  
  list_chains_tol_scenario[[k]] <- store_chains
  
}

set.seed(42)
q1 <- rnorm(1e8)
q2_given_q1 <- rnorm(1e8, mean = pmax(0, q1), sd = 1)
mean(q2_given_q1)

true_mean_q2 <- 
  integrate(function(x) x * true_pdf_q2(x), lower = -Inf, upper = Inf, rel.tol = 1e-13, abs.tol = 1e-13)
true_mean_q2
true_mean_q2_squared <- 
  integrate(function(x) (x ^ 2) * true_pdf_q2(x), lower = -Inf, upper = Inf, rel.tol = 1e-13, abs.tol = 1e-13)
true_mean_q2_squared

integrand_true_mean_q_dot_q <- function(q1, q2) {
  (q1 ^ 2 + q2 ^ 2) * dnorm(q1, mean = 0, sd = 1) * dnorm(q2, mean = max(0, q1), sd = 1)
}

true_mean_q_dot_q <- integrate(function(q1) {
  sapply(q1, function(q1_star) {
    integrate(function(q2) integrand_true_mean_q_dot_q(q1_star, q2), lower = -Inf, upper = Inf, rel.tol = 1e-13, abs.tol = 1e-13)$value
  })
}, lower = -Inf, upper = Inf, rel.tol = 1e-13, abs.tol = 1e-13)

true_mean_q_dot_q

estimated_quantities <- function(list, calculate_wasserstein_2 = FALSE) {
  
  samples_matrix <- list$q_original_samples
  
  mean_q2 <- mean(samples_matrix[, 2])
  mean_q2_squared <- mean(samples_matrix[, 2] ^ 2)
  mean_q_dot_q <- mean(rowSums(samples_matrix * samples_matrix))
  
  if (calculate_wasserstein_2 == TRUE) {
    
    wasserstein_2_integrand <- function(t) {
      
      # quantile_via_cdf <- uniroot(
      #   function(x) true_cdf_q2(x) - t, lower = -10, upper = 10, tol = 1e-13
      # )$root
      #
      # quantile_via_sample <- quantile(list_chains_tol_scenario[[1]][[1]]$q_original_samples[, 2], probs = t)
      #
      # (quantile_via_cdf - quantile_via_sample) ^ 2
      
      quantile_via_cdf <- sapply(
        t,
        function(t_star) {
          uniroot(
            function(x) true_cdf_q2(x) - t_star, lower = -10, upper = 10, tol = 1e-12
          )$root
        }
      )
      
      quantile_via_sample <- quantile(samples_matrix[, 2], probs = t)
      
      (quantile_via_cdf - quantile_via_sample) ^ 2
      
    }
    
    wasserstein_2_squared_value <- 0
    grid <- seq(from = 0, to = 1, by = 0.0001)
    
    for (i in 1:(length(grid) - 1)) {
      # print(i)
      int_run <- integrate(
        wasserstein_2_integrand,
        lower = grid[i],
        upper = grid[i + 1],
        rel.tol = 1e-9,
        abs.tol = 1e-9,
        subdivisions = 100000
      )
      wasserstein_2_squared_value <- wasserstein_2_squared_value + int_run$value
    } 
    
    return(
      c(
        q2 = mean_q2, q2_squared = mean_q2_squared, q_dot_q = mean_q_dot_q, wasserstein_2 = sqrt(wasserstein_2_squared_value)
      )
    )
    
  } else {
    return(
      c(
        q2 = mean_q2, q2_squared = mean_q2_squared, q_dot_q = mean_q_dot_q
      )
    )  
  }
  
}

for (k in 1:length(integrator_tol)) {
  relevant_list <- list_chains_tol_scenario[[k]]
  relevant_matrix <- cbind(do.call(rbind, lapply(relevant_list, estimated_quantities)), 
                           int_tol = rep(integrator_tol[k]))
  if (k == 1) {
    final_result_matrix <- relevant_matrix
  } else {
    final_result_matrix <- rbind(final_result_matrix, relevant_matrix)
  }
}

final_result_matrix <- as.data.frame(final_result_matrix)

diff_result_matrix <- 
  sweep(final_result_matrix, 2, c(true_mean_q2$value, true_mean_q2_squared$value, true_mean_q_dot_q$value, 0), "-")

# Plot of difference between estimates of E(q_2) and true E(q_2) 
par(mfrow = c(1, 3))

plot(
  exp(jitter(log(diff_result_matrix$int_tol), amount = 0.25)), 
  diff_result_matrix$q2, 
  log = "x",
  xaxt = "n",
  xlab = "Integrator tolerance",
  ylab = "Deviation from true value",
  main = expression(E(q[2]))
)

axis(1, at = integrator_tol, labels = integrator_tol)
points(
  integrator_tol, 
  sapply(1:4, function(i) mean(diff_result_matrix$q2[diff_result_matrix$int_tol == integrator_tol[i]])),
  col = "red",
  pch = 15
)
segments(
  integrator_tol, 
  sapply(1:4, function(i) mean(diff_result_matrix$q2[diff_result_matrix$int_tol == integrator_tol[i]]) - 
    sd(diff_result_matrix$q2[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  integrator_tol, 
  sapply(1:4, function(i) mean(diff_result_matrix$q2[diff_result_matrix$int_tol == integrator_tol[i]]) + 
           sd(diff_result_matrix$q2[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  col = "red"
)
cap <- 0.25
segments(
  integrator_tol / (1 + cap), 
  sapply(1:4, function(i) mean(diff_result_matrix$q2[diff_result_matrix$int_tol == integrator_tol[i]]) - 
           sd(diff_result_matrix$q2[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  integrator_tol * (1 + cap), 
  sapply(1:4, function(i) mean(diff_result_matrix$q2[diff_result_matrix$int_tol == integrator_tol[i]]) - 
           sd(diff_result_matrix$q2[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  col = "red"
)
segments(
  integrator_tol / (1 + cap), 
  sapply(1:4, function(i) mean(diff_result_matrix$q2[diff_result_matrix$int_tol == integrator_tol[i]]) + 
           sd(diff_result_matrix$q2[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  integrator_tol * (1 + cap), 
  sapply(1:4, function(i) mean(diff_result_matrix$q2[diff_result_matrix$int_tol == integrator_tol[i]]) + 
           sd(diff_result_matrix$q2[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  col = "red"
)
abline(h = 0, lty = 2, col = "blue")

# Plot of difference between estimates of E(q_2 ^ 2) and true E(q_2 ^ 2) 

plot(
  exp(jitter(log(diff_result_matrix$int_tol), amount = 0.25)), 
  diff_result_matrix$q2_squared, 
  log = "x",
  xaxt = "n",
  xlab = "Integrator tolerance",
  ylab = "Deviation from true value",
  main = expression(E(q[2]^2))
)

axis(1, at = integrator_tol, labels = integrator_tol)
points(
  integrator_tol, 
  sapply(1:4, function(i) mean(diff_result_matrix$q2_squared[diff_result_matrix$int_tol == integrator_tol[i]])),
  col = "red",
  pch = 15
)
segments(
  integrator_tol, 
  sapply(1:4, function(i) mean(diff_result_matrix$q2_squared[diff_result_matrix$int_tol == integrator_tol[i]]) - 
           sd(diff_result_matrix$q2_squared[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  integrator_tol, 
  sapply(1:4, function(i) mean(diff_result_matrix$q2_squared[diff_result_matrix$int_tol == integrator_tol[i]]) + 
           sd(diff_result_matrix$q2_squared[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  col = "red"
)
cap <- 0.25
segments(
  integrator_tol / (1 + cap), 
  sapply(1:4, function(i) mean(diff_result_matrix$q2_squared[diff_result_matrix$int_tol == integrator_tol[i]]) - 
           sd(diff_result_matrix$q2_squared[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  integrator_tol * (1 + cap), 
  sapply(1:4, function(i) mean(diff_result_matrix$q2_squared[diff_result_matrix$int_tol == integrator_tol[i]]) - 
           sd(diff_result_matrix$q2_squared[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  col = "red"
)
segments(
  integrator_tol / (1 + cap), 
  sapply(1:4, function(i) mean(diff_result_matrix$q2_squared[diff_result_matrix$int_tol == integrator_tol[i]]) + 
           sd(diff_result_matrix$q2_squared[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  integrator_tol * (1 + cap), 
  sapply(1:4, function(i) mean(diff_result_matrix$q2_squared[diff_result_matrix$int_tol == integrator_tol[i]]) + 
           sd(diff_result_matrix$q2_squared[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  col = "red"
)
abline(h = 0, lty = 2, col = "blue")

# Plot of difference between estimates of E(q dot q) and true E(q dot q) 

plot(
  exp(jitter(log(diff_result_matrix$int_tol), amount = 0.25)), 
  diff_result_matrix$q_dot_q, 
  log = "x",
  xaxt = "n",
  xlab = "Integrator tolerance",
  ylab = "Deviation from true value",
  main = expression(E(q^T * q))
)

axis(1, at = integrator_tol, labels = integrator_tol)
points(
  integrator_tol, 
  sapply(1:4, function(i) mean(diff_result_matrix$q_dot_q[diff_result_matrix$int_tol == integrator_tol[i]])),
  col = "red",
  pch = 15
)
segments(
  integrator_tol, 
  sapply(1:4, function(i) mean(diff_result_matrix$q_dot_q[diff_result_matrix$int_tol == integrator_tol[i]]) - 
           sd(diff_result_matrix$q_dot_q[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  integrator_tol, 
  sapply(1:4, function(i) mean(diff_result_matrix$q_dot_q[diff_result_matrix$int_tol == integrator_tol[i]]) + 
           sd(diff_result_matrix$q_dot_q[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  col = "red"
)
cap <- 0.25
segments(
  integrator_tol / (1 + cap), 
  sapply(1:4, function(i) mean(diff_result_matrix$q_dot_q[diff_result_matrix$int_tol == integrator_tol[i]]) - 
           sd(diff_result_matrix$q_dot_q[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  integrator_tol * (1 + cap), 
  sapply(1:4, function(i) mean(diff_result_matrix$q_dot_q[diff_result_matrix$int_tol == integrator_tol[i]]) - 
           sd(diff_result_matrix$q_dot_q[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  col = "red"
)
segments(
  integrator_tol / (1 + cap), 
  sapply(1:4, function(i) mean(diff_result_matrix$q_dot_q[diff_result_matrix$int_tol == integrator_tol[i]]) + 
           sd(diff_result_matrix$q_dot_q[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  integrator_tol * (1 + cap), 
  sapply(1:4, function(i) mean(diff_result_matrix$q_dot_q[diff_result_matrix$int_tol == integrator_tol[i]]) + 
           sd(diff_result_matrix$q_dot_q[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  col = "red"
)
abline(h = 0, lty = 2, col = "blue")

par(mfrow = c(1, 1))
