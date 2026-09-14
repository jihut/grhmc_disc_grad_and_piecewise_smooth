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

dim_vec <- 2 ^ {1:11}
n_chains <- 100
n_samples <- 10000

# set.seed(42)
# q1 <- rnorm(1e8)
# q2_given_q1 <- rnorm(1e8, mean = pmax(0, q1), sd = 1)
# mean(q2_given_q1)

true_mean_q2 <- 
  integrate(function(x) x * true_pdf_q2(x), lower = -Inf, upper = Inf, rel.tol = 1e-13, abs.tol = 1e-13)
true_mean_q2
true_mean_q2_squared <- 
  integrate(function(x) (x ^ 2) * true_pdf_q2(x), lower = -Inf, upper = Inf, rel.tol = 1e-13, abs.tol = 1e-13)
true_mean_q2_squared

integrand_true_mean_q1_squared_plus_q2_squared <- function(q1, q2) {
  (q1 ^ 2 + q2 ^ 2) * dnorm(q1, mean = 0, sd = 1) * dnorm(q2, mean = max(0, q1), sd = 1)
}

true_mean_q1_squared_plus_q2_squared  <- integrate(function(q1) {
  sapply(q1, function(q1_star) {
    integrate(function(q2) integrand_true_mean_q1_squared_plus_q2_squared(q1_star, q2), lower = -Inf, upper = Inf, rel.tol = 1e-13, abs.tol = 1e-13)$value
  })
}, lower = -Inf, upper = Inf, rel.tol = 1e-13, abs.tol = 1e-13)

true_mean_q1_squared_plus_q2_squared

estimated_quantities <- function(samples_matrix, calculate_wasserstein_2 = FALSE) {
  
  mean_q2 <- mean(samples_matrix[, 2])
  mean_q2_squared <- mean(samples_matrix[, 2] ^ 2)
  mean_q1_squared_plus_q2_squared <- mean(samples_matrix[, 1] ^ 2 + samples_matrix[, 2] ^ 2)
  
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
        q2 = mean_q2, q2_squared = mean_q2_squared, q1_squared_plus_q2_squared = mean_q1_squared_plus_q2_squared, wasserstein_2 = sqrt(wasserstein_2_squared_value)
      )
    )
    
  } else {
    return(
      c(
        q2 = mean_q2, q2_squared = mean_q2_squared, q1_squared_plus_q2_squared = mean_q1_squared_plus_q2_squared
      )
    )  
  }
  
}

for (k in 10:length(dim_vec)) {
  print(k)
  n_dim <- dim_vec[k]
  stan_run <- readRDS(
    paste0("disc_grad/example_scripts/toy_model/stan/max_model_d_", n_dim, "_stan.RDS")
  )
  # if (k >= 10) {
  #   stan_samples <- rstan::extract(stan_run, pars = c("q1", "q2"), permute = F)
  # } else {
  #   stan_samples <- rstan::extract(stan_run, pars = c("q[1]", "q[2]"), permute = F)
  # }
  stan_samples <- rstan::extract(stan_run, pars = c("q[1]", "q[2]"), permute = F)
  relevant_matrix <- cbind(t(apply(stan_samples, 2, estimated_quantities)), 
                           dim = n_dim)
  individual_ess_matrix <- matrix(nrow = n_chains, ncol = 4)
  individual_ess_matrix[, 4] <- n_dim
  for (l in 1:n_chains) {
    individual_samples_array <- array(dim = c(n_samples, 1, 2)) 
    individual_samples_array[, 1, ] <- stan_samples[, l, 1:2]
    individual_grhmc_monitor <- rstan::monitor(individual_samples_array, warmup = 0)
    individual_ess_matrix[l, 1:2] <- individual_grhmc_monitor$n_eff[1:2]
    individual_ess_matrix[l, 3] <- sum(rstan::get_sampler_params(stan_run, inc_warmup = FALSE)[[l]][, "n_leapfrog__"])
  }
  stan_monitor <- rstan::monitor(stan_samples[, , 1:2], warmup = 0)
  ess_matrix <- c(q = stan_monitor$n_eff[1:2], n_evals_ode = sum(rstan::get_num_leapfrog_per_iteration(stan_run)))
  if (k == 1) {
    final_result_matrix <- relevant_matrix
    final_ess_matrix <- ess_matrix
    final_individual_ess_matrix <- individual_ess_matrix
  } else {
    final_result_matrix <- rbind(final_result_matrix, relevant_matrix)
    final_ess_matrix <- rbind(final_ess_matrix, ess_matrix)
    final_individual_ess_matrix <- rbind(final_individual_ess_matrix, individual_ess_matrix)
  }
}

rm(stan_run)
gc()

final_result_matrix <- as.data.frame(final_result_matrix)
final_ess_matrix <- as.data.frame(final_ess_matrix)
final_individual_ess_matrix <- as.data.frame(final_individual_ess_matrix)
colnames(final_individual_ess_matrix) <- colnames(final_ess_matrix)

diff_result_matrix <- 
  sweep(final_result_matrix, 2, c(true_mean_q2$value, true_mean_q2_squared$value, true_mean_q1_squared_plus_q2_squared$value, 0), "-")

saveRDS(list(final_result_matrix = final_result_matrix, diff_result_matrix = diff_result_matrix, final_ess_matrix = final_ess_matrix, final_individual_ess_matrix = final_individual_ess_matrix),
        "disc_grad/example_scripts/toy_model/max_model_varying_dimensions_stan_results.RDS")
stan_results <- readRDS("disc_grad/example_scripts/toy_model/max_model/max_model_varying_dimensions_stan_results.RDS")
final_result_matrix <- stan_results$final_result_matrix
final_ess_matrix <- stan_results$final_ess_matrix
diff_result_matrix <- stan_results$diff_result_matrix

# Plot of difference between estimates of E(q_2) and true E(q_2) 
par(mfrow = c(1, 3))

plot(
  # exp(jitter(log2(diff_result_matrix$dim), amount = 0.1)), 
  (jitter(log2(diff_result_matrix$dim), amount = 0.25)),
  diff_result_matrix$q2, 
  # log = "x",
  xaxt = "n",
  xlab = "Dimension",
  ylab = "Deviation from true value",
  main = expression(E(q[2]))
)

axis(1, at = log2(dim_vec), labels = dim_vec)
points(
  log2(dim_vec), 
  sapply(1:length(dim_vec), function(i) mean(diff_result_matrix$q2[diff_result_matrix$dim == dim_vec[i]])),
  col = "red",
  pch = 15
)
segments(
  log2(dim_vec), 
  sapply(1:length(dim_vec), function(i) mean(diff_result_matrix$q2[diff_result_matrix$dim == dim_vec[i]]) - 
           sd(diff_result_matrix$q2[diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  log2(dim_vec), 
  sapply(1:length(dim_vec), function(i) mean(diff_result_matrix$q2[diff_result_matrix$dim == dim_vec[i]]) + 
           sd(diff_result_matrix$q2[diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  col = "red"
)
cap <- 0.25
segments(
  log2(dim_vec) - 0.1, 
  sapply(1:length(dim_vec), function(i) mean(diff_result_matrix$q2[diff_result_matrix$dim == dim_vec[i]]) - 
           sd(diff_result_matrix$q2[diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  log2(dim_vec) + 0.1, 
  sapply(1:length(dim_vec), function(i) mean(diff_result_matrix$q2[diff_result_matrix$dim == dim_vec[i]]) - 
           sd(diff_result_matrix$q2[diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  col = "red"
)
segments(
  log2(dim_vec) - 0.1, 
  sapply(1:length(dim_vec), function(i) mean(diff_result_matrix$q2[diff_result_matrix$dim == dim_vec[i]]) + 
           sd(diff_result_matrix$q2[diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  log2(dim_vec) + 0.1, 
  sapply(1:length(dim_vec), function(i) mean(diff_result_matrix$q2[diff_result_matrix$dim == dim_vec[i]]) + 
           sd(diff_result_matrix$q2[diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  col = "red"
)
abline(h = 0, lty = 2, col = "blue")

# Plot of difference between estimates of E(q_2 ^ 2) and true E(q_2 ^ 2) 

plot(
  # exp(jitter(log2(diff_result_matrix$dim), amount = 0.1)), 
  (jitter(log2(diff_result_matrix$dim), amount = 0.25)),
  diff_result_matrix$q2_squared, 
  # log = "x",
  xaxt = "n",
  xlab = "Dimension",
  ylab = "Deviation from true value",
  main = expression(E(q[2]^2))
)

axis(1, at = log2(dim_vec), labels = dim_vec)
points(
  log2(dim_vec), 
  sapply(1:length(dim_vec), function(i) mean(diff_result_matrix$q2_squared[diff_result_matrix$dim == dim_vec[i]])),
  col = "red",
  pch = 15
)
segments(
  log2(dim_vec), 
  sapply(1:length(dim_vec), function(i) mean(diff_result_matrix$q2_squared[diff_result_matrix$dim == dim_vec[i]]) - 
           sd(diff_result_matrix$q2_squared[diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  log2(dim_vec), 
  sapply(1:length(dim_vec), function(i) mean(diff_result_matrix$q2_squared[diff_result_matrix$dim == dim_vec[i]]) + 
           sd(diff_result_matrix$q2_squared[diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  col = "red"
)
cap <- 0.25
segments(
  log2(dim_vec) - 0.1, 
  sapply(1:length(dim_vec), function(i) mean(diff_result_matrix$q2_squared[diff_result_matrix$dim == dim_vec[i]]) - 
           sd(diff_result_matrix$q2_squared[diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  log2(dim_vec) + 0.1, 
  sapply(1:length(dim_vec), function(i) mean(diff_result_matrix$q2_squared[diff_result_matrix$dim == dim_vec[i]]) - 
           sd(diff_result_matrix$q2_squared[diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  col = "red"
)
segments(
  log2(dim_vec) - 0.1, 
  sapply(1:length(dim_vec), function(i) mean(diff_result_matrix$q2_squared[diff_result_matrix$dim == dim_vec[i]]) + 
           sd(diff_result_matrix$q2_squared[diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  log2(dim_vec) + 0.1, 
  sapply(1:length(dim_vec), function(i) mean(diff_result_matrix$q2_squared[diff_result_matrix$dim == dim_vec[i]]) + 
           sd(diff_result_matrix$q2_squared[diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  col = "red"
)
abline(h = 0, lty = 2, col = "blue")

# Plot of difference between estimates of E(q dot q) and true E(q dot q) 

plot(
  # exp(jitter(log2(diff_result_matrix$dim), amount = 0.1)), 
  (jitter(log2(diff_result_matrix$dim), amount = 0.25)),
  diff_result_matrix$q1_squared_plus_q2_squared, 
  # log = "x",
  xaxt = "n",
  xlab = "Dimension",
  ylab = "Deviation from true value",
  main = expression(E(q[1]^2 + q[2]^2))
)

axis(1, at = log2(dim_vec), labels = dim_vec)
points(
  log2(dim_vec), 
  sapply(1:length(dim_vec), function(i) mean(diff_result_matrix$q1_squared_plus_q2_squared[diff_result_matrix$dim == dim_vec[i]])),
  col = "red",
  pch = 15
)
segments(
  log2(dim_vec), 
  sapply(1:length(dim_vec), function(i) mean(diff_result_matrix$q1_squared_plus_q2_squared[diff_result_matrix$dim == dim_vec[i]]) - 
           sd(diff_result_matrix$q1_squared_plus_q2_squared[diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  log2(dim_vec), 
  sapply(1:length(dim_vec), function(i) mean(diff_result_matrix$q1_squared_plus_q2_squared[diff_result_matrix$dim == dim_vec[i]]) + 
           sd(diff_result_matrix$q1_squared_plus_q2_squared[diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  col = "red"
)
cap <- 0.25
segments(
  log2(dim_vec) - 0.1, 
  sapply(1:length(dim_vec), function(i) mean(diff_result_matrix$q1_squared_plus_q2_squared[diff_result_matrix$dim == dim_vec[i]]) - 
           sd(diff_result_matrix$q1_squared_plus_q2_squared[diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  log2(dim_vec) + 0.1, 
  sapply(1:length(dim_vec), function(i) mean(diff_result_matrix$q1_squared_plus_q2_squared[diff_result_matrix$dim == dim_vec[i]]) - 
           sd(diff_result_matrix$q1_squared_plus_q2_squared[diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  col = "red"
)
segments(
  log2(dim_vec) - 0.1, 
  sapply(1:length(dim_vec), function(i) mean(diff_result_matrix$q1_squared_plus_q2_squared[diff_result_matrix$dim == dim_vec[i]]) + 
           sd(diff_result_matrix$q1_squared_plus_q2_squared[diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  log2(dim_vec) + 0.1, 
  sapply(1:length(dim_vec), function(i) mean(diff_result_matrix$q1_squared_plus_q2_squared[diff_result_matrix$dim == dim_vec[i]]) + 
           sd(diff_result_matrix$q1_squared_plus_q2_squared[diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  col = "red"
)
abline(h = 0, lty = 2, col = "blue")

par(mfrow = c(1, 1))

