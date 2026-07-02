# n_eff comparison - Stan vs GRHMC - Longer run

rm(list = ls())

library(MASS)
library(brms)
library(ggplot2)
library(dplyr)
library(rstan)
library(foreach)
library(gridExtra)

n_iterations <- 10
n_samples_per_iteration_stan <- 100000
n_samples_per_iteration_grhmc <- 1000000

head(Boston)

x <- as.matrix(Boston[, -ncol(Boston)]) # 13 covariates
y <- Boston$medv
n_rows <- length(y)
length_beta <- ncol(x)

scale_x <- scale(x)
scale_y <- scale(y)
mean_y <- mean(y)

Boston_scale <- data.frame(scale_x, medv = scale_y)

# Stan

stan_run <- readRDS("disc_grad/example_scripts/section_6/regularized_linear_regression/stan_identity_run_example.RDS")
stan_samples <- rstan::extract(stan_run, permute = FALSE)
original_scale_stan_array <- array(dim = c(dim(stan_samples)[1:2], length_beta + 1))

for (i in 1:n_iterations) {
  
  store_matrix <- matrix(0, nrow = 100000, ncol = length_beta) # store the scaled beta samples 
  
  for (j in 1:length_beta) {
    
    current_beta_vec <- stan_samples[, i, 2 * length_beta + 1 + j]
    
    store_matrix[, j] <- current_beta_vec * sd(y) / sd(x[, j])
    
  }
  # First column: sigma scaled back to original scale, since sigma = exp(0.5 * gamma)
  original_scale_stan_array[, i, ] <- cbind(exp(0.5 * stan_samples[, i, 2 * length_beta + 1]) * sd(y), store_matrix)
}

original_scale_stan_monitor <- rstan::monitor(original_scale_stan_array, warmup = 0)
rownames(original_scale_stan_monitor)
n_eff_original_scale_stan <- original_scale_stan_monitor$n_eff
n_eff_original_scale_stan
n_evals_ode_original_scale_stan <- sum(rstan::get_num_leapfrog_per_iteration(stan_run))
n_evals_ode_original_scale_stan
n_eff_per_eval_ode_original_scale_stan <- n_eff_original_scale_stan / n_evals_ode_original_scale_stan
n_eff_per_eval_ode_original_scale_stan

cbind(round(original_scale_stan_monitor$Q50, 3), round(original_scale_stan_monitor$sd, 3))
cbind(round(n_eff_original_scale_stan / 10000, 3))
cbind(round(n_eff_original_scale_stan * 1000000 / n_evals_ode_original_scale_stan, 3))

original_scale_stan_matrix <- original_scale_stan_array
dim(original_scale_stan_matrix) <- c(n_samples_per_iteration_stan * n_iterations, length_beta + 1)
round(sapply(1:length_beta, function(i) mean(original_scale_stan_matrix[, i + 1] == 0)), 3) # Posterior P(beta = 0)

# GRHMC 

grhmc_disc_grad_run <- readRDS("disc_grad/example_scripts/section_6/regularized_linear_regression/disc_grad_run_identity_with_n_evals_ode_run_example.RDS")

original_scale_grhmc_disc_grad_array <- array(dim = c(n_samples_per_iteration_grhmc, n_iterations, length_beta + 1))

for (i in 1:n_iterations) {
  
  store_matrix <- matrix(0, nrow = n_samples_per_iteration_grhmc, ncol = length_beta) # store the scaled beta samples 
  
  for (j in 1:length_beta) {
    
    current_beta_vec <- pmax(0, grhmc_disc_grad_run[[i]]$q_original_samples[, j]) - pmax(0, grhmc_disc_grad_run[[i]]$q_original_samples[, length_beta + j])
    
    store_matrix[, j] <- current_beta_vec * sd(y) / sd(x[, j])
    
  }
  
  original_scale_grhmc_disc_grad_array[, i, ] <- cbind(exp(0.5 * grhmc_disc_grad_run[[i]]$q_original_samples[, 2 * length_beta + 1]) * sd(y), store_matrix)
}

original_scale_grhmc_disc_grad_monitor <- rstan::monitor(original_scale_grhmc_disc_grad_array, warmup = 0)
n_eff_original_scale_grhmc_disc_grad <- original_scale_grhmc_disc_grad_monitor$n_eff
n_eff_original_scale_grhmc_disc_grad
n_evals_ode_original_scale_grhmc_disc_grad <- 0
for (i in 1:n_iterations) {
  n_evals_ode_original_scale_grhmc_disc_grad <- n_evals_ode_original_scale_grhmc_disc_grad + 
    grhmc_disc_grad_run[[i]]$n_evals_ode
}
n_eff_per_eval_ode_original_scale_grhmc_disc_grad <- n_eff_original_scale_grhmc_disc_grad / n_evals_ode_original_scale_grhmc_disc_grad
n_eff_per_eval_ode_original_scale_grhmc_disc_grad

cbind(round(original_scale_grhmc_disc_grad_monitor$Q50, 3), round(original_scale_grhmc_disc_grad_monitor$sd, 3))
cbind(round(n_eff_original_scale_grhmc_disc_grad / 10000, 3))
cbind(round(n_eff_original_scale_grhmc_disc_grad * 1000000 / n_evals_ode_original_scale_grhmc_disc_grad, 3))

original_scale_grhmc_matrix <- original_scale_grhmc_disc_grad_array
dim(original_scale_grhmc_matrix) <- c(n_samples_per_iteration_grhmc * n_iterations, length_beta + 1)
round(sapply(1:length_beta, function(i) mean(original_scale_grhmc_matrix[, i + 1] == 0)), 3) # Posterior P(beta = 0)

