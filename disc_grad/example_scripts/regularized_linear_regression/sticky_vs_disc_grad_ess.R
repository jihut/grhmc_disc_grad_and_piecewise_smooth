rm(list = ls())

library(MASS)
library(brms)
library(ggplot2)
library(dplyr)
library(rstan)
library(foreach)
library(gridExtra)

n_iterations <- 10
n_samples_per_iteration <- 1000000

head(Boston)

x <- as.matrix(Boston[, -ncol(Boston)]) # 13 covariates
y <- Boston$medv
n_rows <- length(y)
length_beta <- ncol(x)

scale_x <- scale(x)
scale_y <- scale(y)
mean_y <- mean(y)

disc_grad_runs <- readRDS("disc_grad/example_scripts/regularized_linear_regression/disc_grad_run_identity_with_n_evals_ode_run_example.RDS")
sticky_runs <- readRDS("disc_grad/example_scripts/regularized_linear_regression/sticky_run_identity_with_n_evals_ode_run_example.RDS")

# Discontinuous gradient setup

disc_grad_samples_array <- array(dim = c(n_samples_per_iteration, n_iterations, length_beta + 1))

for (i in 1:n_iterations) {
  
  current_samples_matrix <- matrix(nrow = n_samples_per_iteration, ncol = length_beta + 1) # last column for variance parameter
  
  for (j in 1:length_beta) {
    
    current_beta_vec <- pmax(0, disc_grad_runs[[i]]$q_original_samples[, j]) - pmax(0, disc_grad_runs[[i]]$q_original_samples[, length_beta + j])
    current_samples_matrix[, j] <- current_beta_vec
  }
  
  current_samples_matrix[, length_beta + 1] <- disc_grad_runs[[i]]$q_original_samples[, 2 * length_beta + 1]
  
  current_samples_matrix[, 1:length_beta] <- t(t(current_samples_matrix[, 1:length_beta]) / apply(x, 2, sd)) * sd(y) # scale back to original scale
  
  current_samples_matrix[, length_beta + 1] <- exp(0.5 * current_samples_matrix[, length_beta + 1]) * sd(y) # scale back to original scale and sigma
  
  disc_grad_samples_array[, i, ] <- current_samples_matrix
}

# Sticky setup

sticky_samples_array <- array(dim = c(n_samples_per_iteration, n_iterations, length_beta + 1))

for (i in 1:n_iterations) {
  
  current_samples_matrix <- sticky_runs[[i]]$q_original_samples
  
  current_samples_matrix[, 1:length_beta] <- t(t(current_samples_matrix[, 1:length_beta]) / apply(x, 2, sd)) * sd(y) # scale back to original scale
  
  current_samples_matrix[, length_beta + 1] <- exp(0.5 * current_samples_matrix[, length_beta + 1]) * sd(y) # scale back to original scale and sigma
  
  sticky_samples_array[, i, ] <- current_samples_matrix
  
}

# Comparison

disc_grad_monitor <- rstan::monitor(disc_grad_samples_array, warmup = 0)
disc_grad_monitor

sticky_monitor <- rstan::monitor(sticky_samples_array, warmup = 0)
sticky_monitor

sum_n_evals_ode_disc_grad <- sum(sapply(1:n_iterations, function(i) disc_grad_runs[[i]]$n_evals_ode))
sum_n_evals_ode_disc_grad

sum_n_evals_ode_sticky <- sum(sapply(1:n_iterations, function(i) sticky_runs[[i]]$n_evals_ode))
sum_n_evals_ode_sticky

cbind(disc_grad_monitor$n_eff, sticky_monitor$n_eff, disc_grad_monitor$Q50, sticky_monitor$Q50)
cbind(disc_grad_monitor$n_eff / sum_n_evals_ode_disc_grad * 1e6, sticky_monitor$n_eff / sum_n_evals_ode_sticky * 1e6)


