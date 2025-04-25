# Start at mle values, prior prob 0.95, zero vector m and s identity
# Use vectorized version of gradient of log target here. 

rm(list = ls())

library(MASS)
library(brms)
library(ggplot2)
library(dplyr)
library(rstan)
library(foreach)
library(gridExtra)

source("disc_grad/implementation_scripts/general_scripts/grhmc_discontinuous_gradient_transformed_function.R")

n_samples_per_iteration <- 10000
time_period <- 10000

head(Boston)

x <- as.matrix(Boston[, -ncol(Boston)]) # 13 covariates
y <- Boston$medv
n_rows <- length(y)
length_beta <- ncol(x)

scale_x <- scale(x)
scale_y <- scale(y)
mean_y <- mean(y)

Boston_scale <- data.frame(scale_x, medv = scale_y)
og_regmod <- lm(medv ~ ., data = Boston_scale)

# Fix prior Var(beta | beta != 0) to 1, tweak P(beta==0) by tweaking mu and rho.  
cond_prior_var <- 1
# prior_prob_vec <- seq(from = 0.005, to = 0.995, by = 0.005)
prior_prob_vec <- c(0.95)

solve_for_mu_rho_given_p <- function(p) { # p = prior probability of beta
  
  # Given the prior probability of beta and conditional variance --> want to solve for mu and rho
  
  mu_given_rho_and_p <- function(rho, p) { # use this equation to express mu w.r.t rho and p to have a one-dimensional root finding problem
    
    obj <- function(mu) {
      (-1 + 2 * pnorm(sqrt(2) / rho * mu * sqrt(2) / 2) - 1) ^ 2 / 4 - p
    }
    
    return(uniroot(obj, lower = -100, upper = 100, tol = 1e-12)$root)
    
  }
  
  cond_var_given_mu_and_rho <- function(mu, rho) { # the conditional prior variance based on mu and rho
    
    mu_star_1 <- rho * dnorm(mu / rho) + mu * pnorm(mu / rho)
    mu_star_2 <- rho * mu * dnorm(mu / rho) + (mu ^ 2 + rho ^ 2) * pnorm(mu / rho)
    uncond_var_beta <- 2 * (mu_star_2 - mu_star_1 ^ 2)
    cond_var_beta <- uncond_var_beta / (1 - pnorm(-mu / rho) ^ 2)
    return(cond_var_beta)
    
  }
  
  obj_for_rho <- function(rho) { # the root function for rho using the desired conditional prior variance
    cond_var_given_mu_and_rho(mu_given_rho_and_p(rho, p), rho) - cond_prior_var
  }
  
  rho_sol <- uniroot(obj_for_rho, lower = 0.01, upper = 10, tol = 1e-12)$root # solve for rho based on desired conditional prior variance
  mu_sol <- mu_given_rho_and_p(rho_sol, p) # use the result above to obtain the corresponding mu based on the desired prior probability
  
  mu_sol <- ifelse(abs(mu_sol) < 1e-13, 0, mu_sol)
  
  return(list(mu = mu_given_rho_and_p(rho_sol, p), rho = rho_sol))
  
}

mu_rho_sols <- solve_for_mu_rho_given_p(prior_prob_vec)

mu <- mu_rho_sols$mu
rho <- mu_rho_sols$rho

# Also fit using scaled x and y

length_beta <- ncol(x)

# beta_j = max(0, beta_j+) - max(0, beta_j-) with beta_j- and beta_j+ having N(mu, rho^2) as prior

# Also define gamma = log (sigma ^ 2) and set a standard normal distribution prior on gamma as well
# Default prior for sigma in stan_glm is Exp(k) where k is the rate parameter, this can be used here as well
k <- 1 

# q[1:10]: beta+, q[11:20]: beta-, q[21] = gamma = log sigma^2
# Need to introduce gamma as the we need a variable that is defined in the whole real line for this to work

A <- diag(2 * length_beta)
A <- cbind(A, rep(0, 2 * length_beta))
# dim(A)

grhmc_model_list <- list(
  
  n_parameters = 2 * length_beta + 1,
  
  grad_jump_fun = function(q) {
    
    grad_indices <<- 
      sapply(1:(2 * length_beta), function(i) as.integer(q[i] > 0))
    
  },
  
  region_non_lin_root_list = list(root_fun = NULL, event_fun = NULL),
  
  region_lin_root_list = list(A = A, B = rep(0, 2 * length_beta)),
  
  # sim_q0 = function() rnorm(2 * length_beta + 1)
  sim_q0 = function() c(rnorm(2 * length_beta, mean = mu, sd = rho), log(rexp(1, rate = k) ^ 2))
  
)

# grhmc_model_list$log_target_grad <- function(q) { # old version of log target gradient
#   
#   # If grad_indices = 0 for i = 1:length_beta --> max(0, beta_i+) = 0. Similar for beta_i-
#   beta <- q[1:length_beta] * grad_indices[1:length_beta] - q[(length_beta + 1):(2 * length_beta)] * grad_indices[(length_beta + 1):(2 * length_beta)]
#   y_minus_linear_predictor <- scale_y - scale_x %*% beta
#   exp_gamma <- exp(q[grhmc_model_list$n_parameters])
#   
#   log_target_grad_value_vec <- numeric(grhmc_model_list$n_parameters)
#   
#   for (i in 1:(length_beta)) {
#     
#     # If grad_indices = 0 for i = 1:length_beta --> beta_i+ < 0 so only the prior part in this region, no term related to observations and covariates
#     log_target_grad_value_vec[i] <- -(q[i] - mu) / rho ^ 2 + sum(y_minus_linear_predictor * scale_x[, i]) / exp_gamma * grad_indices[i]
#     
#   }
#   
#   for (i in (length_beta + 1):(2 * length_beta)) {
#     # Similar as above, negative sign for the term related to observations and covariates because pmax(0, beta+) - pmax(0, beta-)
#     log_target_grad_value_vec[i] <- -(q[i] - mu) / rho ^ 2 - sum(y_minus_linear_predictor * scale_x[, i - length_beta]) / exp_gamma * grad_indices[i]
#     
#   }
#   
#   log_target_grad_value_vec[grhmc_model_list$n_parameters] <- - n_rows / 2 - q[2 * length_beta + 1] + sum(y_minus_linear_predictor ^ 2) / (2 * exp_gamma) # if prior of gamma = log(sigma ^ 2) is N(0, 1)
#   # log_target_grad_value_vec[grhmc_model_list$n_parameters] <- - 1 / 2 * k * exp(q[2 * length_beta + 1] / 2) + 1 / 2 - n_rows / 2 + sum(y_minus_linear_predictor ^ 2) / (2 * exp_gamma) # gradient of gamma = log(sigma^2) if prior of sigma is Exp(1) # if prior of sigma = exp(0.5 * gamma) is Exp(1)
#   log_target_grad_value_vec
#   
# }

beta_plus_indices <- 1:(length_beta)
beta_minus_indices <- (length_beta + 1):(2 * length_beta)

grhmc_model_list$log_target_grad <- function(q) { # new version, vectorized and therefore slightly faster
  
  # If grad_indices = 0 for i = 1:length_beta --> max(0, beta_i+) = 0. Similar for beta_i-
  
  beta <- q[beta_plus_indices] * grad_indices[beta_plus_indices] - q[beta_minus_indices] * grad_indices[beta_minus_indices]
  y_minus_linear_predictor <- as.numeric(scale_y - (scale_x %*% beta))
  exp_gamma <- exp(q[grhmc_model_list$n_parameters])
  ratio_grad_indices_exp_gamma <- grad_indices / exp_gamma
  log_target_grad_value_vec <- numeric(grhmc_model_list$n_parameters)
  
  prior_beta_minus_plus_term <- -(q - mu) / (rho ^ 2)
  
  y_minus_linear_predictor_times_scale_x <- y_minus_linear_predictor %*% scale_x
  
  # If grad_indices = 0 for i = 1:length_beta --> beta_i+ < 0 so only the prior part in this region, no term related to observations and covariates
  log_target_grad_value_vec[beta_plus_indices] <- prior_beta_minus_plus_term[beta_plus_indices] + y_minus_linear_predictor_times_scale_x * ratio_grad_indices_exp_gamma[beta_plus_indices]
  # Similar as above, negative sign for the term related to observations and covariates because pmax(0, beta+) - pmax(0, beta-)
  log_target_grad_value_vec[beta_minus_indices] <- prior_beta_minus_plus_term[beta_minus_indices] - y_minus_linear_predictor_times_scale_x * ratio_grad_indices_exp_gamma[beta_minus_indices]
  log_target_grad_value_vec[grhmc_model_list$n_parameters] <- - (n_rows / 2) - q[2 * length_beta + 1] + (sum(y_minus_linear_predictor ^ 2) / (2 * exp_gamma)) # if prior of gamma = log(sigma ^ 2) is N(0, 1)
  
  log_target_grad_value_vec
  
}

set.seed(42)

# qbar_initial <- grhmc_model_list$sim_q0()
beta_plus_initial <- ifelse(og_regmod$coefficients[-1] < 0, mu, abs(og_regmod$coefficients[-1])) # start at MLE
beta_minus_initial <- ifelse(og_regmod$coefficients[-1] > 0, mu, abs(og_regmod$coefficients[-1]))
qbar_initial <- c(beta_plus_initial, beta_minus_initial, sum((og_regmod$residuals)^2) / og_regmod$df.residual)
pbar_initial <- rnorm(grhmc_model_list$n_parameters)

grad_indices <<- grhmc_model_list$grad_jump_fun(qbar_initial)

system.time(
  boston_scale_fit_grhmc <- grhmc_discontinuous_gradient_transformed_function(
    model_list = grhmc_model_list,
    lambda = 0.2,
    T = time_period,
    n_samples = n_samples_per_iteration,
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
    last.root.offset.lin.root.finder = 1.0e-8,
    last.root.offset.non.lin.root.finder = 1.0e-8,
    precision_real_root_lin_root_finder = 1.0e-13,
    num_subdiv_non_lin_root_finder = 8L,
    verbose_at_refresh = T
  )   
)

samples_matrix <- boston_scale_fit_grhmc$q_original_samples

boston_scale_grhmc_beta_store_matrix <- matrix(0, nrow = n_samples_per_iteration, ncol = length_beta) # store the scaled beta samples

for (j in 1:length_beta) {
  
  current_beta_vec <- pmax(0, samples_matrix[, j]) - pmax(0, samples_matrix[, length_beta + j])
  
  boston_scale_grhmc_beta_store_matrix[, j] <- current_beta_vec
  
}

boston_grhmc_ratio_samples_zero <- sapply(1:ncol(boston_scale_grhmc_beta_store_matrix), function(i) mean(boston_scale_grhmc_beta_store_matrix[, i] == 0)) # posterior probability of beta equal to zero

output <- c(boston_grhmc_ratio_samples_zero)

names(output) <- 
  c(colnames(x))

output
