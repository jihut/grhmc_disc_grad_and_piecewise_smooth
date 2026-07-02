rm(list = ls())

library(MASS)
library(brms)
library(ggplot2)
library(dplyr)
library(rstan)
library(foreach)
library(gridExtra)

n_cores <- 10
n_iterations <- 10
n_samples_per_iteration <- 100000
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

# Fix prior Var(beta | beta != 0) to 1, tweak P(beta==0) by tweaking mu and rho.  
cond_prior_var <- 1
# prior_prob_vec <- seq(from = 0.005, to = 0.995, by = 0.005)
prior_prob <- 0.95

solve_for_mu_rho_given_p <- function(p) { # p = prior probability of beta == 0
  
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

mu_rho_sols <- solve_for_mu_rho_given_p(prior_prob)

mu <- mu_rho_sols$mu
rho <- mu_rho_sols$rho

stan_model_specification <- function(mu, rho) {
  
  stan_code <- sprintf("
  
  data {
    int <lower=0> n;
    int <lower=0> p;
    vector[n] y;
    matrix[n, p] x;
  }

  parameters {
    vector[p] beta_plus;
    vector[p] beta_minus; 
    real gamma;
  }

  transformed parameters {
    vector[p] beta = fmax(0, beta_plus) - fmax(0, beta_minus);
    real <lower = 0> sigma = exp(0.5 * gamma);
  }

  model {
    real mu;
    real rho; 
    mu = %f;
    rho = %f;
    beta_plus ~ normal(mu, rho);
    beta_minus ~ normal(mu, rho);
    
    gamma ~ normal(0, 1);  
    
    y ~ normal(x * beta, sigma);
    
  }",
  mu,
  rho
  )
  
  stan_code
  
}

stan_code <- stan_model_specification(mu, rho)

stan_run <- rstan::stan(
  model_code = stan_code, 
  data = list(
    n = n_rows, 
    p = length_beta,
    y = as.numeric(scale_y),
    x = scale_x
  ),
  chains = n_iterations, 
  iter = 2 * n_samples_per_iteration,
  cores = n_cores / 2,
  seed = 42 # currently 42, try 420
)

# readRDS(stan_run, "disc_grad/example_scripts/section_6/regularized_linear_regression/stan_run_example.RDS")
stan_run <- readRDS("disc_grad/example_scripts/section_6/regularized_linear_regression/stan_run_example.RDS")

stan_samples <- rstan::extract(stan_run)

boston_scale_stan_beta_store_matrix <- matrix(0, nrow = n_samples_per_iteration * n_iterations, ncol = length_beta)

for (j in 1:length_beta) {
  
  current_beta_vec <- pmax(0, stan_samples$beta_plus[, j]) - pmax(0, stan_samples$beta_minus[, j])
  
  boston_scale_stan_beta_store_matrix[, j] <- current_beta_vec
  
}

boston_stan_ratio_samples_zero <- sapply(1:ncol(boston_scale_stan_beta_store_matrix), function(i) mean(boston_scale_stan_beta_store_matrix[, i] == 0))

boston_stan_ratio_samples_zero

cbind(colnames(x), boston_stan_ratio_samples_zero)
