rm(list = ls())
library(doParallel)
library(MASS)
library(brms)
library(ggplot2)
library(dplyr)
library(rstan)
library(foreach)
library(gridExtra)

source("disc_grad/implementation_scripts/general_scripts/cpp/grhmc_discontinuous_gradient_with_relevant_indices_function_cpp.R")

grhmc_wrapper_file_path <- normalizePath(
  "disc_grad/implementation_scripts/general_scripts/cpp/grhmc_discontinuous_gradient_with_relevant_indices_function_cpp.R"
)
model_file_path <- normalizePath(
  "disc_grad/example_scripts/regularized_linear_regression/regularized_linear_model.cpp"
)
model_cache <- normalizePath(
  "disc_grad/example_scripts/regularized_linear_regression/regularized_linear_model_cache",
  mustWork = F
)
progress_log <- normalizePath("disc_grad/example_scripts/regularized_linear_regression/progress.log", mustWork = FALSE)
if (file.exists(progress_log)) file.remove(progress_log)
dir.create(model_cache, showWarnings = FALSE)
Rcpp::sourceCpp(model_file_path, cacheDir = model_cache, rebuild = T)

trajectory_length <- 100000
n_samples <- 1000000
n_trajectories <- 10
n_cores <- 5
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
gamma_mu <- 0
gamma_sigma <- 1

init_cluster <- parallel::makeCluster(n_cores, outfile = progress_log)
parallel::clusterExport(init_cluster, c("grhmc_wrapper_file_path", "model_file_path", "model_cache"))
for (k in seq_along(init_cluster)) {
  parallel::clusterEvalQ(init_cluster[k], { # load cpp build one worker by one to avoid crash
    library(Rcpp)
    source(grhmc_wrapper_file_path)
    Rcpp::sourceCpp(model_file_path, cacheDir = model_cache, rebuild = FALSE)
    NULL
  })
}
doParallel::registerDoParallel(init_cluster)
doRNG::registerDoRNG(seed = 42)

store_chains <- foreach::foreach(
  i = 1:n_trajectories, 
  .noexport = c("grhmc_linear_model_cpp",
                "grhmc_discontinuous_gradient_function_cpp")
) %dopar% {
  og_regmod <- lm(medv ~ ., data = Boston_scale)
  
  beta_plus_initial <- ifelse(og_regmod$coefficients[-1] < 0, mu, abs(og_regmod$coefficients[-1])) # start at MLE
  beta_minus_initial <- ifelse(og_regmod$coefficients[-1] > 0, mu, abs(og_regmod$coefficients[-1]))
  q_initial <- c(beta_plus_initial, beta_minus_initial, sum((og_regmod$residuals)^2) / og_regmod$df.residual)
  p_initial <- rnorm(length_beta * 2 + 1)
  u_initial <- rexp(1)
  
  test_run_warmup <- grhmc_discontinuous_gradient_function_cpp(
    model_list = function(settings) grhmc_linear_model_cpp(settings, mu = mu, rho = rho, gamma_mu = gamma_mu, gamma_sigma = gamma_sigma, x_matrix = scale_x, y_vector = scale_y),
    lambda = 0.2,
    T = trajectory_length,
    n_samples = n_samples,
    q_initial,
    p_initial,
    u_initial, 
    random_state = NULL,
    rtol = 1e-4,
    atol = 1e-4,
    h_max = 1.0,
    last_root_offset_lin_root_finder = 1e-10,
    precision_real_root_lin_root_finder = 1.0e-13,
    coeff_tol_real_root_lin_root_finder = 1.0e-13,
    push_eps_tol = 1e-4,
    relevant_indices = NULL,
    return_output_from_ode = FALSE,
    trajectory_info_string = paste0("warmup ", i),
    progress_ratio = 0.05
  )
  
  new_p_initial <- rnorm(length_beta * 2 + 1)
  new_u_initial <- rexp(1)
  
  test_run_sampling <- grhmc_discontinuous_gradient_function_cpp(
    model_list = function(settings) grhmc_linear_model_cpp(settings, mu = mu, rho = rho, gamma_mu = gamma_mu, gamma_sigma = gamma_sigma, x_matrix = scale_x, y_vector = scale_y),
    lambda = 0.2,
    T = trajectory_length,
    n_samples = n_samples,
    test_run_warmup$q_final,
    new_p_initial,
    new_u_initial, 
    random_state = NULL,
    rtol = 1e-4,
    atol = 1e-4,
    h_max = 1.0,
    last_root_offset_lin_root_finder = 1e-10,
    precision_real_root_lin_root_finder = 1.0e-13,
    coeff_tol_real_root_lin_root_finder = 1.0e-13,
    push_eps_tol = 1e-4,
    relevant_indices = NULL,
    return_output_from_ode = FALSE,
    trajectory_info_string = paste0("sampling ", i),
    progress_ratio = 0.05
  )
  test_run_sampling
}
parallel::stopCluster(init_cluster)

##################################

parameter_names <- c("sigma", colnames(x))

original_scale_grhmc_disc_grad_array <- array(dim = c(n_samples, n_trajectories, length_beta + 1))

for (i in 1:n_trajectories) {
  
  store_matrix <- matrix(0, nrow = n_samples, ncol = length_beta) # store the scaled beta samples 
  
  for (j in 1:length_beta) {
    
    current_beta_vec <- pmax(0, store_chains[[i]]$q_original_samples[, j]) - pmax(0, store_chains[[i]]$q_original_samples[, length_beta + j])
    
    store_matrix[, j] <- current_beta_vec * sd(y) / sd(x[, j])
    
  }
  
  original_scale_grhmc_disc_grad_array[, i, ] <- cbind(exp(0.5 * store_chains[[i]]$q_original_samples[, 2 * length_beta + 1]) * sd(y), store_matrix)
}

original_scale_grhmc_disc_grad_monitor <- rstan::monitor(original_scale_grhmc_disc_grad_array, warmup = 0)

n_evals_ode_original_scale_grhmc_disc_grad <- 0
for (i in 1:n_trajectories) {
  n_evals_ode_original_scale_grhmc_disc_grad <- n_evals_ode_original_scale_grhmc_disc_grad + 
    store_chains[[i]]$n_evals_ode
}
n_eff_per_eval_ode_original_scale_grhmc_disc_grad <- n_eff_original_scale_grhmc_disc_grad / n_evals_ode_original_scale_grhmc_disc_grad
n_eff_per_eval_ode_original_scale_grhmc_disc_grad
n_eff_per_eval_ode_original_scale_grhmc_disc_grad * 1e6

final_grhmc_individual_rstan_monitor <- NULL

for (i in 1:n_trajectories) {
  grhmc_individual_samples_array <- array(dim = c(n_samples, 1, length_beta + 1))
  grhmc_individual_samples_array[, 1, ] <- original_scale_grhmc_disc_grad_array[, i, ]
  grhmc_individual_rstan_monitor <- rstan::monitor(grhmc_individual_samples_array, warmup = 0)
  rownames(grhmc_individual_rstan_monitor) <- parameter_names
  final_grhmc_individual_rstan_monitor <- rbind(final_grhmc_individual_rstan_monitor, grhmc_individual_rstan_monitor)
}

final_grhmc_individual_rstan_monitor <- as.data.frame(final_grhmc_individual_rstan_monitor)
final_grhmc_individual_rstan_monitor

