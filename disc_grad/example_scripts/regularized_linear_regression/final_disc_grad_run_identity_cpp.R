rm(list = ls())
library(doParallel)
library(MASS)
library(brms)
library(ggplot2)
library(dplyr)
library(rstan)
library(foreach)
library(gridExtra)

source("disc_grad/implementation_scripts/general_scripts/cpp/grhmc_discontinuous_gradient_transformed_with_relevant_indices_function_cpp.R")

grhmc_wrapper_file_path <- normalizePath(
  "disc_grad/implementation_scripts/general_scripts/cpp/grhmc_discontinuous_gradient_transformed_with_relevant_indices_function_cpp.R"
)
model_file_path <- normalizePath(
  "disc_grad/example_scripts/regularized_linear_regression/regularized_linear_model_transformed.cpp"
)
model_cache <- normalizePath(
  "disc_grad/example_scripts/regularized_linear_regression/regularized_linear_model_cache",
  mustWork = F
)
progress_log <- normalizePath("disc_grad/example_scripts/regularized_linear_regression/progress.log", mustWork = FALSE)
if (file.exists(progress_log)) file.remove(progress_log)
dir.create(model_cache, showWarnings = FALSE)
Rcpp::sourceCpp(model_file_path, cacheDir = model_cache, rebuild = T)

lambda <- 2.5
lambda_string <- gsub(".", "_", lambda, fixed = T)

trajectory_length <- 50000 # longer trajectory, but still less gradient evaluations than identity + tol = 1e-4
n_samples <- 1000000
n_trajectories <- 10
n_cores <- 10
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
  qbar_initial <- c(beta_plus_initial, beta_minus_initial, log(sum((og_regmod$residuals)^2) / og_regmod$df.residual))
  # qbar_initial <- c(rnorm(length_beta * 2, mean = mu, sd = rho), rnorm(1, mean = gamma_mu, sd = gamma_sigma))
  pbar_initial <- rnorm(length_beta * 2 + 1)
  u_initial <- rexp(1)
  
  test_run_warmup <- grhmc_discontinuous_gradient_function_cpp(
    model_list = function(settings) grhmc_linear_model_cpp(settings, mu = mu, rho = rho, gamma_mu = gamma_mu, gamma_sigma = gamma_sigma, x_matrix = scale_x, y_vector = scale_y),
    lambda = lambda,
    T = trajectory_length,
    n_samples = n_samples / 1000,
    qbar_initial,
    pbar_initial,
    u_initial, 
    s_elements = NULL,
    m_elements = NULL,
    random_state = NULL,
    rtol = 1e-4, # still using tol 1e-4 here
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
  
  new_pbar_initial <- rnorm(length_beta * 2 + 1)
  new_u_initial <- rexp(1)
  
  test_run_sampling <- grhmc_discontinuous_gradient_function_cpp(
    model_list = function(settings) grhmc_linear_model_cpp(settings, mu = mu, rho = rho, gamma_mu = gamma_mu, gamma_sigma = gamma_sigma, x_matrix = scale_x, y_vector = scale_y),
    lambda = lambda,
    T = trajectory_length,
    n_samples = n_samples,
    test_run_warmup$qbar_final,
    new_pbar_initial,
    new_u_initial, 
    s_elements = NULL,
    m_elements = NULL,
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

saveRDS(
  store_chains,
  "disc_grad/example_scripts/regularized_linear_regression/final_disc_grad_run_identity_cpp.RDS"
)

store_chains <- readRDS(
  "disc_grad/example_scripts/regularized_linear_regression/final_disc_grad_run_identity_cpp.RDS"
)

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

prob_zero_beta <- sapply(1:(length_beta + 1), function(i) mean(original_scale_grhmc_disc_grad_array[, , i] == 0))

median_beta <- sapply(1:(length_beta + 1), function(i) median(original_scale_grhmc_disc_grad_array[, , i]))

n_evals_ode_original_scale_grhmc_disc_grad <- 0
for (i in 1:n_trajectories) {
  n_evals_ode_original_scale_grhmc_disc_grad <- n_evals_ode_original_scale_grhmc_disc_grad + 
    store_chains[[i]]$n_evals_ode
}

original_scale_grhmc_se_mean_q <- sapply(
  1:(length_beta + 1),
  function(i) sd(
    sapply(1:n_trajectories, function(j) mean(original_scale_grhmc_disc_grad_array[, j, i]))
  )
)

original_scale_grhmc_se_prob_zero <- sapply(
  1:(length_beta + 1),
  function(i) sd(
    sapply(1:n_trajectories, function(j) mean(original_scale_grhmc_disc_grad_array[, j, i] == 0))
  )
)

original_scale_grhmc_samples <- # list of samples of q directly
  lapply(1:n_trajectories, function(i)
    (original_scale_grhmc_disc_grad_array[, i, ])) 

original_scale_grhmc_batch_ess_q <- 
  # size chosen by the plot from above
  sapply(seq_len(ncol(original_scale_grhmc_samples[[1]])), function(j) {
    stableGR::n.eff(lapply(original_scale_grhmc_samples, function(x)
      x[, j, drop = FALSE]),
      multivariate = FALSE,
      size = 100000)$n.eff / length(original_scale_grhmc_samples) 
  })

original_scale_grhmc_batch_ess_q * 1e6 / n_evals_ode_original_scale_grhmc_disc_grad

# Examining beta+ and beta-

grhmc_se_mean_q <- sapply(
  1:(2 * length_beta + 1),
  function(i) sd(
    sapply(1:n_trajectories, function(j) mean(store_chains[[j]]$q_original_samples[, i]))
  )
)

mean_beta_plus_minus <- colMeans(store_chains[[i]]$q_original_samples)

grhmc_samples <- # list of samples of q directly
  lapply(1:n_trajectories, function(i)
    (store_chains[[i]]$q_original_samples)) 

grhmc_batch_ess_q <- 
  # size chosen by the plot from above
  sapply(seq_len(ncol(grhmc_samples[[1]])), function(j) {
    stableGR::n.eff(lapply(grhmc_samples, function(x)
      x[, j, drop = FALSE]),
      multivariate = FALSE,
      size = 100000)$n.eff / length(grhmc_samples) 
  })
batch_ess_beta_plus_minus_matrix <- rbind(
  batch_ess_beta_plus_minus_matrix,
  grhmc_batch_ess_q
)

batch_ess_beta_plus_minus_per_million_eval_matrix <- rbind(
  batch_ess_beta_plus_minus_per_million_eval_matrix,
  grhmc_batch_ess_q * 1e6 / n_evals_ode_original_scale_grhmc_disc_grad
)

n_boundary_crossings <- sapply(
  1:n_trajectories,
  function(i) store_chains[[i]]$n_boundary_crossings
)
