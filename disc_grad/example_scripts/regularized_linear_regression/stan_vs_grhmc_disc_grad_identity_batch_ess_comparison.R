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
parameter_names <- c("sigma", colnames(x))

# GRHMC

grhmc_disc_grad_run <- readRDS("disc_grad/example_scripts/regularized_linear_regression/disc_grad_run_identity_with_n_evals_ode_run_example.RDS")

original_scale_grhmc_disc_grad_array <- array(dim = c(n_samples_per_iteration_grhmc, n_iterations, length_beta + 1))

for (i in 1:n_iterations) {
  
  store_matrix <- matrix(0, nrow = n_samples_per_iteration_grhmc, ncol = length_beta) # store the scaled beta samples 
  
  for (j in 1:length_beta) {
    
    current_beta_vec <- pmax(0, grhmc_disc_grad_run[[i]]$q_original_samples[, j]) - pmax(0, grhmc_disc_grad_run[[i]]$q_original_samples[, length_beta + j])
    
    store_matrix[, j] <- current_beta_vec * sd(y) / sd(x[, j])
    
  }
  
  original_scale_grhmc_disc_grad_array[, i, ] <- cbind(exp(0.5 * grhmc_disc_grad_run[[i]]$q_original_samples[, 2 * length_beta + 1]) * sd(y), store_matrix)
}

# First: Preliminary check if the automated chosen batch size (i.e. size = NULL) in stableGR is adequate.
# If the output batch ESS is similar to the value where the batch ESS flattens out with respect to the batch size, then the default is used.
# Otherwise, the output batch ESS is too optimistic and a manual batch size value is chosen close to the elbow of the ESS vs batch size plot or the first part where the plot is flat.

n_evals_ode_original_scale_grhmc_disc_grad <-
  sum(sapply(1:n_iterations, function(i) grhmc_disc_grad_run[[i]]$n_evals_ode))

grhmc_samples <- # list of samples of q directly
  lapply(1:n_iterations, function(i)
    (original_scale_grhmc_disc_grad_array[, i, ])) 

n_splits <- 1000 # want to check for 1000 different batch sizes
n_cores <- 10
batch_size_grid <- seq(from = n_samples_per_iteration_grhmc / n_splits, to = n_samples_per_iteration_grhmc, by = n_samples_per_iteration_grhmc / n_splits)
log_path <- "disc_grad/example_scripts/regularized_linear_regression/batch_ess_progress.log"
if (file.exists(log_path)) file.remove(log_path)
init_cluster <- parallel::makeCluster(n_cores, outfile = log_path)
doParallel::registerDoParallel(init_cluster)

store_uni_ess_q <- foreach::foreach(i = 1:n_splits, .combine = rbind) %dopar% {
  print(i)
  sapply(seq_len(ncol(grhmc_samples[[1]])), function(j) {
    stableGR::n.eff(lapply(grhmc_samples, function(x)
      x[, j, drop = FALSE]),
      multivariate = FALSE,
      size = batch_size_grid[i])$n.eff / length(grhmc_samples)
  })
}
default_uni_ess_q <-
  sapply(seq_len(ncol(grhmc_samples[[1]])), function(j) {
    stableGR::n.eff(lapply(grhmc_samples, function(x)
      x[, j, drop = FALSE]),
      multivariate = FALSE,
      size = NULL)$n.eff / length(grhmc_samples)
  })
final_store_uni_ess_q <- store_uni_ess_q
final_default_uni_ess_q <- default_uni_ess_q

saveRDS(
  list(
    final_store_uni_ess_q = final_store_uni_ess_q,
    final_default_uni_ess_q = final_default_uni_ess_q
  ),
  "disc_grad/example_scripts/regularized_linear_regression/disc_grad_run_identity_batch_ess_analysis.RDS"
)

batch_ESS_size_analysis <- readRDS(
  "disc_grad/example_scripts/regularized_linear_regression/disc_grad_run_identity_batch_ess_analysis.RDS"
)
final_store_uni_ess_q <- batch_ESS_size_analysis$final_store_uni_ess_q
final_default_uni_ess_q <- batch_ESS_size_analysis$final_default_uni_ess_q

for (i in 1:(length_beta + 1)) {
  print(i)
  plot(batch_size_grid, final_store_uni_ess_q[, i], log = "y", main = paste0("i = ", i, ": automated batch size ESS = ", round(final_default_uni_ess_q[i], 1)))
}

# In this case, the batch ESS given by the automated choice of batch size is very large compared to the interval of medium batch sizes between 50000 and 200000 for all coordinates.
# It looks like choosing a batch size between the mentioned interval will give similar answers. 
# Therefore, 100000 is chosen here for all scenarios. 

# Second: Create summary outputs

final_grhmc_individual_rstan_monitor <- NULL
final_individual_result_matrix <- NULL

for (j in 1:n_iterations) {
  print(j)
  grhmc_individual_samples_array <- array(dim = c(n_samples_per_iteration_grhmc, 1, length_beta + 1))
  grhmc_individual_samples_array[, 1, ] <- original_scale_grhmc_disc_grad_array[, j, ]
  grhmc_individual_rstan_monitor <- rstan::monitor(grhmc_individual_samples_array, warmup = 0)
  rownames(grhmc_individual_rstan_monitor) <- parameter_names
  final_grhmc_individual_rstan_monitor <- rbind(final_grhmc_individual_rstan_monitor, grhmc_individual_rstan_monitor)
  individual_result_matrix <- 
    data.frame(
      q = colMeans(original_scale_grhmc_disc_grad_array[, j, ])
    )
  final_individual_result_matrix <- rbind(final_individual_result_matrix, individual_result_matrix)
}
final_grhmc_individual_rstan_monitor <- as.data.frame(final_grhmc_individual_rstan_monitor)
final_grhmc_individual_rstan_monitor

final_grhmc_batch_ess <- NULL

mean_q <- sapply(
  1:(length_beta + 1),
  function(i) 
    mean(final_grhmc_individual_rstan_monitor[grep(parameter_names[i], rownames(final_grhmc_individual_rstan_monitor)), ]$mean)
)
se_mean_q <- sapply(
  1:(length_beta + 1),
  function(i) 
    sd(final_grhmc_individual_rstan_monitor[grep(parameter_names[i], rownames(final_grhmc_individual_rstan_monitor)), ]$mean)
)
sd_q <- sapply(
  1:(length_beta + 1),
  function(i)
    sd(original_scale_grhmc_disc_grad_array[, , i])
)
mean_original_rstan_se_mean <- sapply(
  1:(length_beta + 1),
  function(i)
    sqrt(mean(final_grhmc_individual_rstan_monitor[grep(parameter_names[i], rownames(final_grhmc_individual_rstan_monitor)), ]$se_mean ^ 2))
)
batch_ess_q <- 
  # size chosen by the plot from above
  sapply(seq_len(ncol(grhmc_samples[[1]])), function(j) {
    stableGR::n.eff(lapply(grhmc_samples, function(x)
      x[, j, drop = FALSE]),
      multivariate = FALSE,
      size = 100000)$n.eff / length(grhmc_samples) 
  })
batch_ess_q_per_mill_eval_ode <- batch_ess_q / n_evals_ode_original_scale_grhmc_disc_grad * 1e6 
final_grhmc_batch_ess <- 
  data.frame(
    mean_q = mean_q,
    sd_q = sd_q,
    se_mean_q = se_mean_q,
    mean_original_rstan_se_mean = mean_original_rstan_se_mean,
    batch_ess_q = batch_ess_q,
    batch_ess_q_per_mill_eval_ode = batch_ess_q_per_mill_eval_ode,
    n_evals_ode = n_evals_ode_original_scale_grhmc_disc_grad, 
    parameter_names = parameter_names
  )

saveRDS(
  final_grhmc_batch_ess,
  "disc_grad/example_scripts/regularized_linear_regression/disc_grad_run_identity_batch_ess_results.RDS"
)

final_grhmc_batch_ess <- readRDS(
  "disc_grad/example_scripts/regularized_linear_regression/disc_grad_run_identity_batch_ess_results.RDS"
)

# Stan

stan_run <- readRDS("disc_grad/example_scripts/regularized_linear_regression/stan_identity_run_example.RDS")
stan_samples <- rstan::extract(stan_run, permute = FALSE)
original_scale_stan_array <- array(dim = c(dim(stan_samples)[1:2], length_beta + 1))

for (i in 1:n_iterations) {
  
  store_matrix <- matrix(0, nrow = n_samples_per_iteration_stan, ncol = length_beta) # store the scaled beta samples 
  
  for (j in 1:length_beta) {
    
    current_beta_vec <- stan_samples[, i, 2 * length_beta + 1 + j]
    
    store_matrix[, j] <- current_beta_vec * sd(y) / sd(x[, j])
    
  }
  # First column: sigma scaled back to original scale, since sigma = exp(0.5 * gamma)
  original_scale_stan_array[, i, ] <- cbind(exp(0.5 * stan_samples[, i, 2 * length_beta + 1]) * sd(y), store_matrix)
}

# First: Preliminary check if the automated chosen batch size (i.e. size = NULL) in stableGR is adequate.
# If the output batch ESS is similar to the value where the batch ESS flattens out with respect to the batch size, then the default is used.
# Otherwise, the output batch ESS is too optimistic and a manual batch size value is chosen close to the elbow of the ESS vs batch size plot or the first part where the plot is flat.

n_evals_ode_original_scale_stan <-
  sum(rstan::get_num_leapfrog_per_iteration(stan_run))

stan_samples <- # list of samples of q directly
  lapply(1:n_iterations, function(i)
    (original_scale_stan_array[, i, ])) 

n_splits <- 1000 # want to check for 1000 different batch sizes
n_cores <- 5
batch_size_grid <- seq(from = n_samples_per_iteration_stan / n_splits, to = n_samples_per_iteration_stan, by = n_samples_per_iteration_stan / n_splits)
log_path <- "disc_grad/example_scripts/regularized_linear_regression/batch_ess_progress.log"
if (file.exists(log_path)) file.remove(log_path)
init_cluster <- parallel::makeCluster(n_cores, outfile = log_path)
doParallel::registerDoParallel(init_cluster)

store_uni_ess_q <- foreach::foreach(i = 1:n_splits, .combine = rbind) %dopar% {
  print(i)
  sapply(seq_len(ncol(stan_samples[[1]])), function(j) {
    stableGR::n.eff(lapply(stan_samples, function(x)
      x[, j, drop = FALSE]),
      multivariate = FALSE,
      size = batch_size_grid[i])$n.eff / length(stan_samples)
  })
}
default_uni_ess_q <-
  sapply(seq_len(ncol(stan_samples[[1]])), function(j) {
    stableGR::n.eff(lapply(stan_samples, function(x)
      x[, j, drop = FALSE]),
      multivariate = FALSE,
      size = NULL)$n.eff / length(stan_samples)
  })
final_store_uni_ess_q <- store_uni_ess_q
final_default_uni_ess_q <- default_uni_ess_q

saveRDS(
  list(
    final_store_uni_ess_q = final_store_uni_ess_q,
    final_default_uni_ess_q = final_default_uni_ess_q
  ),
  "disc_grad/example_scripts/regularized_linear_regression/stan_batch_ess_analysis.RDS"
)

batch_ESS_size_analysis <- readRDS(
  "disc_grad/example_scripts/regularized_linear_regression/stan_batch_ess_analysis.RDS"
)
final_store_uni_ess_q <- batch_ESS_size_analysis$final_store_uni_ess_q
final_default_uni_ess_q <- batch_ESS_size_analysis$final_default_uni_ess_q

for (i in 1:(length_beta + 1)) {
  print(i)
  plot(batch_size_grid, final_store_uni_ess_q[, i], log = "y", main = paste0("i = ", i, ": automated batch size ESS = ", round(final_default_uni_ess_q[i], 1)))
}

# From the plots, the automated batch of batch size seems to be comparable with where the ESS vs batch size plot flattens out for certain coordinates, but not for all. 
# This area seems to be around 5000 and 20000 in batch size. 
# Choose batch size of 10000 for these coordinates. 

automated_batch_size_ok <- rep(T, length_beta + 1)
automated_batch_size_ok[c(1, 6, 7, 9, 12, 14)] <- FALSE # coordinates where the automated does not match well with the interval where the plot flattens out

# Second: Create summary outputs

final_stan_individual_rstan_monitor <- NULL
final_individual_result_matrix <- NULL

for (j in 1:n_iterations) {
  print(j)
  stan_individual_samples_array <- array(dim = c(n_samples_per_iteration_stan, 1, length_beta + 1))
  stan_individual_samples_array[, 1, ] <- original_scale_stan_array[, j, ]
  stan_individual_rstan_monitor <- rstan::monitor(stan_individual_samples_array, warmup = 0)
  rownames(stan_individual_rstan_monitor) <- parameter_names
  final_stan_individual_rstan_monitor <- rbind(final_stan_individual_rstan_monitor, stan_individual_rstan_monitor)
  individual_result_matrix <- 
    data.frame(
      q = colMeans(original_scale_stan_array[, j, ])
    )
  final_individual_result_matrix <- rbind(final_individual_result_matrix, individual_result_matrix)
}
final_stan_individual_rstan_monitor <- as.data.frame(final_stan_individual_rstan_monitor)
final_stan_individual_rstan_monitor

final_stan_batch_ess <- NULL

mean_q <- sapply(
  1:(length_beta + 1),
  function(i) 
    mean(final_stan_individual_rstan_monitor[grep(parameter_names[i], rownames(final_stan_individual_rstan_monitor)), ]$mean)
)
se_mean_q <- sapply(
  1:(length_beta + 1),
  function(i) 
    sd(final_stan_individual_rstan_monitor[grep(parameter_names[i], rownames(final_stan_individual_rstan_monitor)), ]$mean)
)
sd_q <- sapply(
  1:(length_beta + 1),
  function(i)
    sd(original_scale_stan_array[, , i])
)
mean_original_rstan_se_mean <- sapply(
  1:(length_beta + 1),
  function(i)
    sqrt(mean(final_stan_individual_rstan_monitor[grep(parameter_names[i], rownames(final_stan_individual_rstan_monitor)), ]$se_mean ^ 2))
)
batch_ess_q <- 
  # size chosen by the plot from above
  sapply(seq_len(ncol(stan_samples[[1]])), function(j) {
    if (!automated_batch_size_ok[j]) {
      stableGR::n.eff(lapply(stan_samples, function(x)
        x[, j, drop = FALSE]),
        multivariate = FALSE,
        size = 10000)$n.eff / length(stan_samples)  
    } else {
      stableGR::n.eff(lapply(stan_samples, function(x)
        x[, j, drop = FALSE]),
        multivariate = FALSE,
        size = NULL)$n.eff / length(stan_samples) 
    }
  })
batch_ess_q_per_mill_eval_ode <- batch_ess_q / n_evals_ode_original_scale_stan * 1e6 
final_stan_batch_ess <- 
  data.frame(
    mean_q = mean_q,
    sd_q = sd_q,
    se_mean_q = se_mean_q,
    mean_original_rstan_se_mean = mean_original_rstan_se_mean,
    batch_ess_q = batch_ess_q,
    batch_ess_q_per_mill_eval_ode = batch_ess_q_per_mill_eval_ode,
    n_evals_ode = n_evals_ode_original_scale_stan, 
    parameter_names = parameter_names
  )

saveRDS(
  final_stan_batch_ess,
  "disc_grad/example_scripts/regularized_linear_regression/stan_identity_batch_ess_results.RDS"
)

final_stan_batch_ess <- readRDS(
  "disc_grad/example_scripts/regularized_linear_regression/stan_identity_batch_ess_results.RDS"
)
