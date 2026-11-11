rm(list = ls())
library(doParallel)
library(dplyr)
library(ggplot2)

set.seed(42)

n_cores <- 5
n_iterations <- 10
n_samples_per_iteration <- 100000
sampling_time_period_per_iteration <- 10000

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

# GRHMC

grhmc_disc_grad_run <- readRDS("disc_grad/example_scripts/bayesian_neural_network/bnn_disc_grad_run_VARI_with_n_evals_ode_run_example.RDS")

grhmc_disc_grad_array <- array(dim = c(n_samples_per_iteration, n_iterations, ncol(grhmc_disc_grad_run[[1]])))

for (i in 1:n_iterations) {
  
  grhmc_disc_grad_array[, i, ] <- grhmc_disc_grad_run[[i]]
  
}

n_evals_ode_grhmc_disc_grad <- 0
for (i in 1:n_iterations) {
  n_evals_ode_grhmc_disc_grad <- n_evals_ode_grhmc_disc_grad + 
    grhmc_disc_grad_run[, 2][[i]]
}

grhmc_samples <- # list of samples of q directly
  lapply(1:n_iterations, function(i)
    (grhmc_disc_grad_array[, i, ])) 

n_splits <- 1000 # want to check for 1000 different batch sizes
n_cores <- 10
batch_size_grid <- seq(from = n_samples_per_iteration / n_splits, to = n_samples_per_iteration, by = n_samples_per_iteration / n_splits)
log_path <- "disc_grad/example_scripts/bayesian_neural_network/batch_ess_progress.log"
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
  "disc_grad/example_scripts/bayesian_neural_network/disc_grad_run_VARI_batch_ess_analysis.RDS"
)

batch_ESS_size_analysis <- readRDS(
  "disc_grad/example_scripts/bayesian_neural_network/disc_grad_run_VARI_batch_ess_analysis.RDS"
)
final_store_uni_ess_q <- batch_ESS_size_analysis$final_store_uni_ess_q
final_default_uni_ess_q <- batch_ESS_size_analysis$final_default_uni_ess_q

for (i in 1:ncol(grhmc_disc_grad_run[[1]])) {
  print(i)
  plot(batch_size_grid, final_store_uni_ess_q[, i], log = "y", main = paste0("i = ", i, ": automated batch size ESS = ", round(final_default_uni_ess_q[i], 1)))
}

# From the first two plots (as the weights themselves are not relevant due to identifiability issues),
# the batch ESS given by the automated choice of batch size is roughly the same as for the batch size values where the ESS flattens out first.
# In this case, no need to manually choose the batch size. 

# Second: Create summary outputs
parameter_names <- c(
  "gamma",     
  "alpha",     
  "w_star[1]", 
  "w_star[2]", 
  "delta[1]",  
  "delta[2]",  
  "beta1[1]",  
  "beta1[2]", 
  "beta2[1]",  
  "beta2[2]",  
  "dummy_parameter"
)
final_grhmc_individual_rstan_monitor <- NULL
final_individual_result_matrix <- NULL

for (j in 1:n_iterations) {
  print(j)
  grhmc_individual_samples_array <- array(dim = c(n_samples_per_iteration, 1, ncol(grhmc_disc_grad_run[[1]])))
  grhmc_individual_samples_array[, 1, ] <- grhmc_disc_grad_array[, j, ]
  grhmc_individual_rstan_monitor <- rstan::monitor(grhmc_individual_samples_array, warmup = 0)
  rownames(grhmc_individual_rstan_monitor) <- parameter_names
  final_grhmc_individual_rstan_monitor <- rbind(final_grhmc_individual_rstan_monitor, grhmc_individual_rstan_monitor)
  individual_result_matrix <- 
    data.frame(
      q = colMeans(grhmc_disc_grad_array[, j, ])
    )
  final_individual_result_matrix <- rbind(final_individual_result_matrix, individual_result_matrix)
}
final_grhmc_individual_rstan_monitor <- as.data.frame(final_grhmc_individual_rstan_monitor)
final_grhmc_individual_rstan_monitor

final_grhmc_batch_ess <- NULL

mean_q <- sapply(
  1:ncol(grhmc_disc_grad_run[[1]]),
  function(i) 
    mean(final_grhmc_individual_rstan_monitor[grep(parameter_names[i], rownames(final_grhmc_individual_rstan_monitor), fixed = T), ]$mean)
)
se_mean_q <- sapply(
  1:ncol(grhmc_disc_grad_run[[1]]),
  function(i) 
    sd(final_grhmc_individual_rstan_monitor[grep(parameter_names[i], rownames(final_grhmc_individual_rstan_monitor), fixed = T), ]$mean)
)
sd_q <- sapply(
  1:ncol(grhmc_disc_grad_run[[1]]),
  function(i)
    sd(grhmc_disc_grad_array[, , i])
)
mean_original_rstan_se_mean <- sapply(
  1:ncol(grhmc_disc_grad_run[[1]]),
  function(i)
    sqrt(mean(final_grhmc_individual_rstan_monitor[grep(parameter_names[i], rownames(final_grhmc_individual_rstan_monitor), fixed = T), ]$se_mean ^ 2))
)
batch_ess_q <- 
  # size chosen by the plot from above
  sapply(seq_len(ncol(grhmc_samples[[1]])), function(j) {
    stableGR::n.eff(lapply(grhmc_samples, function(x)
      x[, j, drop = FALSE]),
      multivariate = FALSE,
      size = NULL)$n.eff / length(grhmc_samples) 
  })
batch_ess_q_per_mill_eval_ode <- batch_ess_q / n_evals_ode_grhmc_disc_grad * 1e6 
final_grhmc_batch_ess <- 
  data.frame(
    mean_q = mean_q,
    sd_q = sd_q,
    se_mean_q = se_mean_q,
    mean_original_rstan_se_mean = mean_original_rstan_se_mean,
    batch_ess_q = batch_ess_q,
    batch_ess_q_per_mill_eval_ode = batch_ess_q_per_mill_eval_ode,
    n_evals_ode = n_evals_ode_grhmc_disc_grad, 
    parameter_names = parameter_names
  )

saveRDS(
  final_grhmc_batch_ess,
  "disc_grad/example_scripts/bayesian_neural_network/disc_grad_run_VARI_batch_ess_results.RDS"
)

final_grhmc_batch_ess <- readRDS(
  "disc_grad/example_scripts/bayesian_neural_network/disc_grad_run_VARI_batch_ess_results.RDS"
)

# Stan

stan_run <- readRDS("disc_grad/example_scripts/bayesian_neural_network/stan_run_example.RDS")
stan_samples_array <- rstan::extract(stan_run, permute = FALSE)
n_evals_ode_stan <- sum(rstan::get_num_leapfrog_per_iteration(stan_run))
stan_samples <- # list of samples of q directly
  lapply(1:n_iterations, function(i)
    (stan_samples_array[, i, ])) 

n_splits <- 1000 # want to check for 1000 different batch sizes
n_cores <- 10
batch_size_grid <- seq(from = n_samples_per_iteration / n_splits, to = n_samples_per_iteration, by = n_samples_per_iteration / n_splits)
log_path <- "disc_grad/example_scripts/bayesian_neural_network/batch_ess_progress.log"
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
  "disc_grad/example_scripts/bayesian_neural_network/stan_batch_ess_analysis.RDS"
)

batch_ESS_size_analysis <- readRDS(
  "disc_grad/example_scripts/bayesian_neural_network/stan_batch_ess_analysis.RDS"
)
final_store_uni_ess_q <- batch_ESS_size_analysis$final_store_uni_ess_q
final_default_uni_ess_q <- batch_ESS_size_analysis$final_default_uni_ess_q

for (i in 1:ncol(grhmc_disc_grad_run[[1]])) {
  print(i)
  plot(batch_size_grid, final_store_uni_ess_q[, i], log = "y", main = paste0("i = ", i, ": automated batch size ESS = ", round(final_default_uni_ess_q[i], 1)))
}

# From the plots, the automated batch of batch size seems to be too optimistic compared to where it flattens out. 
# The plots flatten out when batch size is around 10000. This is therefore chosen as the batch size in this case. 

# Second: Create summary outputs

final_stan_individual_rstan_monitor <- NULL
final_individual_result_matrix <- NULL

for (j in 1:n_iterations) {
  print(j)
  stan_individual_samples_array <- array(dim = c(n_samples_per_iteration, 1, ncol(stan_samples[[1]])))
  stan_individual_samples_array[, 1, ] <- stan_samples_array[, j, ]
  stan_individual_rstan_monitor <- rstan::monitor(stan_individual_samples_array, warmup = 0)
  rownames(stan_individual_rstan_monitor)[1:length(parameter_names)] <- parameter_names
  final_stan_individual_rstan_monitor <- rbind(final_stan_individual_rstan_monitor, stan_individual_rstan_monitor)
  individual_result_matrix <- 
    data.frame(
      q = colMeans(stan_samples_array[, j, ])
    )
  final_individual_result_matrix <- rbind(final_individual_result_matrix, individual_result_matrix)
}
final_stan_individual_rstan_monitor <- as.data.frame(final_stan_individual_rstan_monitor)
final_stan_individual_rstan_monitor

final_stan_batch_ess <- NULL

mean_q <- sapply(
  1:(ncol(stan_samples[[1]]) - 1),
  function(i) 
    mean(final_stan_individual_rstan_monitor[grep(parameter_names[i], rownames(final_stan_individual_rstan_monitor), fixed = T), ]$mean)
)
se_mean_q <- sapply(
  1:(ncol(stan_samples[[1]]) - 1),
  function(i) 
    sd(final_stan_individual_rstan_monitor[grep(parameter_names[i], rownames(final_stan_individual_rstan_monitor), fixed = T), ]$mean)
)
sd_q <- sapply(
  1:(ncol(stan_samples[[1]]) - 1),
  function(i)
    sd(stan_samples_array[, , i])
)
mean_original_rstan_se_mean <- sapply(
  1:(ncol(stan_samples[[1]]) - 1),
  function(i)
    sqrt(mean(final_stan_individual_rstan_monitor[grep(parameter_names[i], rownames(final_stan_individual_rstan_monitor), fixed = T), ]$se_mean ^ 2))
)
batch_ess_q <- 
  sapply(seq_len(ncol(stan_samples[[1]]) - 1), function(j) {
    stableGR::n.eff(lapply(stan_samples, function(x)
      x[, j, drop = FALSE]),
      multivariate = FALSE,
      size = 10000)$n.eff / length(stan_samples) 
  })
batch_ess_q_per_mill_eval_ode <- batch_ess_q / n_evals_ode_stan * 1e6 
final_stan_batch_ess <- 
  data.frame(
    mean_q = mean_q,
    sd_q = sd_q,
    se_mean_q = se_mean_q,
    mean_original_rstan_se_mean = mean_original_rstan_se_mean,
    batch_ess_q = batch_ess_q,
    batch_ess_q_per_mill_eval_ode = batch_ess_q_per_mill_eval_ode,
    n_evals_ode = n_evals_ode_stan, 
    parameter_names = parameter_names
  )

saveRDS(
  final_stan_batch_ess,
  "disc_grad/example_scripts/bayesian_neural_network/stan_var_batch_ess_results.RDS"
)

final_stan_batch_ess <- readRDS(
  "disc_grad/example_scripts/bayesian_neural_network/stan_var_batch_ess_results.RDS"
)
