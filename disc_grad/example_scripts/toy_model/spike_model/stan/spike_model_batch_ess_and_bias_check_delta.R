rm(list = ls())
library(doParallel)
library(dplyr)
library(ggplot2)

delta <- 0.01
delta_string <- sub(".", "_", delta, fixed = T)
eta <- 1
density_kernel <- function(q) {
  exp(-(q ^ 2) / 2 + eta * pmax(0, 1 - abs(q) / delta))
}
norm_const <- integrate(density_kernel, lower = -Inf, upper = Inf, rel.tol = 1e-13, abs.tol = 1e-13)$value
density_fun <- function(q) {
  density_kernel(q) / norm_const
}
true_mean <- integrate(function(q) q * density_fun(q), lower = -Inf, upper = Inf, rel.tol = 1e-13, abs.tol = 1e-13)$value
true_mean
true_var <- integrate(function(q) (q ^  2) * density_fun(q), lower = -Inf, upper = Inf, rel.tol = 1e-13, abs.tol = 1e-13)$value
true_var
true_prob_abs_q_below_delta <- integrate(function(q) density_fun(q), lower = -delta, upper = delta, rel.tol = 1e-13, abs.tol = 1e-13)$value
true_prob_abs_q_below_delta

dim_vec <- 2 ^ {1:10}
n_chains <- 10

n_samples <- 10000

# First: Preliminary check if the automated chosen batch size (i.e. size = NULL) in stableGR is adequate.
# If the output batch ESS is similar to the value where the batch ESS flattens out with respect to the batch size, then the default is used.
# Otherwise, the output batch ESS is too optimistic and a manual batch size value is chosen close to the elbow of the ESS vs batch size plot or the first part where the plot is flat.
final_store_uni_ess_q <- NULL
final_store_uni_ess_q_squared <- NULL
final_store_uni_ess_abs_q_below_delta <- NULL

final_default_uni_ess_q <- NULL
final_default_uni_ess_q_squared <- NULL
final_default_uni_ess_abs_q_below_delta <- NULL

for (i in 1:length(dim_vec)) {
  print(i)
  n_dim <- dim_vec[i]
  stan_run <- readRDS(paste0("disc_grad/example_scripts/toy_model/spike_model/stan/delta_0_01/spike_model_delta_0_01_dim_", n_dim, ".RDS"))
  stan_samples_array <- rstan::extract(stan_run, permute = F)
  stan_samples <-
    lapply(1:n_chains, function(i)
      (stan_samples_array[, i, 1:2]))
  stan_samples_squared <-
    lapply(1:n_chains, function(i)
      (stan_samples_array[, i, 1:2] ^ 2))
  stan_samples_abs_q_below_delta <- lapply(1:n_chains,
                                           function(i)
                                             (abs(stan_samples_array[, i, 1:2]) < delta) * 1L)
  n_splits <- 1000 # want to check for 1000 different batch sizes
  n_cores <- 10
  batch_size_grid <- seq(from = n_samples / n_splits, to = n_samples, by = n_samples / n_splits)
  log_path <- "disc_grad/example_scripts/toy_model/spike_model/stan/batch_ess_progress.log"
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
  final_store_uni_ess_q <- append(final_store_uni_ess_q, list(store_uni_ess_q))
  final_default_uni_ess_q <- append(final_default_uni_ess_q, list(default_uni_ess_q))
  
  store_uni_ess_q_squared <- foreach::foreach(i = 1:n_splits, .combine = rbind) %dopar% {
    print(i)
    sapply(seq_len(ncol(stan_samples_squared[[1]])), function(j) {
      stableGR::n.eff(
        lapply(stan_samples_squared, function(x)
          x[, j, drop = FALSE]),
        multivariate = FALSE,
        size = batch_size_grid[i]
      )$n.eff / length(stan_samples_squared)
    })
  }
  default_uni_ess_q_squared <-
    sapply(seq_len(ncol(stan_samples_squared[[1]])), function(j) {
      stableGR::n.eff(
        lapply(stan_samples_squared, function(x)
          x[, j, drop = FALSE]),
        multivariate = FALSE,
        size = NULL
      )$n.eff / length(stan_samples_squared)
    })
  final_store_uni_ess_q_squared <- append(final_store_uni_ess_q_squared, list(store_uni_ess_q_squared))
  final_default_uni_ess_q_squared <- append(final_default_uni_ess_q_squared, list(default_uni_ess_q_squared))
  
  if (i > 1) {
    store_uni_ess_abs_q_below_delta <- foreach::foreach(i = 1:n_splits, .combine = rbind) %dopar% {
      print(i)
      sapply(seq_len(ncol(stan_samples_abs_q_below_delta[[1]])), function(j) {
        stableGR::n.eff(
          lapply(stan_samples_abs_q_below_delta, function(x)
            x[, j, drop = FALSE]),
          multivariate = FALSE,
          size = batch_size_grid[i]
        )$n.eff / length(stan_samples_abs_q_below_delta)
      })
    }
    default_uni_ess_abs_q_below_delta <-
      sapply(seq_len(ncol(stan_samples_abs_q_below_delta[[1]])), function(j) {
        stableGR::n.eff(
          lapply(stan_samples_abs_q_below_delta, function(x)
            x[, j, drop = FALSE]),
          multivariate = FALSE,
          size = NULL
        )$n.eff / length(stan_samples_abs_q_below_delta)
      })
    final_store_uni_ess_abs_q_below_delta <- append(final_store_uni_ess_abs_q_below_delta, list(store_uni_ess_abs_q_below_delta))
    final_default_uni_ess_abs_q_below_delta <- append(final_default_uni_ess_abs_q_below_delta, list(default_uni_ess_abs_q_below_delta)) 
  }
  
  parallel::stopCluster(init_cluster)
}

saveRDS(
  list(
    final_store_uni_ess_q = final_store_uni_ess_q,
    final_store_uni_ess_q_squared = final_store_uni_ess_q_squared,
    final_store_uni_ess_abs_q_below_delta = final_store_uni_ess_abs_q_below_delta,
    final_default_uni_ess_q = final_default_uni_ess_q,
    final_default_uni_ess_q_squared = final_default_uni_ess_q_squared,
    final_default_uni_ess_abs_q_below_delta = final_default_uni_ess_abs_q_below_delta
  ),
  "disc_grad/example_scripts/toy_model/spike_model/stan/batch_ESS_size_analysis_delta_0_01.RDS"
)

batch_ESS_size_analysis_delta_0_01 <- readRDS(
  "disc_grad/example_scripts/toy_model/spike_model/stan/batch_ESS_size_analysis_delta_0_01.RDS"
)
final_store_uni_ess_q <- batch_ESS_size_analysis_delta_0_01$final_store_uni_ess_q
final_store_uni_ess_q_squared <- batch_ESS_size_analysis_delta_0_01$final_store_uni_ess_q_squared
final_store_uni_ess_abs_q_below_delta <- batch_ESS_size_analysis_delta_0_01$final_store_uni_ess_abs_q_below_delta
final_default_uni_ess_q <- batch_ESS_size_analysis_delta_0_01$final_default_uni_ess_q
final_default_uni_ess_q_squared <- batch_ESS_size_analysis_delta_0_01$final_default_uni_ess_q_squared
final_default_uni_ess_abs_q_below_delta <- batch_ESS_size_analysis_delta_0_01$final_default_uni_ess_abs_q_below_delta

for (i in 1:length(dim_vec)) {
  print(i)
  # print("ess_q")
  plot(batch_size_grid, final_store_uni_ess_q[[i]][, 1], log = "y", main = paste0("i = ", i, ": automated batch size ESS = ", round(final_default_uni_ess_q[[i]][1], 1)))
  # print(final_default_uni_ess_q[[i]])
  # print("ess_q_squared")
  plot(batch_size_grid, final_store_uni_ess_q_squared[[i]][, 1], log = "y", main = paste0("i = ", i, ": automated batch size ESS = ", round(final_default_uni_ess_q_squared[[i]][1], 1)))
  # print(final_default_uni_ess_q_squared[[i]])
  if (i > 1) {
    # print("ess_abs_q_below_delta")
    plot(batch_size_grid, final_store_uni_ess_abs_q_below_delta[[i - 1]][, 1], log = "y", main = paste0("i = ", i, ": automated batch size ESS = ", round(final_default_uni_ess_abs_q_below_delta[[i - 1]][1])))
    # print(final_default_uni_ess_abs_q_below_delta[[i - 1]])  
  }
}

# The batch ESS given by the automated choice of batch size is roughly the same as for the batch size values where the ESS flattens out first.
# In this case, no need to manually choose the batch size. 

# Second: Create summary outputs

full_final_stan_batch_ess <- NULL
full_final_result_matrix <- NULL

for (i in 1:length(dim_vec)) {
  print(i)
  n_dim <- dim_vec[i]
  
  stan_run <- readRDS(paste0("disc_grad/example_scripts/toy_model/spike_model/stan/delta_0_01/spike_model_delta_0_01_dim_", n_dim, ".RDS"))
  stan_samples <- rstan::extract(stan_run, permute = FALSE)
  stan_n_evals_ode <- sum(rstan::get_num_leapfrog_per_iteration(stan_run))
  
  # stan_n_boundary_crossings <- sum(sapply(1:n_chains, function(j) stan_run[[j]]$n_boundary_crossings))
  
  final_stan_individual_rstan_monitor <- NULL
  final_individual_result_matrix <- NULL
  
  parameter_names <- c("q1", "q2")
  
  for (j in 1:n_chains) {
    stan_individual_samples_array <- array(dim = c(n_samples, 1, 2))
    stan_individual_samples_array[, 1, ] <- stan_samples[, j, 1:2]
    stan_individual_rstan_monitor <- rstan::monitor(stan_individual_samples_array, warmup = 0)
    rownames(stan_individual_rstan_monitor) <- parameter_names
    final_stan_individual_rstan_monitor <- rbind(final_stan_individual_rstan_monitor, stan_individual_rstan_monitor)
    individual_result_matrix <- 
      data.frame(
        q1 = mean(stan_samples[, j, 1]),
        q1_squared = mean(stan_samples[, j, 1] ^ 2),
        abs_q1_below_delta = mean(abs(stan_samples[, j, 1]) < delta),
        n_dim = n_dim
      )
    final_individual_result_matrix <- rbind(final_individual_result_matrix, individual_result_matrix)
  }
  final_stan_individual_rstan_monitor <- as.data.frame(final_stan_individual_rstan_monitor)
  final_stan_individual_rstan_monitor
  
  final_stan_batch_ess <- NULL
  
  for (p in 1:(2)) {
    q_pooled_samples <- sapply(1:n_chains, function(k) stan_samples[, k, p])
    q_list_chains <- lapply(1:n_chains, function(k) matrix(stan_samples[, k, p], ncol = 1))
    q_squared_list_chains <- lapply(1:n_chains, function(k) matrix(stan_samples[, k, p] ^ 2, ncol = 1))
    abs_q_below_delta_list_chains <- lapply(1:n_chains, function(k) matrix(as.integer(abs(stan_samples[, k, p]) < delta), ncol = 1))
    mean_q <- mean(final_stan_individual_rstan_monitor[grep(parameter_names[p], rownames(final_stan_individual_rstan_monitor)), ]$mean)
    mean_q_squared <- mean(q_pooled_samples ^ 2)
    mean_abs_q_below_delta <- mean(abs(q_pooled_samples) < delta)
    se_mean_q <- sd(final_stan_individual_rstan_monitor[grep(parameter_names[p], rownames(final_stan_individual_rstan_monitor)), ]$mean)
    se_mean_q_squared <- sd(sapply(1:n_chains, function(k) mean(stan_samples[, k, p] ^ 2)))
    se_mean_abs_q_below_delta <- sd(sapply(1:n_chains, function(k) mean(abs(stan_samples[, k, p]) < delta)))
    sd_q <- sd(q_pooled_samples)
    sd_q_squared <- sd(q_pooled_samples ^ 2)
    sd_abs_q_below_delta <- sd(abs(q_pooled_samples) < delta)
    mean_original_rstan_se_mean <- 
      sqrt(mean(final_stan_individual_rstan_monitor[grep(parameter_names[p], rownames(final_stan_individual_rstan_monitor)), ]$se_mean ^ 2))
    # batch_ess_q <- ((sd_q / se_mean_q) ^ 2) * n_chains
    batch_ess_q <- stableGR::n.eff(q_list_chains, multivariate = F, size = NULL)$n.eff / n_chains # must divide by the number of chains manually here due to a bug in the package
    batch_ess_q_per_mill_eval_ode <- batch_ess_q / stan_n_evals_ode * 1e6 
    # batch_ess_q_squared <- ((sd_q_squared / se_mean_q_squared) ^ 2) * n_chains
    batch_ess_q_squared <- stableGR::n.eff(q_squared_list_chains, multivariate = F, size = NULL)$n.eff / n_chains
    batch_ess_q_squared_per_mill_eval_ode <- batch_ess_q_squared / stan_n_evals_ode * 1e6 
    # batch_ess_abs_q_below_delta <- ((sd_abs_q_below_delta / se_mean_abs_q_below_delta) ^ 2) * n_chains
    if (i > 1) {
      batch_ess_abs_q_below_delta <- stableGR::n.eff(abs_q_below_delta_list_chains, multivariate = F, size = NULL)$n.eff / n_chains
      batch_ess_abs_q_below_delta_per_mill_eval_ode <- batch_ess_abs_q_below_delta / stan_n_evals_ode * 1e6
    } else {
      batch_ess_abs_q_below_delta <- NA
      batch_ess_abs_q_below_delta_per_mill_eval_ode <- NA
    }
    final_stan_batch_ess <- rbind(
      final_stan_batch_ess,
      data.frame(
        mean_q = mean_q,
        mean_q_squared = mean_q_squared,
        mean_abs_q_below_delta = mean_abs_q_below_delta,
        sd_q = sd_q,
        sd_q_squared = sd_q_squared,
        sd_abs_q_below_delta = sd_abs_q_below_delta,
        se_mean_q = se_mean_q,
        se_mean_q_squared = se_mean_q_squared,
        se_mean_abs_q_below_delta = se_mean_abs_q_below_delta,
        mean_original_rstan_se_mean = mean_original_rstan_se_mean,
        batch_ess_q = batch_ess_q,
        batch_ess_q_per_mill_eval_ode = batch_ess_q_per_mill_eval_ode,
        batch_ess_q_squared = batch_ess_q_squared,
        batch_ess_q_squared_per_mill_eval_ode = batch_ess_q_squared_per_mill_eval_ode,
        batch_ess_abs_q_below_delta = batch_ess_abs_q_below_delta,
        batch_ess_abs_q_below_delta_per_mill_eval_ode = batch_ess_abs_q_below_delta_per_mill_eval_ode,
        # n_boundary_crossings = stan_n_boundary_crossings,
        n_evals_ode = stan_n_evals_ode, 
        n_dim = n_dim,
        parameter_names = parameter_names[p]
      )
    )
  }
  rownames(final_stan_batch_ess) <- parameter_names
  full_final_stan_batch_ess <- rbind(full_final_stan_batch_ess, final_stan_batch_ess)
  full_final_result_matrix <- rbind(full_final_result_matrix, final_individual_result_matrix)
  
}

full_final_diff_result_matrix <- 
  sweep(full_final_result_matrix, 2, c(true_mean, true_var, true_prob_abs_q_below_delta, 0), "-")

saveRDS(
  list(
    full_final_stan_batch_ess = full_final_stan_batch_ess,
    full_final_result_matrix = full_final_result_matrix,
    full_final_diff_result_matrix = full_final_diff_result_matrix
  ),
  paste0("disc_grad/example_scripts/toy_model/spike_model/stan/spike_model_results_delta_", delta_string, ".RDS")
)

(full_final_stan_batch_ess[full_final_stan_batch_ess$parameter_names == "q1" & full_final_stan_batch_ess$n_dim == 4, ]$mean_q - true_mean)
(full_final_stan_batch_ess[full_final_stan_batch_ess$parameter_names == "q1" & full_final_stan_batch_ess$n_dim == 4, ]$mean_q - true_mean) - 
  1.96 * full_final_stan_batch_ess[full_final_stan_batch_ess$parameter_names == "q1" & full_final_stan_batch_ess$n_dim == 4, ]$se_mean_q / sqrt(n_chains)
full_final_stan_batch_ess[full_final_stan_batch_ess$parameter_names == "q1" & full_final_stan_batch_ess$n_dim == 4, ]$mean_q + 
  1.96 * full_final_stan_batch_ess[full_final_stan_batch_ess$parameter_names == "q1" & full_final_stan_batch_ess$n_dim == 4, ]$se_mean_q / sqrt(n_chains)

(full_final_stan_batch_ess[full_final_stan_batch_ess$parameter_names == "q1" & full_final_stan_batch_ess$n_dim == 4, ]$mean_q_squared - true_var)
(full_final_stan_batch_ess[full_final_stan_batch_ess$parameter_names == "q1" & full_final_stan_batch_ess$n_dim == 4, ]$mean_q_squared - true_var) - 
  1.96 * full_final_stan_batch_ess[full_final_stan_batch_ess$parameter_names == "q1" & full_final_stan_batch_ess$n_dim == 4, ]$se_mean_q_squared / sqrt(n_chains)
(full_final_stan_batch_ess[full_final_stan_batch_ess$parameter_names == "q1" & full_final_stan_batch_ess$n_dim == 4, ]$mean_q_squared - true_var) + 
  1.96 * full_final_stan_batch_ess[full_final_stan_batch_ess$parameter_names == "q1" & full_final_stan_batch_ess$n_dim == 4, ]$se_mean_q_squared / sqrt(n_chains)

(full_final_stan_batch_ess[full_final_stan_batch_ess$parameter_names == "q1" & full_final_stan_batch_ess$n_dim == 4, ]$mean_abs_q_below_delta - true_prob_abs_q_below_delta)
(full_final_stan_batch_ess[full_final_stan_batch_ess$parameter_names == "q1" & full_final_stan_batch_ess$n_dim == 4, ]$mean_abs_q_below_delta - true_prob_abs_q_below_delta) - 
  1.96 * full_final_stan_batch_ess[full_final_stan_batch_ess$parameter_names == "q1" & full_final_stan_batch_ess$n_dim == 4, ]$se_mean_abs_q_below_delta / sqrt(n_chains)
(full_final_stan_batch_ess[full_final_stan_batch_ess$parameter_names == "q1" & full_final_stan_batch_ess$n_dim == 4, ]$mean_abs_q_below_delta - true_prob_abs_q_below_delta) + 
  1.96 * full_final_stan_batch_ess[full_final_stan_batch_ess$parameter_names == "q1" & full_final_stan_batch_ess$n_dim == 4, ]$se_mean_abs_q_below_delta / sqrt(n_chains)

