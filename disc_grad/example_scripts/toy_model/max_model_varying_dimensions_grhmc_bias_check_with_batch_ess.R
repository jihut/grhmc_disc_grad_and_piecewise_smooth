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

estimated_quantities <- function(list, calculate_wasserstein_2 = FALSE) {
  
  samples_matrix <- list$q_original_samples
  
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

# First: Preliminary check if the automated chosen batch size (i.e. size = NULL) in stableGR is adequate.
# If the output batch ESS is similar to the value where the batch ESS flattens out with respect to the batch size, then the default is used.
# Otherwise, the output batch ESS is too optimistic and a manual batch size value is chosen close to the elbow of the ESS vs batch size plot or the first part where the plot is flat.
final_store_uni_ess_q <- NULL
final_store_uni_ess_q_squared <- NULL
final_store_uni_ess_q1_squared_plus_q2_squared <- NULL

final_default_uni_ess_q <- NULL
final_default_uni_ess_q_squared <- NULL
final_default_uni_ess_q1_squared_plus_q2_squared<- NULL

for (i in 1:length(dim_vec)) {
  print(i)
  n_dim <- dim_vec[i]
  grhmc_run <- readRDS(
    paste0("disc_grad/example_scripts/toy_model/grhmc/max_model_d_", n_dim, "_grhmc.RDS")
  )
  grhmc_samples <- # list of samples of q directly
    lapply(1:n_chains, function(i)
      matrix(grhmc_run[[i]]$q_original_samples[, 2], ncol = 1))
  grhmc_samples_squared <- # list of samples of q^2
    lapply(1:n_chains, function(i)
      matrix(grhmc_run[[i]]$q_original_samples[, 2] ^ 2, ncol = 1))
  grhmc_samples_q1_squared_plus_q2_squared <- lapply(1:n_chains,
                                                     function(i)
                                                       matrix(rowSums(grhmc_run[[i]]$q_original_samples[, 1:2] ^ 2), ncol = 1))
  n_splits <- 1000 # want to check for 1000 different batch sizes
  n_cores <- 10
  batch_size_grid <- seq(from = n_samples / n_splits, to = n_samples, by = n_samples / n_splits)
  log_path <- "disc_grad/example_scripts/toy_model/batch_ess_progress.log"
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
  final_store_uni_ess_q <- append(final_store_uni_ess_q, list(store_uni_ess_q))
  final_default_uni_ess_q <- append(final_default_uni_ess_q, list(default_uni_ess_q))
  
  store_uni_ess_q_squared <- foreach::foreach(i = 1:n_splits, .combine = rbind) %dopar% {
    print(i)
    sapply(seq_len(ncol(grhmc_samples_squared[[1]])), function(j) {
      stableGR::n.eff(
        lapply(grhmc_samples_squared, function(x)
          x[, j, drop = FALSE]),
        multivariate = FALSE,
        size = batch_size_grid[i]
      )$n.eff / length(grhmc_samples_squared)
    })
  }
  default_uni_ess_q_squared <-
    sapply(seq_len(ncol(grhmc_samples_squared[[1]])), function(j) {
      stableGR::n.eff(
        lapply(grhmc_samples_squared, function(x)
          x[, j, drop = FALSE]),
        multivariate = FALSE,
        size = NULL
      )$n.eff / length(grhmc_samples_squared)
    })
  final_store_uni_ess_q_squared <- append(final_store_uni_ess_q_squared, list(store_uni_ess_q_squared))
  final_default_uni_ess_q_squared <- append(final_default_uni_ess_q_squared, list(default_uni_ess_q_squared))
  
  store_uni_ess_q1_squared_plus_q2_squared <- foreach::foreach(i = 1:n_splits, .combine = rbind) %dopar% {
    print(i)
    sapply(seq_len(ncol(grhmc_samples_q1_squared_plus_q2_squared[[1]])), function(j) {
      stableGR::n.eff(
        lapply(grhmc_samples_q1_squared_plus_q2_squared, function(x)
          x[, j, drop = FALSE]),
        multivariate = FALSE,
        size = batch_size_grid[i]
      )$n.eff / length(grhmc_samples_q1_squared_plus_q2_squared)
    })
  }
  default_uni_ess_q1_squared_plus_q2_squared <-
    sapply(seq_len(ncol(grhmc_samples_q1_squared_plus_q2_squared[[1]])), function(j) {
      stableGR::n.eff(
        lapply(grhmc_samples_q1_squared_plus_q2_squared, function(x)
          x[, j, drop = FALSE]),
        multivariate = FALSE,
        size = NULL
      )$n.eff / length(grhmc_samples_q1_squared_plus_q2_squared)
    })
  final_store_uni_ess_q1_squared_plus_q2_squared <- append(final_store_uni_ess_q1_squared_plus_q2_squared, list(store_uni_ess_q1_squared_plus_q2_squared))
  final_default_uni_ess_q1_squared_plus_q2_squared <- append(final_default_uni_ess_q1_squared_plus_q2_squared, list(default_uni_ess_q1_squared_plus_q2_squared))
  
  parallel::stopCluster(init_cluster)
}

saveRDS(
  list(
    final_store_uni_ess_q = final_store_uni_ess_q,
    final_store_uni_ess_q_squared = final_store_uni_ess_q_squared,
    final_store_uni_ess_q1_squared_plus_q2_squared = final_store_uni_ess_q1_squared_plus_q2_squared,
    final_default_uni_ess_q = final_default_uni_ess_q,
    final_default_uni_ess_q_squared = final_default_uni_ess_q_squared,
    final_default_uni_ess_q1_squared_plus_q2_squared = final_default_uni_ess_q1_squared_plus_q2_squared
  ),
  "disc_grad/example_scripts/toy_model/max_model_variying_dimensions_grhmc_bias_check_with_batch_ess_analysis.RDS"
)

batch_ESS_size_analysis <- readRDS(
  "disc_grad/example_scripts/toy_model/max_model_variying_dimensions_grhmc_bias_check_with_batch_ess_analysis.RDS"
)
final_store_uni_ess_q <- batch_ESS_size_analysis$final_store_uni_ess_q
final_store_uni_ess_q_squared <- batch_ESS_size_analysis$final_store_uni_ess_q_squared
final_store_uni_ess_q1_squared_plus_q2_squared <- batch_ESS_size_analysis$final_store_uni_ess_q1_squared_plus_q2_squared
final_default_uni_ess_q <- batch_ESS_size_analysis$final_default_uni_ess_q
final_default_uni_ess_q_squared <- batch_ESS_size_analysis$final_default_uni_ess_q_squared
final_default_uni_ess_q1_squared_plus_q2_squared <- batch_ESS_size_analysis$final_default_uni_ess_q1_squared_plus_q2_squared

for (i in 1:length(dim_vec)) {
  print(i)
  # print("ess_q")
  plot(batch_size_grid, final_store_uni_ess_q[[i]][, 1], log = "y", main = paste0("i = ", i, ": automated batch size ESS = ", round(final_default_uni_ess_q[[i]][1], 1)))
  # print(final_default_uni_ess_q[[i]])
  # print("ess_q_squared")
  plot(batch_size_grid, final_store_uni_ess_q_squared[[i]][, 1], log = "y", main = paste0("i = ", i, ": automated batch size ESS = ", round(final_default_uni_ess_q_squared[[i]][1], 1)))
  # print(final_default_uni_ess_q_squared[[i]])
  # print("ess_q1_squared_plus_q2_squared")
  plot(batch_size_grid, final_store_uni_ess_q1_squared_plus_q2_squared[[i]][, 1], log = "y", main = paste0("i = ", i, ": automated batch size ESS = ", round(final_default_uni_ess_q1_squared_plus_q2_squared[[i]][1])))
  # print(final_default_uni_ess_q1_squared_plus_q2_squared[[i - 1]])  
}

# The batch ESS given by the automated choice of batch size seems to be reasonable for 2 ^ {1:6}.
# At 2^{7} and larger, the automated choice gives too optimistic values based on the plots for all the quantities of interest. 
# Based on the plots, they seem to be flatten out if one set the batch size around 1000. Therefore, this is the chosen value for all quantities when i > 6. 

# Second: Create summary outputs

full_final_grhmc_batch_ess <- NULL
full_final_result_matrix <- NULL

for (i in 1:length(dim_vec)) {
  print(i)
  n_dim <- dim_vec[i]
  
  grhmc_run <- readRDS(
    paste0("disc_grad/example_scripts/toy_model/grhmc/max_model_d_", n_dim, "_grhmc.RDS")
  )
  
  grhmc_n_evals_ode <- sum(sapply(1:n_chains, function(j) grhmc_run[[j]]$n_evals_ode))
  
  # grhmc_n_boundary_crossings <- sum(sapply(1:n_chains, function(j) grhmc_run[[j]]$n_boundary_crossings))
  
  final_grhmc_individual_rstan_monitor <- NULL
  final_individual_result_matrix <- NULL
  
  parameter_names <- c("q1", "q2")
  
  for (j in 1:n_chains) {
    grhmc_individual_samples_array <- array(dim = c(n_samples, 1, 2))
    grhmc_individual_samples_array[, 1, ] <- grhmc_run[[j]]$q_original_samples[, 1:2]
    grhmc_individual_rstan_monitor <- rstan::monitor(grhmc_individual_samples_array, warmup = 0)
    rownames(grhmc_individual_rstan_monitor) <- parameter_names
    final_grhmc_individual_rstan_monitor <- rbind(final_grhmc_individual_rstan_monitor, grhmc_individual_rstan_monitor)
    individual_result_matrix <- 
      data.frame(
        q2 = mean(grhmc_run[[j]]$q_original_samples[, 2]),
        q2_squared = mean(grhmc_run[[j]]$q_original_samples[, 2] ^ 2),
        q1_squared_plus_q2_squared = mean(rowSums(grhmc_run[[j]]$q_original_samples[, 1:2] ^ 2)),
        n_dim = n_dim
      )
    final_individual_result_matrix <- rbind(final_individual_result_matrix, individual_result_matrix)
  }
  final_grhmc_individual_rstan_monitor <- as.data.frame(final_grhmc_individual_rstan_monitor)
  final_grhmc_individual_rstan_monitor
  
  q2_pooled_samples <- sapply(1:n_chains, function(k) grhmc_run[[k]]$q_original_samples[, 2])
  q1_squared_plus_q2_squared_pooled_samples <- 
    sapply(1:n_chains, function(k) rowSums(grhmc_run[[k]]$q_original_samples[, 1:2] ^ 2))
  q2_list_chains <- lapply(1:n_chains, function(k) matrix(grhmc_run[[k]]$q_original_samples[, 2], ncol = 1))
  q2_squared_list_chains <- lapply(1:n_chains, function(k) matrix(grhmc_run[[k]]$q_original_samples[, 2] ^ 2, ncol = 1))
  q1_squared_plus_q2_squared_list_chains <- 
    lapply(1:n_chains, function(k) matrix(rowSums(grhmc_run[[k]]$q_original_samples[, 1:2] ^ 2), ncol = 1))
  mean_q2 <- mean(final_grhmc_individual_rstan_monitor[grep(parameter_names[2], rownames(final_grhmc_individual_rstan_monitor)), ]$mean)
  mean_q2_squared <- mean(q2_pooled_samples ^ 2)
  mean_q1_squared_plus_q2_squared <- mean(q1_squared_plus_q2_squared_pooled_samples)
  se_mean_q2 <- sd(final_grhmc_individual_rstan_monitor[grep(parameter_names[2], rownames(final_grhmc_individual_rstan_monitor)), ]$mean)
  se_mean_q2_squared <- sd(sapply(1:n_chains, function(k) mean(grhmc_run[[k]]$q_original_samples[, 2] ^ 2)))
  se_mean_q1_squared_plus_q2_squared <- 
    sd(sapply(1:n_chains, function(k) mean(rowSums(grhmc_run[[k]]$q_original_samples[, 1:2] ^ 2))))
  sd_q2 <- sd(q2_pooled_samples)
  sd_q2_squared <- sd(q2_pooled_samples ^ 2)
  sd_q1_squared_plus_q2_squared <- sd(q1_squared_plus_q2_squared_pooled_samples)
  mean_original_rstan_se_mean <- 
    sqrt(mean(final_grhmc_individual_rstan_monitor[grep(parameter_names[2], rownames(final_grhmc_individual_rstan_monitor)), ]$se_mean ^ 2))
  # batch_ess_q <- ((sd_q / se_mean_q) ^ 2) * n_chains
  if (i <= 6) {
    batch_ess_q2 <- stableGR::n.eff(q2_list_chains, multivariate = F, size = NULL)$n.eff / n_chains # must divide by the number of chains manually here due to a bug in the package
    batch_ess_q2_per_mill_eval_ode <- batch_ess_q2 / grhmc_n_evals_ode * 1e6 
    batch_ess_q2_squared <- stableGR::n.eff(q2_squared_list_chains, multivariate = F, size = NULL)$n.eff / n_chains
    batch_ess_q2_squared_per_mill_eval_ode <- batch_ess_q2_squared / grhmc_n_evals_ode * 1e6 
    batch_ess_q1_squared_plus_q2_squared <- stableGR::n.eff(q1_squared_plus_q2_squared_list_chains, multivariate = F, size = NULL)$n.eff / n_chains
    batch_ess_q1_squared_plus_q2_squared_per_mill_eval_ode <- batch_ess_q1_squared_plus_q2_squared / grhmc_n_evals_ode * 1e6 
  } else {
    batch_ess_q2 <- stableGR::n.eff(q2_list_chains, multivariate = F, size = 1000)$n.eff / n_chains # must divide by the number of chains manually here due to a bug in the package
    batch_ess_q2_per_mill_eval_ode <- batch_ess_q2 / grhmc_n_evals_ode * 1e6 
    batch_ess_q2_squared <- stableGR::n.eff(q2_squared_list_chains, multivariate = F, size = 1000)$n.eff / n_chains
    batch_ess_q2_squared_per_mill_eval_ode <- batch_ess_q2_squared / grhmc_n_evals_ode * 1e6 
    batch_ess_q1_squared_plus_q2_squared <- stableGR::n.eff(q1_squared_plus_q2_squared_list_chains, multivariate = F, size = 1000)$n.eff / n_chains
    batch_ess_q1_squared_plus_q2_squared_per_mill_eval_ode <- batch_ess_q1_squared_plus_q2_squared / grhmc_n_evals_ode * 1e6 
  }
  final_grhmc_batch_ess <- data.frame(
    mean_q2 = mean_q2,
    mean_q2_squared = mean_q2_squared,
    mean_q1_squared_plus_q2_squared = mean_q1_squared_plus_q2_squared,
    sd_q2 = sd_q2,
    sd_q2_squared = sd_q2_squared,
    sd_q1_squared_plus_q2_squared = sd_q1_squared_plus_q2_squared,
    se_mean_q2 = se_mean_q2,
    se_mean_q2_squared = se_mean_q2_squared,
    se_mean_q1_squared_plus_q2_squared = se_mean_q1_squared_plus_q2_squared,
    mean_original_rstan_se_mean = mean_original_rstan_se_mean,
    batch_ess_q2 = batch_ess_q2,
    batch_ess_q2_per_mill_eval_ode = batch_ess_q2_per_mill_eval_ode,
    batch_ess_q2_squared = batch_ess_q2_squared,
    batch_ess_q2_squared_per_mill_eval_ode = batch_ess_q2_squared_per_mill_eval_ode,
    batch_ess_q1_squared_plus_q2_squared = batch_ess_q1_squared_plus_q2_squared,
    batch_ess_q1_squared_plus_q2_squared_per_mill_eval_ode = batch_ess_q1_squared_plus_q2_squared_per_mill_eval_ode,
    # n_boundary_crossings = grhmc_n_boundary_crossings,
    n_evals_ode = grhmc_n_evals_ode, 
    n_dim = n_dim,
    parameter_names = parameter_names[2]
  )
  
  rownames(final_grhmc_batch_ess) <- parameter_names[2]
  full_final_grhmc_batch_ess <- rbind(full_final_grhmc_batch_ess, final_grhmc_batch_ess)
  full_final_result_matrix <- rbind(full_final_result_matrix, final_individual_result_matrix)
  
}

full_final_diff_result_matrix <- 
  sweep(full_final_result_matrix, 2, c(true_mean_q2$value, true_mean_q2_squared$value, true_mean_q1_squared_plus_q2_squared$value, 0), "-")

saveRDS(
  list(
    full_final_grhmc_batch_ess = full_final_grhmc_batch_ess,
    full_final_result_matrix = full_final_result_matrix,
    full_final_diff_result_matrix = full_final_diff_result_matrix
  ),
  paste0("disc_grad/example_scripts/toy_model/max_model_varying_dimensions_grhmc_results_with_batch_ess_results.RDS")
)

grhmc_results <- readRDS(
  paste0("disc_grad/example_scripts/toy_model/max_model_varying_dimensions_grhmc_results_with_batch_ess_results.RDS")
)

full_final_grhmc_batch_ess <- grhmc_results$full_final_grhmc_batch_ess
full_final_grhmc_result_matrix <- grhmc_results$full_final_result_matrix
full_final_grhmc_diff_result_matrix <- grhmc_results$full_final_diff_result_matrix
