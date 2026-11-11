rm(list = ls())
library(doParallel)

source("disc_grad/implementation_scripts/general_scripts/cpp/grhmc_discontinuous_gradient_with_relevant_indices_function_cpp.R")

grhmc_wrapper_file_path <- normalizePath(
  "disc_grad/implementation_scripts/general_scripts/cpp/grhmc_discontinuous_gradient_with_relevant_indices_function_cpp.R"
)
model_file_path <- normalizePath(
  "disc_grad/example_scripts/toy_model/spike_model/grhmc/spike_model.cpp"
)
model_cache <- normalizePath(
  "disc_grad/example_scripts/toy_model/spike_model/grhmc/spike_model_cache",
  mustWork = F
)
progress_log <- normalizePath("disc_grad/example_scripts/toy_model/spike_model/grhmc/progress.log", mustWork = FALSE)
if (file.exists(progress_log)) file.remove(progress_log)
dir.create(model_cache, showWarnings = FALSE)
Rcpp::sourceCpp(model_file_path, cacheDir = model_cache, rebuild = T)

set.seed(42)
d <- 2 ^ 5
delta <- 0.01
eta <- 1
delta_string <- sub(".", "_", delta, fixed = T)
n_cores <- 10
n_trajectories <- 10
trajectory_length <- 100000
n_samples <- 1000000

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
  .noexport = c("grhmc_spike_model_cpp",
                "grhmc_discontinuous_gradient_function_cpp")
) %dopar% {
  n_dim <- d
  q_initial <- rnorm(n_dim)
  p_initial <- rnorm(n_dim)
  u_initial <- rexp(1)
  
  test_run_warmup <- grhmc_discontinuous_gradient_function_cpp(
    model_list = function(settings) grhmc_spike_model_cpp(settings, delta = delta, eta = 1),
    lambda = 0.2,
    T = trajectory_length,
    n_samples = n_samples,
    q_initial,
    p_initial,
    u_initial, 
    random_state = NULL,
    rtol = 1e-6,
    atol = 1e-6,
    h_max = 1.0,
    last_root_offset_lin_root_finder = 1e-10,
    precision_real_root_lin_root_finder = 1.0e-13,
    coeff_tol_real_root_lin_root_finder = 1.0e-13,
    push_eps_tol = 1e-4,
    relevant_indices = c(1, n_dim),
    return_output_from_ode = FALSE,
    trajectory_info_string = paste0("warmup ", i),
    progress_ratio = 0.05
  )
  
  new_p_initial <- rnorm(n_dim)
  new_u_initial <- rexp(1)
  
  test_run_sampling <- grhmc_discontinuous_gradient_function_cpp(
    model_list = function(settings) grhmc_spike_model_cpp(settings, delta = delta, eta = 1),
    lambda = 0.2,
    T = trajectory_length,
    n_samples = n_samples,
    test_run_warmup$q_final,
    new_p_initial,
    new_u_initial, 
    random_state = NULL,
    rtol = 1e-6,
    atol = 1e-6,
    h_max = 1.0,
    last_root_offset_lin_root_finder = 1e-10,
    precision_real_root_lin_root_finder = 1.0e-13,
    coeff_tol_real_root_lin_root_finder = 1.0e-13,
    push_eps_tol = 1e-4,
    relevant_indices = c(1, n_dim),
    return_output_from_ode = FALSE,
    trajectory_info_string = paste0("sampling ", i),
    progress_ratio = 0.05
  )
  test_run_sampling
}

saveRDS(
  store_chains, 
  paste0("disc_grad/example_scripts/toy_model/spike_model/grhmc/delta_", delta_string, "/spike_model_delta_", delta_string, "_dim_", d, ".RDS")
)

