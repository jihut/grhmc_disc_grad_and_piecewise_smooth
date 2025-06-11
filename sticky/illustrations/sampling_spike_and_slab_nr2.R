source("sticky/implementation_scripts/general_scripts/grhmc_sticky_transformed_function.R")

# Here: Only five coordinates follow spike and slab distribution, the remaining follow a normal distribution

n_samples_per_iteration <- 100000
time_period <- 100000

set.seed(42)

p <- 10

sigma_slab <- 4

grhmc_model_list <- list(
  
  n_parameters = p,
  
  sticky_indices = seq(from = 2, to = p, by = 2), 
  
  # sim_q0 = function() rnorm(2 * length_beta + 1)
  sim_q0 = function() c(rnorm(p, mean = mu, sd = rho))
  
)

grhmc_model_list$weight_slab <- rep(0.25, length(grhmc_model_list$sticky_indices))
grhmc_model_list$sigma_slab <- rep(sigma_slab, length(grhmc_model_list$sticky_indices))

grhmc_model_list$log_target_grad <- function(q) { # new version, vectorized and therefore slightly faster
  
  log_target_grad_vec <- numeric(grhmc_model_list$n_parameters)
  log_target_grad_vec[grhmc_model_list$sticky_indices] <- 
    -q[grhmc_model_list$sticky_indices] / grhmc_model_list$sigma_slab ^ 2
  log_target_grad_vec[-grhmc_model_list$sticky_indices] <- 
    -(q[-grhmc_model_list$sticky_indices] - 4) / (5 ^ 2)
  log_target_grad_vec
  
}

set.seed(42)

qbar_initial <- rnorm(grhmc_model_list$n_parameters)
pbar_initial <- rnorm(grhmc_model_list$n_parameters)

system.time(
  
  boston_scale_fit_grhmc <- grhmc_sticky_transformed_function(
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

output <- boston_scale_fit_grhmc$output_from_ode_solver

plot(boston_scale_fit_grhmc$q_original_samples[, 1])
plot(boston_scale_fit_grhmc$q_original_samples[, 2])

sapply(1:p, function(i) mean(boston_scale_fit_grhmc$q_original_samples[, i] == 0))
sapply(1:p, function(i) var(boston_scale_fit_grhmc$q_original_samples[boston_scale_fit_grhmc$q_original_samples[, i] != 0, i]))
