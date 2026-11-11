rm(list = ls())

n_chains <- 10
n_cores <- 10
n_samples <- 10000

delta <- 0.01
delta_string <- sub(".", "_", delta, fixed = T)
eta <- 1
n_dim <- 2 ^ 5

store_chains <- rstan::stan(
  file = "disc_grad/example_scripts/toy_model/spike_model/stan/spike_model.stan",
  data = list(
    d = n_dim,
    eta = eta, 
    delta = delta
  ),
  iter = 2 * n_samples,
  seed = 42,
  chains = n_chains,
  cores = n_cores,
  control = list(
    metric = "unit_e",
    max_treedepth = 15
  ),
  pars = c("q1", "q2")
)
store_chains
sum(rstan::get_num_leapfrog_per_iteration(store_chains)) / 1e6

saveRDS(
  store_chains, 
  paste0("disc_grad/example_scripts/toy_model/spike_model/stan/delta_", delta_string, "/spike_model_delta_", delta_string, "_dim_", d, ".RDS")
) 
