rm(list = ls())

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

# Stan 

stan_run <- readRDS("disc_grad/example_scripts/bayesian_neural_network/stan_run_example.RDS")
stan_monitor <- rstan::monitor(stan_run)
rownames(stan_monitor)
n_eff_stan <- stan_monitor[, "n_eff"]
n_eff_stan
n_evals_ode_stan <- sum(rstan::get_num_leapfrog_per_iteration(stan_run))
n_evals_ode_stan
n_eff_per_eval_ode_stan <- n_eff_stan / n_evals_ode_stan
n_eff_per_eval_ode_stan

# GRHMC 

grhmc_disc_grad_run <- readRDS("disc_grad/example_scripts/bayesian_neural_network/bnn_disc_grad_run_VARI_with_n_evals_ode_run_example.RDS")

grhmc_disc_grad_array <- array(dim = c(n_samples_per_iteration, n_iterations, ncol(grhmc_disc_grad_run[[1]])))

for (i in 1:n_iterations) {
  
  grhmc_disc_grad_array[, i, ] <- grhmc_disc_grad_run[[i]]
  
}

grhmc_disc_grad_monitor <- rstan::monitor(grhmc_disc_grad_array, warmup = 0)
n_eff_grhmc_disc_grad <- grhmc_disc_grad_monitor$n_eff
n_eff_grhmc_disc_grad
n_evals_ode_grhmc_disc_grad <- 0
for (i in 1:n_iterations) {
  n_evals_ode_grhmc_disc_grad <- n_evals_ode_grhmc_disc_grad + 
    grhmc_disc_grad_run[, 2][[i]]
}
n_eff_per_eval_ode_grhmc_disc_grad <- n_eff_grhmc_disc_grad / n_evals_ode_grhmc_disc_grad
n_eff_per_eval_ode_grhmc_disc_grad
