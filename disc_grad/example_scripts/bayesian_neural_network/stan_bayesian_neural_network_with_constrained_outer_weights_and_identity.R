rm(list = ls())

library(rstan)

set.seed(42)

n_cores <- 10
n_iterations <- 10
n_samples_per_iteration <- 100000
time_period <- 10000

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

# STAN

bnn_stan_run <- rstan::stan(
  "disc_grad/example_scripts/bayesian_neural_network/stan_bnn_constrained_outer_weights.stan",
  data = list(
    n = n, 
    p = ncol(x),
    x = x,
    y = y
  ),
  chains = n_iterations, 
  iter = 2 * n_samples_per_iteration,
  cores = n_cores / 2,
  seed = 42 # currently 42, try 420
)

saveRDS(bnn_stan_run, "disc_grad/example_scripts/bayesian_neural_network/stan_identity_run_example.RDS")
