rm(list = ls())

library(rstan)
library(dplyr)
library(ggplot2)

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

n_dim_vec <- 2 ^ {1:11}
n_chains <- 100
n_samples <- 10000
n_cores <- 10

dir.create("disc_grad/example_scripts/toy_model/stan", showWarnings = F)

for (k in 1:length(dim_vec)) {
  
  store_chains <- rstan::stan(
    file = "disc_grad/example_scripts/toy_model/max_model.stan",
    data = list(
      d = n_dim,
      c = 1
    ),
    iter = 2 * n_samples,
    seed = 42,
    chains = n_chains,
    cores = n_cores,
    control = list(
      metric = "unit_e"
    )
  )

  saveRDS(
    store_chains, 
    paste0("disc_grad/example_scripts/toy_model/stan/max_model_d_", n_dim, "_stan.RDS")
  )  
  
}

