rm(list = ls())

library(doParallel)
library(progress)
library(dplyr)
library(ggplot2)

source("disc_grad/illustrations/illustration_global_error/general_scripts/grhmc_discontinuous_gradient_fixed_step_size_transformed_function.R")

set.seed(42)

alpha <- 0
delta_1 <- 0.5 
delta_2 <- -0.5
beta_1 <- c(1, 0)
beta_2 <- c(-0.1, 1)
w1 <- w2 <- 1
sigma <- 0.1
gamma <- log(sigma ^ 2)
# n <- 3
n <- 100
n_test <- 100
x <- matrix(rnorm(n * 2), ncol = 2)
h_vec <- 1 / seq(1000, 100000, by = 1000)

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

q_initial <- c(gamma, alpha, log(c(w1, w2)), c(delta_1, delta_2), beta_1, beta_2)
p_initial <- rep(c(0.1, -0.1), 5)
z_initial <- c(q_initial, p_initial)
time_period <- 0.5

## Use LSODAR with a very high precision to estimate the "true" solution at final time point ##

lsodar_ode_func <- function(t, z, parms) {
  
  # z[1] <- gamma = log(sigma^2)
  # z[2] <- alpha
  # z[3:4] <- w* = log(w) in order to get w >= 0
  # z[5:6] <- delta_j
  # z[7:8] <- beta_1
  # z[9:10] <- beta_2
  # z[11:20] <- p
  
  neuron_1 <- (z[5] + x %*% z[7:8]) * lsodar_neuron_indices_list[[1]]
  neuron_2 <- (z[6] + x %*% z[9:10]) * lsodar_neuron_indices_list[[2]]
  exp_gamma <- exp(z[1])
  w <- exp(z[3:4])
  mu <- z[2] + cbind(neuron_1, neuron_2) %*% w
  
  diff_squared_full <- (y - mu) ^ 2
  diff_full <- y - mu
  diff_with_x1_full <- (y - mu) * x[, 1]
  diff_with_x2_full <- (y - mu) * x[, 2]
  diff_with_neuron_1 <- (y - mu) * neuron_1
  diff_with_neuron_2 <- (y - mu) * neuron_2
  
  sum_diff_squared_full <- sum(diff_squared_full)
  sum_diff_full <- sum(diff_full)
  sum_diff_with_x1_full <- sum(diff_with_x1_full)
  sum_diff_with_x2_full <- sum(diff_with_x2_full)
  sum_diff_neuron_1 <- sum(diff_with_neuron_1)
  sum_diff_neuron_2 <- sum(diff_with_neuron_2)
  
  ret <- c(
    z[11:20],
    - n / 2 + (1 / (2 * exp_gamma)) * sum_diff_squared_full - sqrt(exp_gamma) / 2 + 1 / 2, # grad wrt gamma
    - (z[2] - prior_mu[1]) / (prior_sigma[1] ^ 2) + (1 / exp_gamma) * sum_diff_full, # grad wrt alpha
    - (z[3:4] - prior_mu[2:3]) / (prior_sigma[2:3] ^ 2) + c(sum_diff_neuron_1, sum_diff_neuron_2) / exp_gamma * w[1:2], # grad wrt w
    - (z[5:6] - prior_mu[4:5]) / (prior_sigma[4:5] ^ 2)  + c(sum(diff_full * lsodar_neuron_indices_list[[1]]), sum(diff_full * lsodar_neuron_indices_list[[2]])) / exp_gamma * w, # grad wrt delta_j
    - (z[7:8] - prior_mu[6:7]) / (prior_sigma[6:7] ^ 2) + c(sum(diff_with_x1_full * lsodar_neuron_indices_list[[1]]), sum(diff_with_x2_full * lsodar_neuron_indices_list[[1]])) / exp_gamma * w[1], # grad wrt beta_1
    - (z[9:10] - prior_mu[8:9]) / (prior_sigma[8:9] ^ 2)  + c(sum(diff_with_x1_full * lsodar_neuron_indices_list[[2]]), sum(diff_with_x2_full * lsodar_neuron_indices_list[[2]])) / exp_gamma * w[2] # grad wrt beta_2
  )
  
  list(ret)
  
}

lsodar_root_func <- function(t, z, parms) {
  c(
    z[5] + x %*% z[7:8],
    z[6] + x %*% z[9:10] 
  )
}

lsodar_event_func <- function(t, z, parms) {
  
  if (t > 0) {
    old_q <- z[1:10]
    p <- z[11:20]
    new_q <- old_q + 1e-5 * p # explictly push q to the new region
    new_z <- c(new_q, p)
    lsodar_neuron_indices_list <<- list(
      as.integer(new_z[5] + x %*% new_z[7:8] > 0), # observations that will contribute to the first neuron
      as.integer(new_z[6] + x %*% new_z[9:10] > 0) # observations that will contribute to the second neuron
    ) 
  }
  
  z
}

neuron_1_initial <- pmax(0, z_initial[5] + x %*% z_initial[7:8])
neuron_2_initial <- pmax(0, z_initial[6] + x %*% z_initial[9:10])
neuron_1_indices_initial <- as.integer(neuron_1_initial > 0)
neuron_2_indices_initial <- as.integer(neuron_2_initial > 0)
lsodar_neuron_indices_list <<- list(
  neuron_1_indices_initial,
  neuron_2_indices_initial
)

lsodar_neuron_indices_list

lsodar_run <- deSolve::lsodar(
  y = z_initial, 
  times = seq(from = 0, to = time_period, by = 0.01),
  func = lsodar_ode_func,
  events = list(func = lsodar_event_func, root = TRUE),
  rootfunc = lsodar_root_func,
  rtol = 1e-12,
  atol = 1e-12,
  maxsteps = 1000000
)

z_true <- lsodar_run[nrow(lsodar_run), 2:21]
z_true

## Standard leapfrog - adapted from Neal et al. (2011) ## 

# No detecting event, gradient can change during a single leapfrog step

leapfrog <- function(grad_U, h, L, q_initial, p_initial) {
  
  q <- q_initial
  p <- p_initial
  
  p <- p - h * grad_U(q) / 2 # half step 
  for (i in 1:L)
  {
    q <- q + h * p
    if (i != L) {
      p <- p - h * grad_U(q) # full step
    }
  }
  p <- p - h * grad_U(q) / 2 # half step
  list(q = q, p = p)
  
}

leapfrog_grad_U <- function(q) { # U = negative log target, therefore negative log gradient here
  # q[1] <- gamma = log(sigma^2)
  # q[2] <- alpha
  # q[3:4] <- w* = log(w) in order to get w >= 0
  # q[5:6] <- delta_j
  # q[7:8] <- beta_1
  # q[9:10] <- beta_2

  neuron_1 <- pmax(0, q[5] + x %*% q[7:8])
  neuron_2 <- pmax(0, q[6] + x %*% q[9:10])
  neuron_1_indices <- which(neuron_1 > 0)
  neuron_2_indices <- which(neuron_2 > 0)
  
  exp_gamma <- exp(q[1])
  w <- exp(q[3:4])
  mu <- q[2] + cbind(neuron_1, neuron_2) %*% w
  
  diff_squared_full <- (y - mu) ^ 2
  diff_full <- y - mu
  diff_with_x1_full <- (y - mu) * x[, 1]
  diff_with_x2_full <- (y - mu) * x[, 2]
  diff_with_neuron_1 <- (y - mu) * neuron_1
  diff_with_neuron_2 <- (y - mu) * neuron_2
  
  sum_diff_squared_full <- sum(diff_squared_full)
  sum_diff_full <- sum(diff_full)
  sum_diff_with_x1_full <- sum(diff_with_x1_full)
  sum_diff_with_x2_full <- sum(diff_with_x2_full)
  sum_diff_neuron_1 <- sum(diff_with_neuron_1)
  sum_diff_neuron_2 <- sum(diff_with_neuron_2)
  
  ret <- -c(
    - n / 2 + (1 / (2 * exp_gamma)) * sum_diff_squared_full - sqrt(exp_gamma) / 2 + 1 / 2, # grad wrt gamma
    - (q[2] - prior_mu[1]) / (prior_sigma[1] ^ 2) + (1 / exp_gamma) * sum_diff_full, # grad wrt alpha
    - (q[3:4] - prior_mu[2:3]) / (prior_sigma[2:3] ^ 2) + c(sum_diff_neuron_1, sum_diff_neuron_2) / exp_gamma * w[1:2], # grad wrt w
    - (q[5:6] - prior_mu[4:5]) / (prior_sigma[4:5] ^ 2)  + c(sum(diff_full[neuron_1_indices]), sum(diff_full[neuron_2_indices])) / exp_gamma * w, # grad wrt delta_j
    - (q[7:8] - prior_mu[6:7]) / (prior_sigma[6:7] ^ 2) + c(sum(diff_with_x1_full[neuron_1_indices]), sum(diff_with_x2_full[neuron_1_indices])) / exp_gamma * w[1], # grad wrt beta_1
    - (q[9:10] - prior_mu[8:9]) / (prior_sigma[8:9] ^ 2)  + c(sum(diff_with_x1_full[neuron_2_indices]), sum(diff_with_x2_full[neuron_2_indices])) / exp_gamma * w[2] # grad wrt beta_2
  )
  
  ret
}

leapfrog_output <- function(grad_U, time_period, q_initial, p_initial, z_true, h) {
  
  L <- time_period / h
  
  # if (abs(L - floor(L)) > 1.0e-9) {
  #   stop("Please provide epsilon so that time_period / h is an integer!")
  # }
  
  if (abs(L - round(L)) > 1.0e-9) {
    stop("Please provide epsilon so that time_period / h is an integer!")
  }
  
  run_leapfrog <- leapfrog(grad_U, h, L, q_initial, p_initial)
  
  q_leapfrog_solution <- run_leapfrog$q
  p_leapfrog_solution <- run_leapfrog$p
  
  list(q = q_leapfrog_solution, p = p_leapfrog_solution, l2_error = sqrt(sum((c(q_leapfrog_solution, p_leapfrog_solution) - z_true) ^ 2)))
  
}

# Test 

leapfrog_output(grad_U = leapfrog_grad_U, time_period = time_period, q_initial = q_initial, p_initial = p_initial, z_true = z_true, h = 0.01)  
leapfrog_output(grad_U = leapfrog_grad_U, time_period = time_period, q_initial = q_initial, p_initial = p_initial, z_true = z_true, h = 0.005)  

# Run for a range of values for the step size h

leapfrog_store_error <- numeric(length(h_vec))

# for (i in 1:length(h_vec)) {
#   
#   leapfrog_store_error[i] <- leapfrog_output(grad_U = leapfrog_grad_U, time_period = time_period, q_initial = q_initial, p_initial = p_initial, z_true = z_true, h = h_vec[i])$l2_error  
#   
# }

init_cluster <- parallel::makeCluster(10) # if 10 cores, change if not
doParallel::registerDoParallel(init_cluster)
parallel::clusterExport(
  init_cluster,
  c(
    "x",
    "y",
    "n",
    "prior_mu",
    "prior_sigma"
  )
)

leapfrog_store_error <- foreach::foreach(i = 1:length(h_vec), .combine = "c") %dopar% {
  leapfrog_output(grad_U = leapfrog_grad_U, time_period = time_period, q_initial = q_initial, p_initial = p_initial, z_true = z_true, h = h_vec[i])$l2_error
}

leapfrog_store_matrix <- data.frame(method = rep("Leapfrog", length(h_vec)), h = h_vec, l2_error = leapfrog_store_error)

## GRHMC - Gradient can change during a single integration step ##

A <- matrix(nrow = 2 * n, ncol = 10)

for (i in 1:n) { # related to first neuron 
  
  A[i, ] <- c(rep(0, 4), 1, 0, x[i,  ], rep(0, 2))
  
}

for (i in 1:n) { # related to second neuron
  
  A[n + i, ] <- c(rep(0, 4), 0, 1, rep(0, 2), x[i, ])
  
}

A

B <- rep(0, 2 * n)

grhmc_model_nr1 <- list(
  
  n_parameters = 10,
  
  grad_jump_fun = function(q) {
    
    NULL
    
  },
  
  region_non_lin_root_list = list(root_fun = NULL, event_fun = NULL),
  
  region_lin_root_list = list(A = A, B = B),
  
  sim_q0 = function() rep(1, 10)
  
)

grhmc_model_nr1$log_target_grad <- function(q) {
  # q[1] <- gamma = log(sigma^2)
  # q[2] <- alpha
  # q[3:4] <- w* = log(w) in order to get w >= 0
  # q[5:6] <- delta_j
  # q[7:8] <- beta_1
  # q[9:10] <- beta_2
  
  neuron_1 <- pmax(0, q[5] + x %*% q[7:8])
  neuron_2 <- pmax(0, q[6] + x %*% q[9:10])
  neuron_1_indices <- which(neuron_1 > 0)
  neuron_2_indices <- which(neuron_2 > 0)
  
  exp_gamma <- exp(q[1])
  w <- exp(q[3:4])
  mu <- q[2] + cbind(neuron_1, neuron_2) %*% w
  
  diff_squared_full <- (y - mu) ^ 2
  diff_full <- y - mu
  diff_with_x1_full <- (y - mu) * x[, 1]
  diff_with_x2_full <- (y - mu) * x[, 2]
  diff_with_neuron_1 <- (y - mu) * neuron_1
  diff_with_neuron_2 <- (y - mu) * neuron_2
  
  sum_diff_squared_full <- sum(diff_squared_full)
  sum_diff_full <- sum(diff_full)
  sum_diff_with_x1_full <- sum(diff_with_x1_full)
  sum_diff_with_x2_full <- sum(diff_with_x2_full)
  sum_diff_neuron_1 <- sum(diff_with_neuron_1)
  sum_diff_neuron_2 <- sum(diff_with_neuron_2)
  
  c(
    - n / 2 + (1 / (2 * exp_gamma)) * sum_diff_squared_full - sqrt(exp_gamma) / 2 + 1 / 2, # grad wrt gamma
    - (q[2] - prior_mu[1]) / (prior_sigma[1] ^ 2) + (1 / exp_gamma) * sum_diff_full, # grad wrt alpha
    - (q[3:4] - prior_mu[2:3]) / (prior_sigma[2:3] ^ 2) + c(sum_diff_neuron_1, sum_diff_neuron_2) / exp_gamma * w[1:2], # grad wrt w
    - (q[5:6] - prior_mu[4:5]) / (prior_sigma[4:5] ^ 2)  + c(sum(diff_full[neuron_1_indices]), sum(diff_full[neuron_2_indices])) / exp_gamma * w, # grad wrt delta_j
    - (q[7:8] - prior_mu[6:7]) / (prior_sigma[6:7] ^ 2) + c(sum(diff_with_x1_full[neuron_1_indices]), sum(diff_with_x2_full[neuron_1_indices])) / exp_gamma * w[1], # grad wrt beta_1
    - (q[9:10] - prior_mu[8:9]) / (prior_sigma[8:9] ^ 2)  + c(sum(diff_with_x1_full[neuron_2_indices]), sum(diff_with_x2_full[neuron_2_indices])) / exp_gamma * w[2] # grad wrt beta_2
  )
  
}

grhmc_output <- function(model_list, time_period, q_initial, p_initial, z_true, h) {
  
  grad_indices <- model_list$grad_jump_fun(q_initial)
  
  grhmc_run <- grhmc_discontinuous_gradient_fixed_step_size_transformed_function(
    model_list = model_list,
    T = time_period,
    n_samples = 100,
    diag_s_elements_initial = rep(1, 10),
    m_initial = rep(0, 10),
    qbar_initial = q_initial,
    pbar_initial = p_initial,
    h = h,
    last.root.offset.lin.root.finder = 1.0e-8,
    last.root.offset.non.lin.root.finder = 1.0e-8,
    precision_real_root_lin_root_finder = 1.0e-13,
    num_subdiv_non_lin_root_finder = 8L
  )
  
  q_grhmc_solution <- grhmc_run$q_final 
  p_grhmc_solution <- grhmc_run$p_final
  z_grhmc_solution <- c(q_grhmc_solution, p_grhmc_solution)
  
  list(q = q_grhmc_solution, p = p_grhmc_solution, l2_error = sqrt(sum((z_grhmc_solution - z_true) ^ 2)))
  
}

grhmc_output(model_list = grhmc_model_nr1, time_period = time_period, q_initial = q_initial, p_initial = p_initial, z_true = z_true, h = 0.01)
grhmc_output(model_list = grhmc_model_nr1, time_period = time_period, q_initial = q_initial, p_initial = p_initial, z_true = z_true, h = 0.005)

# Run for a range of values for the step size h

grhmc_nr1_store_error <- numeric(length(h_vec))

# for (i in 1:length(h_vec)) {
#   print(i)
#   grhmc_nr1_store_error[i] <- grhmc_output(model_list = grhmc_model_nr1, time_period = time_period, q_initial = q_initial, p_initial = p_initial, z_true = z_true, h = h_vec[i])$l2_error  
# }

grhmc_nr1_store_error <- foreach::foreach(i = 1:length(h_vec), .combine = "c") %dopar% {
  grhmc_output(model_list = grhmc_model_nr1, time_period = time_period, q_initial = q_initial, p_initial = p_initial, z_true = z_true, h = h_vec[i])$l2_error
}

grhmc_nr1_store_matrix <- data.frame(method = rep("GRHMC nr.1", length(h_vec)), h = h_vec, l2_error = grhmc_nr1_store_error)

## GRHMC - Gradient fixed during a single integration step (Proposed method) ##

grhmc_model_nr2 <- list(
  
  n_parameters = 10,
  
  grad_jump_fun = function(q) {
    
    neuron_indices_list <<- list(
      as.integer(q[5] + x %*% q[7:8] > 0), # observations that will contribute to the first neuron
      as.integer(q[6] + x %*% q[9:10] > 0) # observations that will contribute to the second neuron
    )
    
  },
  
  region_non_lin_root_list = list(root_fun = NULL, event_fun = NULL),
  
  region_lin_root_list = list(A = A, B = B),
  
  sim_q0 = function() rep(1, 10)
  
)

grhmc_model_nr2$log_target_grad <- function(q) {
  # q[1] <- gamma = log(sigma^2)
  # q[2] <- alpha
  # q[3:4] <- w* = log(w) in order to get w >= 0
  # q[5:6] <- delta_j
  # q[7:8] <- beta_1
  # q[9:10] <- beta_2
  
  neuron_1 <- (q[5] + x %*% q[7:8]) * neuron_indices_list[[1]]
  neuron_2 <- (q[6] + x %*% q[9:10]) * neuron_indices_list[[2]]
  exp_gamma <- exp(q[1])
  w <- exp(q[3:4])
  mu <- q[2] + cbind(neuron_1, neuron_2) %*% w
  
  diff_squared_full <- (y - mu) ^ 2
  diff_full <- y - mu
  diff_with_x1_full <- (y - mu) * x[, 1]
  diff_with_x2_full <- (y - mu) * x[, 2]
  diff_with_neuron_1 <- (y - mu) * neuron_1
  diff_with_neuron_2 <- (y - mu) * neuron_2
  
  sum_diff_squared_full <- sum(diff_squared_full)
  sum_diff_full <- sum(diff_full)
  sum_diff_with_x1_full <- sum(diff_with_x1_full)
  sum_diff_with_x2_full <- sum(diff_with_x2_full)
  sum_diff_neuron_1 <- sum(diff_with_neuron_1)
  sum_diff_neuron_2 <- sum(diff_with_neuron_2)
  
  c(
    - n / 2 + (1 / (2 * exp_gamma)) * sum_diff_squared_full - sqrt(exp_gamma) / 2 + 1 / 2, # grad wrt gamma
    - (q[2] - prior_mu[1]) / (prior_sigma[1] ^ 2) + (1 / exp_gamma) * sum_diff_full, # grad wrt alpha
    - (q[3:4] - prior_mu[2:3]) / (prior_sigma[2:3] ^ 2) + c(sum_diff_neuron_1, sum_diff_neuron_2) / exp_gamma * w[1:2], # grad wrt w
    - (q[5:6] - prior_mu[4:5]) / (prior_sigma[4:5] ^ 2)  + c(sum(diff_full * neuron_indices_list[[1]]), sum(diff_full * neuron_indices_list[[2]])) / exp_gamma * w, # grad wrt delta_j
    - (q[7:8] - prior_mu[6:7]) / (prior_sigma[6:7] ^ 2) + c(sum(diff_with_x1_full * neuron_indices_list[[1]]), sum(diff_with_x2_full * neuron_indices_list[[1]])) / exp_gamma * w[1], # grad wrt beta_1
    - (q[9:10] - prior_mu[8:9]) / (prior_sigma[8:9] ^ 2)  + c(sum(diff_with_x1_full * neuron_indices_list[[2]]), sum(diff_with_x2_full * neuron_indices_list[[2]])) / exp_gamma * w[2] # grad wrt beta_2
    
  )
  
}

grhmc_output(model_list = grhmc_model_nr2, time_period = time_period, q_initial = q_initial, p_initial = p_initial, z_true = z_true, h = 0.01)
grhmc_output(model_list = grhmc_model_nr2, time_period = time_period, q_initial = q_initial, p_initial = p_initial, z_true = z_true, h = 0.005)

# Run for a range of values for the step size h

# grhmc_nr2_store_error <- numeric(length(h_vec))
# 
# for (i in 1:length(h_vec)) {
#   print(i)
#   grhmc_nr2_store_error[i] <- grhmc_output(model_list = grhmc_model_nr2, time_period = time_period, q_initial = q_initial, p_initial = p_initial, z_true = z_true, h = h_vec[i])$l2_error  
# }

grhmc_nr2_store_error <- foreach::foreach(i = 1:length(h_vec), .combine = "c") %dopar% {
  grhmc_output(model_list = grhmc_model_nr2, time_period = time_period, q_initial = q_initial, p_initial = p_initial, z_true = z_true, h = h_vec[i])$l2_error
}
parallel::stopCluster(init_cluster)

grhmc_nr2_store_matrix <- data.frame(method = rep("GRHMC nr.2", length(h_vec)), h = h_vec, l2_error = grhmc_nr2_store_error)

## Comparison plot ##

summary(lm(l2_error ~ I(h ^ 3) + I(h ^ 2) + h, data = grhmc_nr2_store_matrix))
range(c(leapfrog_store_error, grhmc_nr1_store_error, grhmc_nr2_store_error))

plot(h_vec, leapfrog_store_error, log = "xy", ylim = c(1e-10, 1), ylab = "L2 error at final state", xlab = "h (step size)")
points(h_vec, grhmc_nr1_store_error, col = "red")
points(h_vec, grhmc_nr2_store_error, col = "blue", pch = 3)
lines(h_vec, 1.2e7 * h_vec ^ 3, lty = 2, col = "green")
legend(
  "bottomright",
  legend = c(
    "Leapfrog",
    "GRHMC nr.1",
    "GRHMC nr.2",
    latex2exp::TeX("$\\propto h^3$")
  ),
  pch = c(1, 1, 3, NA),
  lty = c(NA, NA, NA, 2),
  col = c("black", "red", "blue", "green")
)

plot(h_vec, leapfrog_store_error, log = "xy", ylim = c(1e-10, 1), xlim = c(1e-05, 1.5e-03), ylab = "L2 error at final state", xlab = "h (step size)")
points(h_vec, grhmc_nr1_store_error, col = "red")
points(h_vec, grhmc_nr2_store_error, col = "blue", pch = 3)
lines(h_vec, 1.2e7 * h_vec ^ 3, lty = 2, col = "green")
legend(
  "bottomright",
  legend = c(
    "Leapfrog",
    "GRHMC nr.1",
    "GRHMC nr.2",
    latex2exp::TeX("$\\propto h^3$")
  ),
  pch = c(1, 1, 3, NA),
  lty = c(NA, NA, NA, 2),
  col = c("black", "red", "blue", "green")
)

pdf("disc_grad/example_scripts/section_6/bayesian_neural_network/illustration_global_error/bayesian_neural_network_comparison_plot.pdf",
    width = 11.7, height = 8.3)

relevant_h_indices <- h_vec >= 2e-5
plot(
  h_vec[relevant_h_indices], 
  leapfrog_store_error[relevant_h_indices], 
  log = "xy", 
  ylim = c(1e-10, 1), 
  xlim = c(2e-05, 1.5e-03), 
  ylab = "L2 error at final state", 
  xlab = "h (step size)",
  cex.lab = 1.5,
  cex.axis = 1.5
)
points(h_vec[relevant_h_indices], grhmc_nr1_store_error[relevant_h_indices], col = "red")
points(h_vec[relevant_h_indices], grhmc_nr2_store_error[relevant_h_indices], col = "blue", pch = 3)
lines(h_vec, 1.2e7 * h_vec ^ 3, lty = 2, col = "green")
legend(
  "bottomright",
  legend = c(
    "Leapfrog",
    "GRHMC without boundary handling",
    "GRHMC with boundary handling",
    latex2exp::TeX("$\\propto h^3$")
  ),
  pch = c(1, 1, 3, NA),
  lty = c(NA, NA, NA, 2),
  col = c("black", "red", "blue", "green"), 
  cex = 1.5
)

dev.off()

final_store_matrix <- rbind(
  leapfrog_store_matrix, 
  grhmc_nr1_store_matrix,
  grhmc_nr2_store_matrix
)

final_store_matrix %>%
  ggplot(aes(x = h, y = l2_error)) + 
  geom_point(aes(colour = method, shape = method)) + 
  scale_shape_manual(values = c(1, 3, 1), guide = "none") + 
  geom_function(aes(colour = "Line"), fun = function(x) 0.25 * x ^ 3, linetype = 2) + 
  scale_color_manual(
    values = c("red", "blue", "black", "green"),
    labels = c("GRHMC nr.1", "GRHMC nr.2", "Leapfrog", expression("" %prop% h^3))
  ) + 
  guides(
    color = guide_legend(override.aes = list(shape = c(1, 3, 1, NA), linetype = c(NA, NA, NA, 2)))
  ) + 
  scale_x_continuous(trans = "log") + 
  scale_y_continuous(trans = "log") + 
  labs(colour = "") + 
  xlab("h (step size)") + 
  ylab("L2 error at final state") + 
  theme(
    legend.text.align = 0
  )
