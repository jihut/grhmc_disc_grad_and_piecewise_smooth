rm(list = ls())

library(doParallel)
library(dplyr)
library(ggplot2)
source("disc_grad/implementation_scripts/ISG/full_ISG_grhmc_discontinuous_gradient_transformed_function.R")

# q1 ~ N(0,1)
# q2 ~ N(max(0, q1), 1)
max_model <- list(
  
  n_parameters = 2,
  
  grad_jump_fun = function(q) {
    
    grad_indices <<- c(
      as.integer(q[1] > 0)
    )
    
  },
  
  additional_non_lin_root_list = NULL, # usually: list(root_fun = ..., event_fun = ...)
  
  additional_lin_root_list = NULL, # usually: list(A = ..., B = ..., event_fun = ...)
  
  region_lin_root_list = list(A = matrix(c(1, 0), nrow = 1), B = 0), # also need to specify the break point here again
  
  sim_q0 = function() rnorm(2)
  
)

max_model$log_target_grad <- function(q) { # define the log target grad
  # grad_indices = 1 if q1 > 0
  c(
    -q[1] + (q[2] - q[1]) * grad_indices,
    -q[2] + q[1] * grad_indices
  )
  
}

model_list <- max_model

set.seed(1)
qbar_initial <- model_list$sim_q0()
pbar_initial <- rnorm(2)
u_initial <- rexp(1)

grad_indices <- model_list$grad_jump_fun(qbar_initial)

system.time(
  test_run <- ISG_grhmc_discontinuous_gradient_transformed_function(
    model_list = model_list,
    lambda_initial = 0.2,
    time_period_adaptive = 10000,
    time_period_generating_samples = 100000,
    n_generated_samples = 100000,
    n_adaptive_samples = 10000,
    diag_s_elements_initial = rep(1, 2),
    m_initial = rep(0, 2),
    qbar_initial = qbar_initial,
    pbar_initial = pbar_initial, 
    random_state = NULL, 
    rtol = NULL,
    atol = NULL,
    verbose_at_refresh = T
  )
)

test_run$adaptive_step_run

test_run$sample_step_run$output_from_ode_solver
test_run$sample_step_run$n_evals_ode

hist(test_run$sample_step_run$q_original_samples[, 1])
range(test_run$sample_step_run$q_original_samples[, 1])
mean(test_run$sample_step_run$q_original_samples[, 1])
sd(test_run$sample_step_run$q_original_samples[, 1])

hist(test_run$sample_step_run$q_original_samples[, 2])
range(test_run$sample_step_run$q_original_samples[, 2])
mean(test_run$sample_step_run$q_original_samples[, 2])
sd(test_run$sample_step_run$q_original_samples[, 2])

# Check against using rnorm #

q1 <- rnorm(100000)
q2_given_q1 <- rnorm(100000, mean = pmax(0, q1), sd = 1)

par(mfrow = c(1, 2))
hist(test_run$sample_step_run$q_original_samples[, 2])
hist(q2_given_q1)
par(mfrow = c(1, 1))

#####################

# Run 10 independent trajectories and look at the samples generated all together. 
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

init_cluster <- parallel::makeCluster(5)
doParallel::registerDoParallel(init_cluster)
doRNG::registerDoRNG(seed = 42)

store_samples_matrix <- foreach::foreach(i = 1:10, .combine = rbind) %dopar% {
  
  max_model <- list(
    
    n_parameters = 2,
    
    grad_jump_fun = function(q) {
      
      grad_indices <<- c(
        as.integer(q[1] > 0)
      )
      
    },
    
    additional_non_lin_root_list = NULL, # usually: list(root_fun = ..., event_fun = ...)
    
    additional_lin_root_list = NULL, # usually: list(A = ..., B = ..., event_fun = ...)
    
    region_lin_root_list = list(A = matrix(c(1, 0), nrow = 1), B = 0), # also need to specify the break point here again
    
    sim_q0 = function() rnorm(2)
    
  )
  
  max_model$log_target_grad <- function(q) { # define the log target grad
    
    c(
      -q[1] + (q[2] - q[1]) * grad_indices,
      -q[2] + q[1] * grad_indices
    )
    
  }
  
  model_list <- max_model
  
  qbar_initial <- model_list$sim_q0()
  pbar_initial <- rnorm(2)
  
  test_run <- ISG_grhmc_discontinuous_gradient_transformed_function(
    model_list = model_list,
    lambda_initial = 0.2,
    time_period_adaptive = 10000,
    time_period_generating_samples = 100000,
    n_generated_samples = 50000,
    n_adaptive_samples = 10000,
    diag_s_elements_initial = rep(1, 2),
    m_initial = rep(0, 2),
    qbar_initial = qbar_initial,
    pbar_initial = pbar_initial, 
    random_state = NULL, 
    rtol = NULL,
    atol = NULL,
    verbose_at_refresh = T
  )
  
  test_run$sample_step_run$q_original_samples
  
}

parallel::stopCluster(init_cluster)

plot(store_samples_matrix)

mean(store_samples_matrix[, 2])

hist(store_samples_matrix[, 2], probability = T, breaks = 50, xlab = expression(q[2]), 
     main = expression(paste("Histogram of ", q[2])))
lines(seq(from = - 5, to = 7, by = 0.01), true_pdf_q2(seq(from = - 5, to = 7, by = 0.01)), col = "red")

ecdf_q2 <- ecdf(store_samples_matrix[, 2])
ecdf_q2_vec <- ecdf_q2(store_samples_matrix[, 2])
ecdf_q2_df <- data.frame(x = store_samples_matrix[, 2], y = ecdf_q2_vec)

ecdf_q2_df %>%
  ggplot(aes(x = x, y = y)) + 
  geom_point(size = 1) + 
  geom_line(stat = "function", fun = true_cdf_q2, color = "red", linetype = "dashed", linewidth = 1) + 
  xlab(bquote(q[2])) + 
  ylab("CDF") + 
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 20)
  )
