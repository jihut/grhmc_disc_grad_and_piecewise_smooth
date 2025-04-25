rm(list = ls())
source("disc_grad/illustrations/illustration_trajectories/general_scripts/grhmc_discontinuous_gradient_fixed_step_size_transformed_function.R")
library(ggplot2)
library(dplyr)

# Consider q1 ~ N(0, 1), q2 ~ N(max(0, c * q1), 1) with c = 1.
# Start point: q = (-0.5, 1), p = (1.0, -0.25).
# Evolve for 1 time unit.
# During this time period, q1 will cross q1 = 0, the boundary point of discontinuous gradient.

q_initial <- c(-0.5, 1.0)
p_initial <- c(1.0, -0.25)
time_period <- 0.75
z_initial <- c(q_initial, p_initial)
h <- 0.125
c <- 10
grad_nr1 <- function(z) { # gradient when q1 < 0
  c(
    -z[1],
    -z[2]
  )
}

grad_nr2 <- function(z) { # gradient when q1 > 0
  c(
    -z[1] + (z[2] - c * z[1]) * c,
    -z[2] + c * z[1]
  )
}

## Use LSODAR with a very high precision to estimate the "true" solution at final time point ##

lsodar_grad <<- grad_nr1 # start out when q1 is less than zero

lsodar_ode_func <- function(t, z, parms) {
  
  # z[1:2] = q
  # z[3:4] = p
  
  list(
    c(
      z[3:4],
      lsodar_grad(z)
    )
  )
  
}

lsodar_event_func <- function(t, z, parms) {
  
  if(t > 0.0001){
    lsodar_grad <<- grad_nr2 # if q1 crosses zero --> change to the part of the gradient for q1 > 0. 
  }
  
  z
  
}

lsodar_root_func <- function(t, z, parms) {
  
  z[1] # q1 = 0 is the boundary
  
}

lsodar_run <- deSolve::lsodar(
  y = z_initial,
  times = seq(from = 0, to = time_period, by = 0.01),
  func = lsodar_ode_func, 
  events = list(func = lsodar_event_func, root = TRUE),
  rootfunc = lsodar_root_func,
  rtol = 1e-12,
  atol = 1e-12
)

z_true <- lsodar_run[nrow(lsodar_run), 2:5]
z_true

## GRHMC - No truncation when hitting at boundary ##

grhmc_model_nr1 <- list(
  
  n_parameters = 2, 
  
  grad_jump_fun = function(q) {
    NULL
  },
  
  additional_non_lin_root_list = NULL, # usually: list(root_fun = ..., event_fun = ...)
  
  additional_lin_root_list = NULL, # usually: list(A = ..., B = ..., event_fun = ...)
  
  region_lin_root_list = NULL,
  
  log_target_grad = function(q) {
    
    if (q[1] > 0) {
      c(
        -q[1] + (q[2] - c * q[1]) * c,
        -q[2] + c * q[1]
      )
    } else {
      -q
    }
    
  }
  
)

grad_indices <- grhmc_model_nr1$grad_jump_fun(q_initial)

grhmc_run_nr1 <- grhmc_discontinuous_gradient_fixed_step_size_transformed_function(
  model_list = grhmc_model_nr1,
  T = time_period,
  n_samples = 1000,
  diag_s_elements_initial = rep(1, 2),
  m_initial = rep(0, 2),
  qbar_initial = q_initial,
  pbar_initial = p_initial,
  h = h,
  last.root.offset.lin.root.finder = 1.0e-8,
  last.root.offset.non.lin.root.finder = 1.0e-8,
  precision_real_root_lin_root_finder = 1.0e-13,
  num_subdiv_non_lin_root_finder = 8L
)

grhmc_run_nr1$sim_output$z_end

## GRHMC - Gradient fixed during a single integration step (Proposed method) ##

grhmc_model_nr2 <- list(
  
  n_parameters = 2,
  
  grad_jump_fun = function(q) {
    
    grad_indices <<- as.integer(q[1] > 0)
    
  },
  
  additional_non_lin_root_list = NULL, # usually: list(root_fun = ..., event_fun = ...)
  
  additional_lin_root_list = NULL, # usually: list(A = ..., B = ..., event_fun = ...)
  
  region_lin_root_list = list(A = matrix(c(1, 0), nrow = 1), B = 0), # also need to specify the break point here again
  
  log_target_grad = function(q) { # explicitly state the gradient of both q1 and q2 in the two cases
    
    c(
      -q[1] + (q[2] - c * q[1]) * c * grad_indices,
      -q[2] + c * q[1] * grad_indices
    )
    
  }
  
)

grad_indices <- grhmc_model_nr2$grad_jump_fun(q_initial)

grhmc_run_nr2 <- grhmc_discontinuous_gradient_fixed_step_size_transformed_function(
  model_list = grhmc_model_nr2,
  T = time_period,
  n_samples = 1000,
  diag_s_elements_initial = rep(1, 2),
  m_initial = rep(0, 2),
  qbar_initial = q_initial,
  pbar_initial = p_initial,
  h = h,
  last.root.offset.lin.root.finder = 1.0e-8,
  last.root.offset.non.lin.root.finder = 1.0e-8,
  precision_real_root_lin_root_finder = 1.0e-13,
  num_subdiv_non_lin_root_finder = 8L
)

## GRHMC - Use gradient of region q1 < 0 all the time ##

grhmc_model_nr3 <- list(
  
  n_parameters = 2,
  
  grad_jump_fun = function(q) {
    
    NULL
    
  },
  
  additional_non_lin_root_list = NULL, # usually: list(root_fun = ..., event_fun = ...)
  
  additional_lin_root_list = NULL, # usually: list(A = ..., B = ..., event_fun = ...)
  
  region_lin_root_list = NULL, # also need to specify the break point here again
  
  log_target_grad = function(q) { # explicitly state the gradient of both q1 and q2 in the two cases
    
    c(
      -q[1],
      -q[2]
    )
    
  }
  
)

grad_indices <- grhmc_model_nr3$grad_jump_fun(q_initial)

grhmc_run_nr3 <- grhmc_discontinuous_gradient_fixed_step_size_transformed_function(
  model_list = grhmc_model_nr3,
  T = time_period,
  n_samples = 1000,
  diag_s_elements_initial = rep(1, 2),
  m_initial = rep(0, 2),
  qbar_initial = q_initial,
  pbar_initial = p_initial,
  h = h,
  last.root.offset.lin.root.finder = 1.0e-8,
  last.root.offset.non.lin.root.finder = 1.0e-8,
  precision_real_root_lin_root_finder = 1.0e-13,
  num_subdiv_non_lin_root_finder = 8L
)

########

grhmc_run_nr1_sim_output <- grhmc_run_nr1$sim_output$z_end
grhmc_run_nr1_sim_output <- na.omit(grhmc_run_nr1_sim_output)
grhmc_run_nr1_sim_output

grhmc_run_nr2_sim_output <- grhmc_run_nr2$sim_output$z_end
grhmc_run_nr2_sim_output <- na.omit(grhmc_run_nr2_sim_output)
grhmc_run_nr2_sim_output

grhmc_run_nr3_sim_output <- grhmc_run_nr3$sim_output$z_end
grhmc_run_nr3_sim_output <- na.omit(grhmc_run_nr3_sim_output)
grhmc_run_nr3_sim_output

# pdf("illustration_trajectories/final_trajectories_plot.pdf", width = 11.7, height = 8.3)
plot(grhmc_run_nr1$q_original_samples[grhmc_run_nr1$q_original_samples[, 1] > -0.23704524, ], type = "l", lty = 1, col = "grey", xlab = expression(q[1]), ylab = expression(q[2]), main = "Trajectories") # remove first two integrator steps for more "zoomed in" plot
abline(v = 0, col = "red")
points(grhmc_run_nr1_sim_output[, 1:2], col = "black")
lines(grhmc_run_nr2$q_original_samples, type = "l", lty = 2, col = "blue")
points(grhmc_run_nr2_sim_output[, 1:2], col = "blue", pch = 3)
lines(grhmc_run_nr3$q_original_samples[, 1:2], col = "green")
points(grhmc_run_nr3_sim_output[, 1:2], col = "green")
legend(
  "bottomleft",
  legend = c(
    "Standard GRHMC trajectory",
    "Standard GRHMC integrator steps",
    "Standard GRHMC trajectory (No switch)",
    "Standard GRHMC integrator steps (No switch)",
    "Proposed GRHMC trajectory",
    "Proposed GRHMC integrator steps",
    "Discontinuous gradient boundary"
  ),
  pch = c(NA, 1, NA, 1, NA, 3, NA),
  lty = c(1, NA, 1, NA, 2, NA, 1),
  col = c("grey", "black", "green", "green", "blue", "blue", "red")
)
# dev.off()

sqrt(sum((grhmc_run_nr1_sim_output[nrow(grhmc_run_nr1_sim_output), ] - z_true) ^ 2))
sqrt(sum((grhmc_run_nr2_sim_output[nrow(grhmc_run_nr2_sim_output), ] - z_true) ^ 2))
