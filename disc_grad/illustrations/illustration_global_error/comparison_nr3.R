# rm(list = ls())
source("disc_grad/illustrations/illustration_global_error/general_scripts/grhmc_discontinuous_gradient_fixed_step_size_transformed_function.R")
library(ggplot2)
library(dplyr)

# Consider q1 ~ N(0, 1), q2 ~ N(max(0, c * q1), 1) with c = 10.
# Start point: q = (-0.5, 1), p = (1.0, -0.25).
# Evolve for 0.75 time units (in order to have one single crossing).
# During this time period, q1 will cross q1 = 0, the boundary point of discontinuous gradient.

q_initial <- c(-0.5, 1.0)
p_initial <- c(1.0, -0.25)
time_period <- 0.75
z_initial <- c(q_initial, p_initial)
# h_vec <- 2 / (2 ^ seq(from = 4, to = 10, by = 1) * 10) # reduce step size by a half every time, starting from h = 0.0125 for this example (larger step sizes lead to instability)
h_vec <- 40 ^ {-1} / 1:25

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
  if (q[1] <= 0) {
    -grad_nr1(q)
  } else {
    -grad_nr2(q)
  }
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

for (i in 1:length(h_vec)) {
  
  leapfrog_store_error[i] <- leapfrog_output(grad_U = leapfrog_grad_U, time_period = time_period, q_initial = q_initial, p_initial = p_initial, z_true = z_true, h = h_vec[i])$l2_error  
  
}

leapfrog_store_matrix <- data.frame(method = rep("Leapfrog", length(h_vec)), h = h_vec, l2_error = leapfrog_store_error)

## GRHMC - Gradient can change during a single integration step ##

# Example of this: Only explicitly changing the gradient with respect to q1 when crossing boundary
# Set the gradient of q2 to max(0, q1)

grhmc_model_nr1 <- list(
  
  n_parameters = 2,
  
  grad_jump_fun = function(q) {
    
    grad_indices <<- c(
      which(c(q[1] <= 0, q[1] > 0))
    )
    
  },
  
  additional_non_lin_root_list = NULL, # usually: list(root_fun = ..., event_fun = ...)
  
  additional_lin_root_list = NULL, # usually: list(A = ..., B = ..., event_fun = ...)
  
  region_lin_root_list = list(A = matrix(c(1, 0), nrow = 1), B = 0), # also need to specify the break point here again
  
  log_target_grad_list = list(
    
    list(
      function(q) -q[1], # gradient w.r.t q1 if q1 < 0
      function(q) -q[1] + (q[2] - c * q[1]) * c # gradient w.r.t q1 if q1 > 0
    ),
    
    list(
      function(q) -q[2] + max(0, c * q[1]) # gradient w.r.t q2 not specified explicitly for the two cases above
    )
    
  )
  
)

grhmc_model_nr1$log_target_grad <- function(q) { # define the log target grad
  
  c(
    grhmc_model_nr1$log_target_grad_list[[1]][[grad_indices]](q),
    grhmc_model_nr1$log_target_grad_list[[2]][[1]](q)
  )
  
}

# Test

grhmc_output <- function(model_list, time_period, q_initial, p_initial, z_true, h) {
  
  grad_indices <- model_list$grad_jump_fun(q_initial)
  
  grhmc_run <- grhmc_discontinuous_gradient_fixed_step_size_transformed_function(
    model_list = model_list,
    T = time_period,
    n_samples = 100,
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
  
  q_grhmc_solution <- grhmc_run$q_final 
  p_grhmc_solution <- grhmc_run$p_final
  z_grhmc_solution <- c(q_grhmc_solution, p_grhmc_solution)
  
  list(q = q_grhmc_solution, p = p_grhmc_solution, l2_error = sqrt(sum((z_grhmc_solution - z_true) ^ 2)))
  
}

grhmc_output(model_list = grhmc_model_nr1, time_period = time_period, q_initial = q_initial, p_initial = p_initial, z_true = z_true, h = 0.01)
grhmc_output(model_list = grhmc_model_nr1, time_period = time_period, q_initial = q_initial, p_initial = p_initial, z_true = z_true, h = 0.005)

# Run for a range of values for the step size h

grhmc_nr1_store_error <- numeric(length(h_vec))

for (i in 1:length(h_vec)) {
  
  grhmc_nr1_store_error[i] <- grhmc_output(model_list = grhmc_model_nr1, time_period = time_period, q_initial = q_initial, p_initial = p_initial, z_true = z_true, h = h_vec[i])$l2_error  
  
}

grhmc_nr1_store_matrix <- data.frame(method = rep("GRHMC nr.1", length(h_vec)), h = h_vec, l2_error = grhmc_nr1_store_error)

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

grhmc_output(model_list = grhmc_model_nr2, time_period = time_period, q_initial = q_initial, p_initial = p_initial, z_true = z_true, h = 0.01)
grhmc_output(model_list = grhmc_model_nr2, time_period = time_period, q_initial = q_initial, p_initial = p_initial, z_true = z_true, h = 0.005)

# Run for a range of values for the step size h

grhmc_nr2_store_error <- numeric(length(h_vec))

for (i in 1:length(h_vec)) {
  
  grhmc_nr2_store_error[i] <- grhmc_output(model_list = grhmc_model_nr2, time_period = time_period, q_initial = q_initial, p_initial = p_initial, z_true = z_true, h = h_vec[i])$l2_error  
  
}

grhmc_nr2_store_matrix <- data.frame(method = rep("GRHMC nr.2", length(h_vec)), h = h_vec, l2_error = grhmc_nr2_store_error)

## Comparison plot ##

summary(lm(l2_error ~ I(h ^ 3), data = grhmc_nr2_store_matrix))
range(c(leapfrog_store_error, grhmc_nr1_store_error, grhmc_nr2_store_error))

h_vec_c_10 <- h_vec
leapfrog_store_error_c_10 <- leapfrog_store_error
leapfrog_store_matrix_c_10 <- leapfrog_store_matrix
grhmc_nr1_store_error_c_10 <- grhmc_nr1_store_error
grhmc_nr1_store_matrix_c_10 <- grhmc_nr1_store_matrix
grhmc_nr2_store_error_c_10 <- grhmc_nr2_store_error
grhmc_nr2_store_matrix_c_10 <- grhmc_nr2_store_matrix

plot(h_vec_c_10, leapfrog_store_error_c_10, log = "xy", ylim = c(1e-10, 1e-1), ylab = "L2 error at final state", xlab = "h (step size)", main = "c = 10")
points(h_vec_c_10, grhmc_nr1_store_error_c_10, col = "red")
points(h_vec_c_10, grhmc_nr2_store_error_c_10, col = "blue", pch = 3)
lines(h_vec_c_10, 105 * h_vec_c_10 ^ 3, lty = 2, col = "green")
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
  col = c("black", "red", "blue", "green")
)

final_store_matrix_c_10 <- rbind(
  leapfrog_store_matrix_c_10, 
  grhmc_nr1_store_matrix_c_10,
  grhmc_nr2_store_matrix_c_10
)

final_store_matrix_c_10 %>%
  ggplot(aes(x = h, y = l2_error)) + 
  geom_point(aes(colour = method, shape = method)) + 
  scale_shape_manual(values = c(1, 3, 1), guide = "none") + 
  geom_function(aes(colour = "Line"), fun = function(x) 105 * x ^ 3, linetype = 2) + 
  scale_color_manual(
    values = c("red", "blue", "black", "green"),
    labels = c("GRHMC without boundary handling", "GRHMC with boundary handling", "Leapfrog", expression("" %prop% h^3))
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
