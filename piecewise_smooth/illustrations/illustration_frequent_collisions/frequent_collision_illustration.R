rm(list = ls())
library(doParallel)
library(dplyr)
source("piecewise_smooth/implementation_scripts/general_scripts/grhmc_piecewise_smooth_density_transformed_function.R")

c1 <- exp(-1 / 2)
c2 <- exp(-1 / 8) # both are normalising constant in order to get 1 when integrating over the whole R^2.  


marginal_q1_pdf <- function(q1) {
  
  ifelse(
    abs(q1) > 1,
    c1 / c2 * dnorm(q1, mean = 0, sd = 2),
    c1 / c2 * dnorm(q1, mean = 0, sd = 2) * 2 * (1 - pnorm(sqrt(1 - q1 ^ 2) / 2)) + 
      dnorm(q1) * (2 * pnorm(sqrt(1 - q1 ^ 2)) - 1)
  )
  
}

# Deterministic

deterministic_model_list <- list(
  
  n_parameters = 2, 
  
  target_jump_fun = function(q) {
    
    region_id <<- as.integer(sum(q ^ 2) <= 1) + 1
    
  },
  
  region_lin_root_list = NULL,
  
  additional_lin_root_list = NULL, 
  
  log_target_grad_list = list(
    
    list(
      function(q) -diag(rep(0.25, 2)) %*% q,
      function(q) -q
    )
    
  ),
  
  target_list = list(
    
    list(
      function(q) mvtnorm::dmvnorm(q, sigma = diag(rep(4, 2))) * c1 / c2,
      function(q) mvtnorm::dmvnorm(q)
    )
    
  ),
  
  log_target_list = list(
    
    list(
      function(q) log(mvtnorm::dmvnorm(q, sigma = diag(rep(4, 2))) * c1 / c2),
      function(q) mvtnorm::dmvnorm(q, log = T)
    )
    
  ),
  
  sim_q0 = function() rnorm(2)
  
)

deterministic_additional_non_lin_root_fun <- function(t, y, parms, m_vector, s_vector, adaptive_yes_or_no) {
  
  # print(paste0("m_vector:", m_vector))
  # print(paste0("s_vector:", s_vector))
  d <- deterministic_model_list$n_parameters
  # print(paste0("qbar: ", y[4:(4 + d - 1)]))
  # print(paste0("pbar: ", y[(4 + 2 * d):(4 + 2 * d - 1)]))
  
  q <- m_vector + s_vector * y[4:(4 + d - 1)]
  
  -sum(q ^ 2) + 1
  
}

adaptive_num_of_y_elements <- 4 + 7 * deterministic_model_list$n_parameters - 1
sampling_num_of_y_elements <- 4 + 6 * deterministic_model_list$n_parameters - 1


deterministic_additional_non_lin_event_fun <- function(t, y, parms, m_vector, s_vector, adaptive_yes_or_no) {
  
  d <- deterministic_model_list$n_parameters
  
  # print("#######################")
  # print("Boundary crossing")
  # print(paste0("Current region: ", region_id))
  # print(paste0("t: ", t))
  
  old_qbar <- y[4:(4 + d - 1)]
  old_pbar <- y[(4 + d):(4 + 2 * d - 1)]
  
  # print(paste0("Current qbar: ", old_qbar))
  # print(paste0("Current pbar: ", old_pbar))
  
  old_q <- m_vector + s_vector * old_qbar
  old_q[abs(old_q) <= 1e-13] <- 0
  # print(paste0("Current q: ", old_q))
  # print(paste0("Current q length squared: ", sum(old_q^2)))
  
  normal_vec <- c(-2 * old_q[1], -2 * old_q[2]) * s_vector
  # print(paste0("normal vec: ", normal_vec))
  
  # Next: Find the projection of pbar onto the normal vec
  old_pbar_perpendicular <- sum(old_pbar * normal_vec) / sum(normal_vec ^ 2) * normal_vec
  # print(paste0("pbar_perpendicular: ", old_pbar_perpendicular))
  old_pbar_parallel <- old_pbar - old_pbar_perpendicular
  # print(paste0("pbar_parallel: ", old_pbar_parallel))
  
  norm_squared_old_pbar_perpendicular <- sum(old_pbar_perpendicular ^ 2)
  # print(paste0("norm_squared_old_pbar_perpendicular: ", norm_squared_old_pbar_perpendicular))
  
  # For this restriction, the method of adding a tiny amount of the momentum does not always work as the boundary is curved
  # To ensure that we evaluate the potential energy of a new region, a small amount of the vector from the origin to the current qbar itself is added/subtracted
  
  if (region_id == 2) {
    new_qbar_eps <- old_qbar - normal_vec * 1e-5
    old_qbar_eps <- old_qbar + normal_vec * 1e-5
  } else {
    old_qbar_eps <- old_qbar - normal_vec * 1e-5
    new_qbar_eps <- old_qbar + normal_vec * 1e-5
  }
  
  new_q_eps <- m_vector + s_vector * new_qbar_eps
  # print(paste0("New q: ", new_q_eps))
  # print(paste0("New q length squared: ", sum(new_q_eps^2)))
  
  old_q_eps <- m_vector + s_vector * old_qbar_eps
  # print(paste0("Old q: ", old_q_eps))
  # print(paste0("Old q length squared: ", sum(old_q_eps^2)))
  
  old_potential_energy <- -deterministic_model_list$log_target_fun(old_q)
  # print(paste0("old_potential_energy: ", old_potential_energy))
  
  deterministic_model_list$target_jump_fun(new_q_eps)
  
  new_potential_energy <- -deterministic_model_list$log_target_fun(old_q) # log det(S) cancel out? 
  # print(paste0("new_potential_energy: ", new_potential_energy))
  
  delta_U <- new_potential_energy - old_potential_energy
  
  # print(paste0("delta_U: ", delta_U))
  
  if (delta_U == 0) {
    stop("Numerical issue - tweak either offset of root finder or increment")
  }
  
  if (norm_squared_old_pbar_perpendicular > 2 * delta_U) {
    
    # print("Cross boundary")
    
    new_pbar_perpendicular <- sqrt(norm_squared_old_pbar_perpendicular - 2 * delta_U) * old_pbar_perpendicular / sqrt(norm_squared_old_pbar_perpendicular)
    
    new_pbar <- old_pbar_parallel + new_pbar_perpendicular
    # print(paste0("new pbar:", new_pbar))
    
    if (adaptive_yes_or_no) {
      
      new.state <- c(
        y[1:3],
        old_qbar,
        # new_qbar_eps,
        new_pbar,
        y[(4 + 2 * d):adaptive_num_of_y_elements]
      )
      
    } else {
      
      if (sampling_compute_temporal_averages_of_moments) {
        new.state <- c(
          y[1:3],
          old_qbar,
          # new_qbar_eps,
          new_pbar,
          y[(4 + 2 * d):sampling_num_of_y_elements]
        )
      } else {
        new.state <- c(
          y[1:3],
          old_qbar,
          # new_qbar_eps,
          new_pbar
        )
      }
      
      
      
    }
    
    # print(paste0("New region: ", region_id))
    
  } else {
    # print("No cross boundary")
    new_pbar_perpendicular <- -old_pbar_perpendicular
    
    new_pbar <- old_pbar_parallel + new_pbar_perpendicular
    new_qbar <- old_qbar_eps
    # print(paste0("new pbar:", new_pbar))
    
    deterministic_model_list$target_jump_fun(old_q_eps) # ensure that the region just before hit is relevant again
    # print(paste0("New region: ", region_id))
    
    deterministic_num_of_reflections <<- deterministic_num_of_reflections + 1
    deterministic_store_q_collisions[deterministic_num_of_reflections, ] <<- old_qbar
    
    if (adaptive_yes_or_no) {
      
      if (!nut_time_found_indicator) {
        
        sum_all_nut_times <<- sum_all_nut_times + (t - t_at_current_qbar)
        
      } else {
        
        nut_time_found_indicator <<- FALSE
        
      }
      
      t_at_current_qbar <<- t
      current_qbar <<- old_qbar
      
      if (t > proportion_time_until_adaptive_start * time_period_adaptive) {
        
        lambda_adaptive <<- n_uncensored_nut_times / sum_all_nut_times
        
        if (lambda_adaptive < lambda_lower_limit) {
          lambda_adaptive <<- lambda_lower_limit
        }
        
      }
      
      new.state <- c(
        y[1:3],
        old_qbar,
        # new_qbar,
        new_pbar,
        y[(4 + 2 * d):(4 + 7 * d - 1)]
      )
    } else {
      
      if (sampling_compute_temporal_averages_of_moments) {
        new.state <- c(
          y[1:3],
          old_qbar,
          # new_qbar,
          new_pbar,
          y[(4 + 2 * d):(4 + 6 * d - 1)]
        )
      } else {
        new.state <- c(
          y[1:3],
          old_qbar,
          # new_qbar,
          new_pbar
        )
      }
      
      
    }
    
  }
  
  new.state
  
}

deterministic_model_list$additional_non_lin_root_list = list(root_fun = deterministic_additional_non_lin_root_fun, event_fun = deterministic_additional_non_lin_event_fun)


deterministic_model_list$log_target_grad <- function(q) {
  
  deterministic_model_list$log_target_grad_list[[1]][[region_id]](q)
  
}

deterministic_model_list$target_fun <- function(q) {
  
  deterministic_model_list$target_list[[1]][[region_id]](q)
  
}

deterministic_model_list$log_target_fun <- function(q) {
  
  deterministic_model_list$log_target_list[[1]][[region_id]](q)
  
}

set.seed(1)
qbar_initial <- c(0.95, 0.02)
pbar_initial <- c(0.05, 2)
u_initial <- rexp(1)
deterministic_num_of_reflections <- 0
deterministic_store_q_collisions <- matrix(nrow = 100000, ncol = 2)

system.time(
  deterministic_test_run <- grhmc_piecewise_smooth_density_transformed_function(
    model_list = deterministic_model_list,
    lambda = 0.001,
    # T = 100,
    T = 5,
    n_samples = 10000,
    diag_s_elements_initial = rep(1, 2),
    m_initial = rep(0, 2),
    qbar_initial = qbar_initial,
    pbar_initial = pbar_initial,
    Lambda_initial = 0,
    u_initial = u_initial,
    random_state = NULL,
    rtol = NULL,
    atol = NULL,
    h_max = 1.0,
    reflection_type = NULL, 
    verbose_at_refresh = T,
    last.root.offset.lin.root.finder = 1e-10,
    last.root.offset.non.lin.root.finder = 1e-5,
    precision_real_root_lin_root_finder = 1.0e-13,
    num_subdiv_non_lin_root_finder = 20L
  )
)

# Randomized 

randomized_model_list <- list(
  
  n_parameters = 2, 
  
  target_jump_fun = function(q) {
    
    region_id <<- as.integer(sum(q ^ 2) <= 1) + 1
    
  },
  
  region_lin_root_list = NULL,
  
  additional_lin_root_list = NULL, 
  
  log_target_grad_list = list(
    
    list(
      function(q) -diag(rep(0.25, 2)) %*% q,
      function(q) -q
    )
    
  ),
  
  target_list = list(
    
    list(
      function(q) mvtnorm::dmvnorm(q, sigma = diag(rep(4, 2))) * c1 / c2,
      function(q) mvtnorm::dmvnorm(q)
    )
    
  ),
  
  log_target_list = list(
    
    list(
      function(q) log(mvtnorm::dmvnorm(q, sigma = diag(rep(4, 2))) * c1 / c2),
      function(q) mvtnorm::dmvnorm(q, log = T)
    )
    
  ),
  
  sim_q0 = function() rnorm(2)
  
)

randomized_additional_non_lin_root_fun <- function(t, y, parms, m_vector, s_vector, adaptive_yes_or_no) {
  
  # print(paste0("m_vector:", m_vector))
  # print(paste0("s_vector:", s_vector))
  d <- randomized_model_list$n_parameters
  # print(paste0("qbar: ", y[4:(4 + d - 1)]))
  # print(paste0("pbar: ", y[(4 + 2 * d):(4 + 2 * d - 1)]))
  
  q <- m_vector + s_vector * y[4:(4 + d - 1)]
  
  -sum(q ^ 2) + 1
  
}

randomized_adaptive_num_of_y_elements <- 4 + 7 * randomized_model_list$n_parameters - 1
randomized_sampling_num_of_y_elements <- 4 + 6 * randomized_model_list$n_parameters - 1


randomized_additional_non_lin_event_fun <- function(t, y, parms, m_vector, s_vector, adaptive_yes_or_no) {
  
  d <- randomized_model_list$n_parameters
  
  # print("#######################")
  # print("Boundary crossing")
  # print(paste0("Current region: ", region_id))
  # print(paste0("t: ", t))
  
  old_qbar <- y[4:(4 + d - 1)]
  old_pbar <- y[(4 + d):(4 + 2 * d - 1)]
  
  # print(paste0("Current qbar: ", old_qbar))
  # print(paste0("Current pbar: ", old_pbar))
  
  old_q <- m_vector + s_vector * old_qbar
  old_q[abs(old_q) <= 1e-13] <- 0
  # print(paste0("Current q: ", old_q))
  # print(paste0("Current q length squared: ", sum(old_q^2)))
  
  normal_vec <- c(-2 * old_q[1], -2 * old_q[2]) * s_vector
  # print(paste0("normal vec: ", normal_vec))
  
  # Next: Find the projection of pbar onto the normal vec
  old_pbar_perpendicular <- sum(old_pbar * normal_vec) / sum(normal_vec ^ 2) * normal_vec
  # print(paste0("pbar_perpendicular: ", old_pbar_perpendicular))
  old_pbar_parallel <- old_pbar - old_pbar_perpendicular
  # print(paste0("pbar_parallel: ", old_pbar_parallel))
  
  norm_squared_old_pbar_perpendicular <- sum(old_pbar_perpendicular ^ 2)
  # print(paste0("norm_squared_old_pbar_perpendicular: ", norm_squared_old_pbar_perpendicular))
  
  # For this restriction, the method of adding a tiny amount of the momentum does not always work as the boundary is curved
  # To ensure that we evaluate the potential energy of a new region, a small amount of the vector from the origin to the current qbar itself is added/subtracted
  
  if (region_id == 2) {
    new_qbar_eps <- old_qbar - normal_vec * 1e-5
    old_qbar_eps <- old_qbar + normal_vec * 1e-5
  } else {
    old_qbar_eps <- old_qbar - normal_vec * 1e-5
    new_qbar_eps <- old_qbar + normal_vec * 1e-5
  }
  
  new_q_eps <- m_vector + s_vector * new_qbar_eps
  # print(paste0("New q: ", new_q_eps))
  # print(paste0("New q length squared: ", sum(new_q_eps^2)))
  
  old_q_eps <- m_vector + s_vector * old_qbar_eps
  # print(paste0("Old q: ", old_q_eps))
  # print(paste0("Old q length squared: ", sum(old_q_eps^2)))
  
  old_potential_energy <- -randomized_model_list$log_target_fun(old_q)
  # print(paste0("old_potential_energy: ", old_potential_energy))
  
  randomized_model_list$target_jump_fun(new_q_eps)
  
  new_potential_energy <- -randomized_model_list$log_target_fun(old_q) # log det(S) cancel out? 
  # print(paste0("new_potential_energy: ", new_potential_energy))
  
  delta_U <- new_potential_energy - old_potential_energy
  
  # print(paste0("delta_U: ", delta_U))
  
  if (delta_U == 0) {
    stop("Numerical issue - tweak either offset of root finder or increment")
  }
  
  if (norm_squared_old_pbar_perpendicular > 2 * delta_U) {
    
    # print("Cross boundary")
    
    new_pbar_perpendicular <- sqrt(norm_squared_old_pbar_perpendicular - 2 * delta_U) * old_pbar_perpendicular / sqrt(norm_squared_old_pbar_perpendicular)
    
    new_pbar <- old_pbar_parallel + new_pbar_perpendicular
    # print(paste0("new pbar:", new_pbar))
    
    if (adaptive_yes_or_no) {
      
      new.state <- c(
        y[1:3],
        old_qbar,
        # new_qbar_eps,
        new_pbar,
        y[(4 + 2 * d):adaptive_num_of_y_elements]
      )
      
    } else {
      
      if (sampling_compute_temporal_averages_of_moments) {
        new.state <- c(
          y[1:3],
          old_qbar,
          # new_qbar_eps,
          new_pbar,
          y[(4 + 2 * d):sampling_num_of_y_elements]
        )
      } else {
        new.state <- c(
          y[1:3],
          old_qbar,
          # new_qbar_eps,
          new_pbar
        )
      }
      
      
      
    }
    
    # print(paste0("New region: ", region_id))
    
  } else {
    # print("No cross boundary")
    z <- rnorm(d)
    
    new_pbar <- z - sum((old_pbar + z) * normal_vec) / sum(normal_vec ^ 2) * normal_vec
    new_qbar <- old_qbar_eps
    # print(paste0("new pbar:", new_pbar))
    
    randomized_model_list$target_jump_fun(old_q_eps) # ensure that the region just before hit is relevant again
    # print(paste0("New region: ", region_id))
    
    randomized_num_of_reflections <<- randomized_num_of_reflections + 1
    randomized_store_q_collisions[randomized_num_of_reflections, ] <<- old_qbar
    
    if (adaptive_yes_or_no) {
      
      if (!nut_time_found_indicator) {
        
        sum_all_nut_times <<- sum_all_nut_times + (t - t_at_current_qbar)
        
      } else {
        
        nut_time_found_indicator <<- FALSE
        
      }
      
      t_at_current_qbar <<- t
      current_qbar <<- old_qbar
      
      if (t > proportion_time_until_adaptive_start * time_period_adaptive) {
        
        lambda_adaptive <<- n_uncensored_nut_times / sum_all_nut_times
        
        if (lambda_adaptive < lambda_lower_limit) {
          lambda_adaptive <<- lambda_lower_limit
        }
        
      }
      
      new.state <- c(
        y[1:3],
        old_qbar,
        # new_qbar,
        new_pbar,
        y[(4 + 2 * d):(4 + 7 * d - 1)]
      )
    } else {
      if (sampling_compute_temporal_averages_of_moments) {
        new.state <- c(
          y[1:3],
          old_qbar,
          # new_qbar,
          new_pbar,
          y[(4 + 2 * d):(4 + 6 * d - 1)]
        )
      } else {
        new.state <- c(
          y[1:3],
          old_qbar,
          # new_qbar,
          new_pbar
        )
      }
      
    }
    
  }
  
  new.state
  
}

randomized_model_list$additional_non_lin_root_list = list(root_fun = randomized_additional_non_lin_root_fun, event_fun = randomized_additional_non_lin_event_fun)


randomized_model_list$log_target_grad <- function(q) {
  
  randomized_model_list$log_target_grad_list[[1]][[region_id]](q)
  
}

randomized_model_list$target_fun <- function(q) {
  
  randomized_model_list$target_list[[1]][[region_id]](q)
  
}

randomized_model_list$log_target_fun <- function(q) {
  
  randomized_model_list$log_target_list[[1]][[region_id]](q)
  
}

par(pty = "s")
plot(deterministic_test_run$q_original_samples, type = "l", lwd = 0.1, ylab = expression(q[2]), xlab = expression(q[1]))
lines(cos(seq(from = 0, to = 2 * pi, by = 0.01)), sin(seq(from = 0, to = 2 * pi, by = 0.01)), col = "red", lwd = 2)
points(deterministic_store_q_collisions, col = "green", pch = 19)
for (i in 1:100) {
  
  set.seed(1)
  qbar_initial <- c(0.95, 0.02)
  pbar_initial <- c(0.05, 2)
  u_initial <- rexp(1)
  randomized_num_of_reflections <- 0
  randomized_store_q_collisions <- matrix(nrow = 100000, ncol = 2)
  
  system.time(
    randomized_test_run <- grhmc_piecewise_smooth_density_transformed_function(
      model_list = randomized_model_list,
      lambda = 0.001,
      T = 1,
      n_samples = 10000,
      diag_s_elements_initial = rep(1, 2),
      m_initial = rep(0, 2),
      qbar_initial = qbar_initial,
      pbar_initial = pbar_initial,
      Lambda_initial = 0,
      u_initial = u_initial,
      random_state = i,
      rtol = NULL,
      atol = NULL,
      h_max = 1.0,
      reflection_type = NULL, 
      verbose_at_refresh = T,
      last.root.offset.lin.root.finder = 1e-10,
      last.root.offset.non.lin.root.finder = 1e-5,
      precision_real_root_lin_root_finder = 1.0e-13,
      num_subdiv_non_lin_root_finder = 20L
    )
  )
  lines(randomized_test_run$q_original_samples, col = "blue", lwd = 1, lty = 2)
}
arrows(qbar_initial[1], qbar_initial[2], qbar_initial[1] + 0.1 * pbar_initial[1], qbar_initial[2] + 0.1 * pbar_initial[2], col = "pink", lwd = 5)
points(matrix(qbar_initial, nrow = 1), col = "orange", pch = 19, lwd = 5)

