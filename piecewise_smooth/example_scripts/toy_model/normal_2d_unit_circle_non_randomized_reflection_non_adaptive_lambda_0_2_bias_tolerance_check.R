rm(list = ls())

library(doParallel)
library(dplyr)
library(ggplot2)
source("piecewise_smooth/implementation_scripts/ISG/full_ISG_grhmc_piecewise_smooth_density_transformed_function.R")

# For a range of integrator tolerance (absolute and relative)
# Run 100 independent trajectories and look at the samples generated all together. 
# Compare marginal of q_1 to verify that everything is correct. 
# Example of piecewise smooth target
# Standard normal distribution restricted to a circle with radius 1


c1 <- exp(-1 / 2)
c2 <- exp(-1 / 8) # both are normalising constant in order to get 1 when integrating over the whole R^2. 

marginal_q1_pdf <- function(q1) {
  
  ifelse(
    abs(q1) > 1,
    c1 / c2 * dnorm(q1, mean = 0, sd = 2), # marginal density of q1 if outside the unit circle
    c1 / c2 * dnorm(q1, mean = 0, sd = 2) * 2 * (1 - pnorm(sqrt(1 - q1 ^ 2) / 2)) + # marginal density of q1 if inside the unit circle
      dnorm(q1) * (2 * pnorm(sqrt(1 - q1 ^ 2)) - 1)
  )
  
}

init_cluster <- parallel::makeCluster(10)
doParallel::registerDoParallel(init_cluster)

integrator_tol <- c(1e-3, 1e-4, 1e-5, 1e-6)
n_trajectories <- 100

list_chains_tol_scenario <- list()

for (k in 1:length(integrator_tol)) {
  doRNG::registerDoRNG(seed = 42)
  
  store_chains <- foreach::foreach(i = 1:n_trajectories) %dopar% {
    
    model_list <- list(
      
      n_parameters = 2, 
      
      target_jump_fun = function(q) {
        
        region_id <<- as.integer(sum(q ^ 2) <= 1) + 1
        
      },
      
      region_lin_root_list = NULL,
      
      additional_lin_root_list = NULL, 
      
      log_target_grad_list = list(
        
        list(
          function(q) -diag(rep(0.25, 2)) %*% q, # the gradient of log target if outside the unit circle
          function(q) -q # the gradient of log target if inside the unit circle
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
          function(q) log(mvtnorm::dmvnorm(q, sigma = diag(rep(4, 2))) * c1 / c2), # log target outside the unit circle
          function(q) mvtnorm::dmvnorm(q, log = T) # log target inside the unit circle
        )
        
      ),
      
      sim_q0 = function() rnorm(2)
      
    )
    
    additional_non_lin_root_fun <- function(t, y, parms, m_vector, s_vector, adaptive_yes_or_no) {
      
      # print(paste0("m_vector:", m_vector))
      # print(paste0("s_vector:", s_vector))
      d <- model_list$n_parameters
      # print(paste0("qbar: ", y[4:(4 + d - 1)]))
      # print(paste0("pbar: ", y[(4 + 2 * d):(4 + 2 * d - 1)]))
      
      q <- m_vector + s_vector * y[4:(4 + d - 1)]
      
      -sum(q ^ 2) + 1
      
    }
    
    adaptive_num_of_y_elements <- 4 + 7 * model_list$n_parameters - 1
    sampling_num_of_y_elements <- 4 + 6 * model_list$n_parameters - 1
    
    
    additional_non_lin_event_fun <- function(t, y, parms, m_vector, s_vector, adaptive_yes_or_no) {
      
      d <- model_list$n_parameters
      
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
      
      old_potential_energy <- -model_list$log_target_fun(old_q) # potential energy just before hitting the boundary
      # print(paste0("old_potential_energy: ", old_potential_energy))
      
      model_list$target_jump_fun(new_q_eps) # make sure that one is in the new region
      
      new_potential_energy <- -model_list$log_target_fun(old_q) # log det(S) cancel out? - potential energy right after switching region
      # print(paste0("new_potential_energy: ", new_potential_energy))
      
      delta_U <- new_potential_energy - old_potential_energy # difference between the potential energy in the two regions
      
      # print(paste0("delta_U: ", delta_U))
      
      if (delta_U == 0) {
        stop("Numerical issue - tweak either offset of root finder or increment")
      }
      
      if (norm_squared_old_pbar_perpendicular > 2 * delta_U) { # if the momentum squared along the normal direction is larger than two times the potential energy --> transition to new region
        
        # print("Cross boundary")
        
        new_pbar_perpendicular <- sqrt(norm_squared_old_pbar_perpendicular - 2 * delta_U) * old_pbar_perpendicular / sqrt(norm_squared_old_pbar_perpendicular) # refract the momentum along the normal direction after moving to a different region
        
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
        
      } else { # if not, then reflect the component of the momentum along the normal direction, either by deterministic or randomized reflection
        # print("No cross boundary")
        new_pbar_perpendicular <- -old_pbar_perpendicular
        
        new_pbar <- old_pbar_parallel + new_pbar_perpendicular
        new_qbar <- old_qbar_eps
        # print(paste0("new pbar:", new_pbar))
        
        model_list$target_jump_fun(old_q_eps) # ensure that the region just before hit is relevant again
        # print(paste0("New region: ", region_id))
        
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
    
    model_list$additional_non_lin_root_list = list(root_fun = additional_non_lin_root_fun, event_fun = additional_non_lin_event_fun)
    
    
    model_list$log_target_grad <- function(q) {
      
      model_list$log_target_grad_list[[1]][[region_id]](q)
      
    }
    
    model_list$target_fun <- function(q) {
      
      model_list$target_list[[1]][[region_id]](q)
      
    }
    
    model_list$log_target_fun <- function(q) {
      
      model_list$log_target_list[[1]][[region_id]](q)
      
    }
    
    qbar_initial <- model_list$sim_q0()
    pbar_initial <- rnorm(2)
    u_initial <- rexp(1)
    
    test_run_warmup <- grhmc_piecewise_smooth_density_transformed_function(
      model_list = model_list,
      lambda = 0.2,
      T = 10000,
      n_samples = 10000,
      diag_s_elements_initial = rep(1, 2),
      m_initial = rep(0, 2),
      qbar_initial = qbar_initial,
      pbar_initial = pbar_initial,
      Lambda_initial = 0,
      u_initial = u_initial,
      random_state = NULL,
      rtol = integrator_tol[k],
      atol = integrator_tol[k],
      h_max = 1.0,
      reflection_type = NULL, 
      verbose_at_refresh = T,
      last.root.offset.lin.root.finder = 1e-10,
      last.root.offset.non.lin.root.finder = 1e-5,
      precision_real_root_lin_root_finder = 1.0e-13,
      num_subdiv_non_lin_root_finder = 20L
    )
    
    new_pbar_initial <- rnorm(2)
    new_u_initial <- rexp(1)
    
    test_run_sampling <- grhmc_piecewise_smooth_density_transformed_function(
      model_list = model_list,
      lambda = 0.2,
      T = 10000,
      n_samples = 10000,
      diag_s_elements_initial = rep(1, 2),
      m_initial = rep(0, 2),
      qbar_initial = test_run_warmup$q_original_samples[10000, ],
      pbar_initial = new_pbar_initial,
      Lambda_initial = 0,
      u_initial = new_u_initial,
      random_state = NULL,
      rtol = integrator_tol[k],
      atol = integrator_tol[k],
      h_max = 1.0,
      reflection_type = NULL, 
      verbose_at_refresh = T,
      last.root.offset.lin.root.finder = 1e-10,
      last.root.offset.non.lin.root.finder = 1e-5,
      precision_real_root_lin_root_finder = 1.0e-13,
      num_subdiv_non_lin_root_finder = 20L
    )
    
    test_run_sampling
    
  }
  
  list_chains_tol_scenario[[k]] <- store_chains
  
}

true_mean_q1 <- 
  integrate(function(x) x * marginal_q1_pdf(x), lower = -Inf, upper = Inf, rel.tol = 1e-13, abs.tol = 1e-13)
true_mean_q1
true_mean_q1_squared <- 
  integrate(function(x) (x ^ 2) * marginal_q1_pdf(x), lower = -Inf, upper = Inf, rel.tol = 1e-13, abs.tol = 1e-13)
true_mean_q1_squared

integrand_true_mean_q_dot_q <- function(q1, q2) {
  (q1 ^ 2 + q2 ^ 2) * (
    dnorm(q1, mean = 0, sd = 1) * dnorm(q2, mean = 0, sd = 1) * (q1 ^ 2 + q2 ^ 2 <= 1) + 
      c1 / c2 * dnorm(q1, mean = 0, sd = 2) * dnorm(q2, mean = 0, sd = 2) * (q1 ^ 2 + q2 ^ 2 > 1) 
  )
}

true_mean_q_dot_q <- integrate(function(q1) {
  sapply(q1, function(q1_star) {
    integrate(function(q2) integrand_true_mean_q_dot_q(q1_star, q2), lower = -Inf, upper = Inf, rel.tol = 1e-9, abs.tol = 1e-9, subdivisions = 100000)$value
  })
}, lower = -Inf, upper = Inf, rel.tol = 1e-8, abs.tol = 1e-8, subdivisions = 100000)

true_mean_q_dot_q

estimated_quantities <- function(list, calculate_wasserstein_2 = FALSE) {
  
  samples_matrix <- list$q_original_samples
  
  mean_q1 <- mean(samples_matrix[, 1])
  mean_q1_squared <- mean(samples_matrix[, 1] ^ 2)
  mean_q_dot_q <- mean(rowSums(samples_matrix * samples_matrix))
  
  if (calculate_wasserstein_2 == TRUE) {
    wasserstein_2_integrand <- function(t) {

      # quantile_via_cdf <- uniroot(
      #   function(x) true_cdf_q2(x) - t, lower = -10, upper = 10, tol = 1e-13
      # )$root
      #
      # quantile_via_sample <- quantile(list_chains_tol_scenario[[1]][[1]]$q_original_samples[, 2], probs = t)
      #
      # (quantile_via_cdf - quantile_via_sample) ^ 2

      quantile_via_cdf <- sapply(
        t,
        function(t_star) {
          uniroot(
            function(x) true_cdf_q1(x) - t_star, lower = -10, upper = 10, tol = 1e-12
          )$root
        }
      )

      quantile_via_sample <- quantile(samples_matrix[, 2], probs = t)

      (quantile_via_cdf - quantile_via_sample) ^ 2

    }

    wasserstein_2_squared_value <- 0
    grid <- seq(from = 0, to = 1, by = 0.0001)
    
    for (i in 1:(length(grid) - 1)) {
      # print(i)
      int_run <- integrate(
        wasserstein_2_integrand,
        lower = grid[i],
        upper = grid[i + 1],
        rel.tol = 1e-9,
        abs.tol = 1e-9,
        subdivisions = 100000
      )
      wasserstein_2_squared_value <- wasserstein_2_squared_value + int_run$value
    }
    return(
      c(
        q1 = mean_q1, q1_squared = mean_q1_squared, q_dot_q = mean_q_dot_q, wasserstein_2 = sqrt(wasserstein_2_squared_value)
      ) 
    )
  } else {
    return(
      c(
        q1 = mean_q1, q1_squared = mean_q1_squared, q_dot_q = mean_q_dot_q
      ) 
    )
  }

}

for (k in 1:length(integrator_tol)) {
  relevant_list <- list_chains_tol_scenario[[k]]
  relevant_matrix <- cbind(do.call(rbind, lapply(relevant_list, estimated_quantities)), 
                           int_tol = rep(integrator_tol[k]))
  if (k == 1) {
    final_result_matrix <- relevant_matrix
  } else {
    final_result_matrix <- rbind(final_result_matrix, relevant_matrix)
  }
}

final_result_matrix <- as.data.frame(final_result_matrix)

diff_result_matrix <- 
  sweep(final_result_matrix, 2, c(true_mean_q1$value, true_mean_q1_squared$value, true_mean_q_dot_q$value, 0), "-")

# Plot of difference between estimates of E(q_2) and true E(q_2) 

plot(
  exp(jitter(log(diff_result_matrix$int_tol), amount = 0.25)), 
  diff_result_matrix$q1, 
  log = "x",
  xaxt = "n",
  xlab = "Integrator tolerance",
  ylab = "Deviation from true value",
  main = expression(E(q[1]))
)

axis(1, at = integrator_tol, labels = integrator_tol)
points(
  integrator_tol, 
  sapply(1:4, function(i) mean(diff_result_matrix$q1[diff_result_matrix$int_tol == integrator_tol[i]])),
  col = "red",
  pch = 15
)
segments(
  integrator_tol, 
  sapply(1:4, function(i) mean(diff_result_matrix$q1[diff_result_matrix$int_tol == integrator_tol[i]]) - 
           sd(diff_result_matrix$q1[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  integrator_tol, 
  sapply(1:4, function(i) mean(diff_result_matrix$q1[diff_result_matrix$int_tol == integrator_tol[i]]) + 
           sd(diff_result_matrix$q1[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  col = "red"
)
cap <- 0.25
segments(
  integrator_tol / (1 + cap), 
  sapply(1:4, function(i) mean(diff_result_matrix$q1[diff_result_matrix$int_tol == integrator_tol[i]]) - 
           sd(diff_result_matrix$q1[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  integrator_tol * (1 + cap), 
  sapply(1:4, function(i) mean(diff_result_matrix$q1[diff_result_matrix$int_tol == integrator_tol[i]]) - 
           sd(diff_result_matrix$q1[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  col = "red"
)
segments(
  integrator_tol / (1 + cap), 
  sapply(1:4, function(i) mean(diff_result_matrix$q1[diff_result_matrix$int_tol == integrator_tol[i]]) + 
           sd(diff_result_matrix$q1[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  integrator_tol * (1 + cap), 
  sapply(1:4, function(i) mean(diff_result_matrix$q1[diff_result_matrix$int_tol == integrator_tol[i]]) + 
           sd(diff_result_matrix$q1[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  col = "red"
)
abline(h = 0, lty = 2, col = "blue")

# Plot of difference between estimates of E(q_2 ^ 2) and true E(q_2 ^ 2) 

plot(
  exp(jitter(log(diff_result_matrix$int_tol), amount = 0.25)),
  diff_result_matrix$q1_squared, 
  log = "x",
  xaxt = "n",
  xlab = "Integrator tolerance",
  ylab = "Deviation from true value",
  main = expression(E(q[1]^2))
)

axis(1, at = integrator_tol, labels = integrator_tol)
points(
  integrator_tol, 
  sapply(1:4, function(i) mean(diff_result_matrix$q1_squared[diff_result_matrix$int_tol == integrator_tol[i]])),
  col = "red",
  pch = 15
)
segments(
  integrator_tol, 
  sapply(1:4, function(i) mean(diff_result_matrix$q1_squared[diff_result_matrix$int_tol == integrator_tol[i]]) - 
           sd(diff_result_matrix$q1_squared[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  integrator_tol, 
  sapply(1:4, function(i) mean(diff_result_matrix$q1_squared[diff_result_matrix$int_tol == integrator_tol[i]]) + 
           sd(diff_result_matrix$q1_squared[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  col = "red"
)
cap <- 0.25
segments(
  integrator_tol / (1 + cap), 
  sapply(1:4, function(i) mean(diff_result_matrix$q1_squared[diff_result_matrix$int_tol == integrator_tol[i]]) - 
           sd(diff_result_matrix$q1_squared[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  integrator_tol * (1 + cap), 
  sapply(1:4, function(i) mean(diff_result_matrix$q1_squared[diff_result_matrix$int_tol == integrator_tol[i]]) - 
           sd(diff_result_matrix$q1_squared[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  col = "red"
)
segments(
  integrator_tol / (1 + cap), 
  sapply(1:4, function(i) mean(diff_result_matrix$q1_squared[diff_result_matrix$int_tol == integrator_tol[i]]) + 
           sd(diff_result_matrix$q1_squared[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  integrator_tol * (1 + cap), 
  sapply(1:4, function(i) mean(diff_result_matrix$q1_squared[diff_result_matrix$int_tol == integrator_tol[i]]) + 
           sd(diff_result_matrix$q1_squared[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  col = "red"
)
abline(h = 0, lty = 2, col = "blue")

# Plot of difference between estimates of E(q dot q) and true E(q dot q) 

plot(
  exp(jitter(log(diff_result_matrix$int_tol), amount = 0.25)), 
  diff_result_matrix$q_dot_q, 
  log = "x",
  xaxt = "n",
  xlab = "Integrator tolerance",
  ylab = "Deviation from true value",
  main = expression(E(q^T * q))
)

axis(1, at = integrator_tol, labels = integrator_tol)
points(
  integrator_tol, 
  sapply(1:4, function(i) mean(diff_result_matrix$q_dot_q[diff_result_matrix$int_tol == integrator_tol[i]])),
  col = "red",
  pch = 15
)
segments(
  integrator_tol, 
  sapply(1:4, function(i) mean(diff_result_matrix$q_dot_q[diff_result_matrix$int_tol == integrator_tol[i]]) - 
           sd(diff_result_matrix$q_dot_q[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  integrator_tol, 
  sapply(1:4, function(i) mean(diff_result_matrix$q_dot_q[diff_result_matrix$int_tol == integrator_tol[i]]) + 
           sd(diff_result_matrix$q_dot_q[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  col = "red"
)
cap <- 0.25
segments(
  integrator_tol / (1 + cap), 
  sapply(1:4, function(i) mean(diff_result_matrix$q_dot_q[diff_result_matrix$int_tol == integrator_tol[i]]) - 
           sd(diff_result_matrix$q_dot_q[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  integrator_tol * (1 + cap), 
  sapply(1:4, function(i) mean(diff_result_matrix$q_dot_q[diff_result_matrix$int_tol == integrator_tol[i]]) - 
           sd(diff_result_matrix$q_dot_q[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  col = "red"
)
segments(
  integrator_tol / (1 + cap), 
  sapply(1:4, function(i) mean(diff_result_matrix$q_dot_q[diff_result_matrix$int_tol == integrator_tol[i]]) + 
           sd(diff_result_matrix$q_dot_q[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  integrator_tol * (1 + cap), 
  sapply(1:4, function(i) mean(diff_result_matrix$q_dot_q[diff_result_matrix$int_tol == integrator_tol[i]]) + 
           sd(diff_result_matrix$q_dot_q[diff_result_matrix$int_tol == integrator_tol[i]]) * 1.96 / sqrt(n_trajectories)),
  col = "red"
)
abline(h = 0, lty = 2, col = "blue")
