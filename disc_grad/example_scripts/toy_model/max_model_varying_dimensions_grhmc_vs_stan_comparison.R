rm(list = ls())

library(doParallel)
library(dplyr)
library(ggplot2)
source("disc_grad/implementation_scripts/ISG/full_ISG_grhmc_discontinuous_gradient_transformed_function.R")

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

dim_vec <- 2 ^ {1:11}
n_chains <- 100
n_samples <- 10000

# Stan

stan_results <- readRDS("disc_grad/example_scripts/toy_model/max_model_varying_dimensions_stan_results.RDS")
stan_final_result_matrix <- stan_results$final_result_matrix
stan_final_ess_matrix <- stan_results$final_ess_matrix
stan_diff_result_matrix <- stan_results$diff_result_matrix
stan_final_individual_ess_matrix <- stan_results$final_individual_ess_matrix
colnames(stan_final_individual_ess_matrix)[4] <- "dim"

# Plot of difference between estimates of E(q_2) and true E(q_2) 
par(mfrow = c(1, 3), oma = c(0, 0, 2, 0))

plot(
  # exp(jitter(log2(diff_result_matrix$dim), amount = 0.1)), 
  (jitter(log2(stan_diff_result_matrix$dim), amount = 0.25)),
  stan_diff_result_matrix$q2, 
  # log = "x",
  xaxt = "n",
  xlab = "Dimension",
  ylab = "Deviation from true value",
  main = expression(E(q[2]))
)

axis(1, at = log2(dim_vec), labels = dim_vec)
points(
  log2(dim_vec), 
  sapply(1:length(dim_vec), function(i) mean(stan_diff_result_matrix$q2[stan_diff_result_matrix$dim == dim_vec[i]])),
  col = "red",
  pch = 15
)
segments(
  log2(dim_vec), 
  sapply(1:length(dim_vec), function(i) mean(stan_diff_result_matrix$q2[stan_diff_result_matrix$dim == dim_vec[i]]) - 
           sd(stan_diff_result_matrix$q2[stan_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  log2(dim_vec), 
  sapply(1:length(dim_vec), function(i) mean(stan_diff_result_matrix$q2[stan_diff_result_matrix$dim == dim_vec[i]]) + 
           sd(stan_diff_result_matrix$q2[stan_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  col = "red"
)
cap <- 0.25
segments(
  log2(dim_vec) - 0.1, 
  sapply(1:length(dim_vec), function(i) mean(stan_diff_result_matrix$q2[stan_diff_result_matrix$dim == dim_vec[i]]) - 
           sd(stan_diff_result_matrix$q2[stan_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  log2(dim_vec) + 0.1, 
  sapply(1:length(dim_vec), function(i) mean(stan_diff_result_matrix$q2[stan_diff_result_matrix$dim == dim_vec[i]]) - 
           sd(stan_diff_result_matrix$q2[stan_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  col = "red"
)
segments(
  log2(dim_vec) - 0.1, 
  sapply(1:length(dim_vec), function(i) mean(stan_diff_result_matrix$q2[stan_diff_result_matrix$dim == dim_vec[i]]) + 
           sd(stan_diff_result_matrix$q2[stan_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  log2(dim_vec) + 0.1, 
  sapply(1:length(dim_vec), function(i) mean(stan_diff_result_matrix$q2[stan_diff_result_matrix$dim == dim_vec[i]]) + 
           sd(stan_diff_result_matrix$q2[stan_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  col = "red"
)
abline(h = 0, lty = 2, col = "blue")

# Plot of difference between estimates of E(q_2 ^ 2) and true E(q_2 ^ 2) 

plot(
  # exp(jitter(log2(stan_diff_result_matrix$dim), amount = 0.1)), 
  (jitter(log2(stan_diff_result_matrix$dim), amount = 0.25)),
  stan_diff_result_matrix$q2_squared, 
  # log = "x",
  xaxt = "n",
  xlab = "Dimension",
  ylab = "Deviation from true value",
  main = expression(E(q[2]^2))
)

axis(1, at = log2(dim_vec), labels = dim_vec)
points(
  log2(dim_vec), 
  sapply(1:length(dim_vec), function(i) mean(stan_diff_result_matrix$q2_squared[stan_diff_result_matrix$dim == dim_vec[i]])),
  col = "red",
  pch = 15
)
segments(
  log2(dim_vec), 
  sapply(1:length(dim_vec), function(i) mean(stan_diff_result_matrix$q2_squared[stan_diff_result_matrix$dim == dim_vec[i]]) - 
           sd(stan_diff_result_matrix$q2_squared[stan_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  log2(dim_vec), 
  sapply(1:length(dim_vec), function(i) mean(stan_diff_result_matrix$q2_squared[stan_diff_result_matrix$dim == dim_vec[i]]) + 
           sd(stan_diff_result_matrix$q2_squared[stan_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  col = "red"
)
cap <- 0.25
segments(
  log2(dim_vec) - 0.1, 
  sapply(1:length(dim_vec), function(i) mean(stan_diff_result_matrix$q2_squared[stan_diff_result_matrix$dim == dim_vec[i]]) - 
           sd(stan_diff_result_matrix$q2_squared[stan_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  log2(dim_vec) + 0.1, 
  sapply(1:length(dim_vec), function(i) mean(stan_diff_result_matrix$q2_squared[stan_diff_result_matrix$dim == dim_vec[i]]) - 
           sd(stan_diff_result_matrix$q2_squared[stan_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  col = "red"
)
segments(
  log2(dim_vec) - 0.1, 
  sapply(1:length(dim_vec), function(i) mean(stan_diff_result_matrix$q2_squared[stan_diff_result_matrix$dim == dim_vec[i]]) + 
           sd(stan_diff_result_matrix$q2_squared[stan_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  log2(dim_vec) + 0.1, 
  sapply(1:length(dim_vec), function(i) mean(stan_diff_result_matrix$q2_squared[stan_diff_result_matrix$dim == dim_vec[i]]) + 
           sd(stan_diff_result_matrix$q2_squared[stan_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  col = "red"
)
abline(h = 0, lty = 2, col = "blue")

# Plot of difference between estimates of E(q dot q) and true E(q dot q) 

plot(
  # exp(jitter(log2(stan_diff_result_matrix$dim), amount = 0.1)), 
  (jitter(log2(stan_diff_result_matrix$dim), amount = 0.25)),
  stan_diff_result_matrix$q1_squared_plus_q2_squared, 
  # log = "x",
  xaxt = "n",
  xlab = "Dimension",
  ylab = "Deviation from true value",
  main = expression(E(q[1]^2 + q[2]^2))
)

axis(1, at = log2(dim_vec), labels = dim_vec)
points(
  log2(dim_vec), 
  sapply(1:length(dim_vec), function(i) mean(stan_diff_result_matrix$q1_squared_plus_q2_squared[stan_diff_result_matrix$dim == dim_vec[i]])),
  col = "red",
  pch = 15
)
segments(
  log2(dim_vec), 
  sapply(1:length(dim_vec), function(i) mean(stan_diff_result_matrix$q1_squared_plus_q2_squared[stan_diff_result_matrix$dim == dim_vec[i]]) - 
           sd(stan_diff_result_matrix$q1_squared_plus_q2_squared[stan_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  log2(dim_vec), 
  sapply(1:length(dim_vec), function(i) mean(stan_diff_result_matrix$q1_squared_plus_q2_squared[stan_diff_result_matrix$dim == dim_vec[i]]) + 
           sd(stan_diff_result_matrix$q1_squared_plus_q2_squared[stan_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  col = "red"
)
cap <- 0.25
segments(
  log2(dim_vec) - 0.1, 
  sapply(1:length(dim_vec), function(i) mean(stan_diff_result_matrix$q1_squared_plus_q2_squared[stan_diff_result_matrix$dim == dim_vec[i]]) - 
           sd(stan_diff_result_matrix$q1_squared_plus_q2_squared[stan_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  log2(dim_vec) + 0.1, 
  sapply(1:length(dim_vec), function(i) mean(stan_diff_result_matrix$q1_squared_plus_q2_squared[stan_diff_result_matrix$dim == dim_vec[i]]) - 
           sd(stan_diff_result_matrix$q1_squared_plus_q2_squared[stan_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  col = "red"
)
segments(
  log2(dim_vec) - 0.1, 
  sapply(1:length(dim_vec), function(i) mean(stan_diff_result_matrix$q1_squared_plus_q2_squared[stan_diff_result_matrix$dim == dim_vec[i]]) + 
           sd(stan_diff_result_matrix$q1_squared_plus_q2_squared[stan_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  log2(dim_vec) + 0.1, 
  sapply(1:length(dim_vec), function(i) mean(stan_diff_result_matrix$q1_squared_plus_q2_squared[stan_diff_result_matrix$dim == dim_vec[i]]) + 
           sd(stan_diff_result_matrix$q1_squared_plus_q2_squared[stan_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  col = "red"
)
abline(h = 0, lty = 2, col = "blue")
mtext("Stan", outer = TRUE, cex = 1, font = 1, line = 0.5)

par(mfrow = c(1, 1))

# GRHMC

grhmc_results <- readRDS("disc_grad/example_scripts/toy_model/max_model_varying_dimensions_grhmc_results.RDS")
grhmc_final_result_matrix <- grhmc_results$final_result_matrix
grhmc_final_ess_matrix <- grhmc_results$final_ess_matrix
grhmc_diff_result_matrix <- grhmc_results$diff_result_matrix
grhmc_final_individual_ess_matrix <- grhmc_results$final_individual_ess_matrix
colnames(grhmc_final_individual_ess_matrix)[4] <- "dim"

# Plot of difference between estimates of E(q_2) and true E(q_2) 
par(mfrow = c(1, 3), oma = c(0, 0, 2, 0))

plot(
  # exp(jitter(log2(grhmc_diff_result_matrix$dim), amount = 0.1)), 
  (jitter(log2(grhmc_diff_result_matrix$dim), amount = 0.25)),
  grhmc_diff_result_matrix$q2, 
  # log = "x",
  xaxt = "n",
  xlab = "Dimension",
  ylab = "Deviation from true value",
  main = expression(E(q[2]))
)

axis(1, at = log2(dim_vec), labels = dim_vec)
points(
  log2(dim_vec), 
  sapply(1:length(dim_vec), function(i) mean(grhmc_diff_result_matrix$q2[grhmc_diff_result_matrix$dim == dim_vec[i]])),
  col = "red",
  pch = 15
)
segments(
  log2(dim_vec), 
  sapply(1:length(dim_vec), function(i) mean(grhmc_diff_result_matrix$q2[grhmc_diff_result_matrix$dim == dim_vec[i]]) - 
           sd(grhmc_diff_result_matrix$q2[grhmc_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  log2(dim_vec), 
  sapply(1:length(dim_vec), function(i) mean(grhmc_diff_result_matrix$q2[grhmc_diff_result_matrix$dim == dim_vec[i]]) + 
           sd(grhmc_diff_result_matrix$q2[grhmc_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  col = "red"
)
cap <- 0.25
segments(
  log2(dim_vec) - 0.1, 
  sapply(1:length(dim_vec), function(i) mean(grhmc_diff_result_matrix$q2[grhmc_diff_result_matrix$dim == dim_vec[i]]) - 
           sd(grhmc_diff_result_matrix$q2[grhmc_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  log2(dim_vec) + 0.1, 
  sapply(1:length(dim_vec), function(i) mean(grhmc_diff_result_matrix$q2[grhmc_diff_result_matrix$dim == dim_vec[i]]) - 
           sd(grhmc_diff_result_matrix$q2[grhmc_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  col = "red"
)
segments(
  log2(dim_vec) - 0.1, 
  sapply(1:length(dim_vec), function(i) mean(grhmc_diff_result_matrix$q2[grhmc_diff_result_matrix$dim == dim_vec[i]]) + 
           sd(grhmc_diff_result_matrix$q2[grhmc_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  log2(dim_vec) + 0.1, 
  sapply(1:length(dim_vec), function(i) mean(grhmc_diff_result_matrix$q2[grhmc_diff_result_matrix$dim == dim_vec[i]]) + 
           sd(grhmc_diff_result_matrix$q2[grhmc_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  col = "red"
)
abline(h = 0, lty = 2, col = "blue")

# Plot of difference between estimates of E(q_2 ^ 2) and true E(q_2 ^ 2) 

plot(
  # exp(jitter(log2(grhmc_diff_result_matrix$dim), amount = 0.1)), 
  (jitter(log2(grhmc_diff_result_matrix$dim), amount = 0.25)),
  grhmc_diff_result_matrix$q2_squared, 
  # log = "x",
  xaxt = "n",
  xlab = "Dimension",
  ylab = "Deviation from true value",
  main = expression(E(q[2]^2))
)

axis(1, at = log2(dim_vec), labels = dim_vec)
points(
  log2(dim_vec), 
  sapply(1:length(dim_vec), function(i) mean(grhmc_diff_result_matrix$q2_squared[grhmc_diff_result_matrix$dim == dim_vec[i]])),
  col = "red",
  pch = 15
)
segments(
  log2(dim_vec), 
  sapply(1:length(dim_vec), function(i) mean(grhmc_diff_result_matrix$q2_squared[grhmc_diff_result_matrix$dim == dim_vec[i]]) - 
           sd(grhmc_diff_result_matrix$q2_squared[grhmc_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  log2(dim_vec), 
  sapply(1:length(dim_vec), function(i) mean(grhmc_diff_result_matrix$q2_squared[grhmc_diff_result_matrix$dim == dim_vec[i]]) + 
           sd(grhmc_diff_result_matrix$q2_squared[grhmc_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  col = "red"
)
cap <- 0.25
segments(
  log2(dim_vec) - 0.1, 
  sapply(1:length(dim_vec), function(i) mean(grhmc_diff_result_matrix$q2_squared[grhmc_diff_result_matrix$dim == dim_vec[i]]) - 
           sd(grhmc_diff_result_matrix$q2_squared[grhmc_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  log2(dim_vec) + 0.1, 
  sapply(1:length(dim_vec), function(i) mean(grhmc_diff_result_matrix$q2_squared[grhmc_diff_result_matrix$dim == dim_vec[i]]) - 
           sd(grhmc_diff_result_matrix$q2_squared[grhmc_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  col = "red"
)
segments(
  log2(dim_vec) - 0.1, 
  sapply(1:length(dim_vec), function(i) mean(grhmc_diff_result_matrix$q2_squared[grhmc_diff_result_matrix$dim == dim_vec[i]]) + 
           sd(grhmc_diff_result_matrix$q2_squared[grhmc_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  log2(dim_vec) + 0.1, 
  sapply(1:length(dim_vec), function(i) mean(grhmc_diff_result_matrix$q2_squared[grhmc_diff_result_matrix$dim == dim_vec[i]]) + 
           sd(grhmc_diff_result_matrix$q2_squared[grhmc_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  col = "red"
)
abline(h = 0, lty = 2, col = "blue")

# Plot of difference between estimates of E(q dot q) and true E(q dot q) 

plot(
  # exp(jitter(log2(grhmc_diff_result_matrix$dim), amount = 0.1)), 
  (jitter(log2(grhmc_diff_result_matrix$dim), amount = 0.25)),
  grhmc_diff_result_matrix$q1_squared_plus_q2_squared, 
  # log = "x",
  xaxt = "n",
  xlab = "Dimension",
  ylab = "Deviation from true value",
  main = expression(E(q[1]^2 + q[2]^2))
)

axis(1, at = log2(dim_vec), labels = dim_vec)
points(
  log2(dim_vec), 
  sapply(1:length(dim_vec), function(i) mean(grhmc_diff_result_matrix$q1_squared_plus_q2_squared[grhmc_diff_result_matrix$dim == dim_vec[i]])),
  col = "red",
  pch = 15
)
segments(
  log2(dim_vec), 
  sapply(1:length(dim_vec), function(i) mean(grhmc_diff_result_matrix$q1_squared_plus_q2_squared[grhmc_diff_result_matrix$dim == dim_vec[i]]) - 
           sd(grhmc_diff_result_matrix$q1_squared_plus_q2_squared[grhmc_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  log2(dim_vec), 
  sapply(1:length(dim_vec), function(i) mean(grhmc_diff_result_matrix$q1_squared_plus_q2_squared[grhmc_diff_result_matrix$dim == dim_vec[i]]) + 
           sd(grhmc_diff_result_matrix$q1_squared_plus_q2_squared[grhmc_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  col = "red"
)
cap <- 0.25
segments(
  log2(dim_vec) - 0.1, 
  sapply(1:length(dim_vec), function(i) mean(grhmc_diff_result_matrix$q1_squared_plus_q2_squared[grhmc_diff_result_matrix$dim == dim_vec[i]]) - 
           sd(grhmc_diff_result_matrix$q1_squared_plus_q2_squared[grhmc_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  log2(dim_vec) + 0.1, 
  sapply(1:length(dim_vec), function(i) mean(grhmc_diff_result_matrix$q1_squared_plus_q2_squared[grhmc_diff_result_matrix$dim == dim_vec[i]]) - 
           sd(grhmc_diff_result_matrix$q1_squared_plus_q2_squared[grhmc_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  col = "red"
)
segments(
  log2(dim_vec) - 0.1, 
  sapply(1:length(dim_vec), function(i) mean(grhmc_diff_result_matrix$q1_squared_plus_q2_squared[grhmc_diff_result_matrix$dim == dim_vec[i]]) + 
           sd(grhmc_diff_result_matrix$q1_squared_plus_q2_squared[grhmc_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  log2(dim_vec) + 0.1, 
  sapply(1:length(dim_vec), function(i) mean(grhmc_diff_result_matrix$q1_squared_plus_q2_squared[grhmc_diff_result_matrix$dim == dim_vec[i]]) + 
           sd(grhmc_diff_result_matrix$q1_squared_plus_q2_squared[grhmc_diff_result_matrix$dim == dim_vec[i]]) * 1.96 / sqrt(n_chains)),
  col = "red"
)
abline(h = 0, lty = 2, col = "blue")
mtext("GRHMC", outer = TRUE, cex = 1, font = 1, line = 0.5)

par(mfrow = c(1, 1))

# Stan vs GRHMC - ESS

par(mfrow = c(1, 2))

# q1

plot(
  log2(dim_vec),
  log(stan_final_ess_matrix$q1 / stan_final_ess_matrix$n_evals_ode * 1e6), 
  yaxt = "n",
  xaxt = "n",
  xlab = "Dimension",
  ylab = "ESS per 1 million gradient evaluations",
  ylim = c(1, 15),
  # main = expression(E(q[1]^2 + q[2]^2)),
  main = expression(q[1])
)


axis(1, at = log2(dim_vec), labels = dim_vec)
axis(2, at = log(stan_final_ess_matrix$q1 / stan_final_ess_matrix$n_evals_ode * 1e6), labels = round(stan_final_ess_matrix$q1 / stan_final_ess_matrix$n_evals_ode * 1e6), 2)

points(
  log2(dim_vec), 
  log(grhmc_final_ess_matrix$q1 / grhmc_final_ess_matrix$n_evals_ode * 1e6),
  col = "red",
  pch = 15
)

legend(
  "bottomleft",
  legend = c("Stan", "GRHMC"),
  pch = c(1, 15),
  col = c("black", "red"),
  cex = 0.75
)

# q2 and similarly for the remaining q_j

plot(
  log2(dim_vec),
  # log(grhmc_final_ess_matrix$q2 / grhmc_final_ess_matrix$n_evals_ode * 1e6), 
  log(stan_final_ess_matrix$q2 / stan_final_ess_matrix$n_evals_ode * 1e6),
  yaxt = "n",
  xaxt = "n",
  xlab = "Dimension",
  ylab = "ESS per 1 million gradient evaluations",
  main = expression(q[2])
  # main = expression(E(q[1]^2 + q[2]^2)),
)


axis(1, at = log2(dim_vec), labels = dim_vec)
axis(2, at = log(stan_final_ess_matrix$q2 / stan_final_ess_matrix$n_evals_ode * 1e6), labels = round((stan_final_ess_matrix$q2 / stan_final_ess_matrix$n_evals_ode * 1e6), 1))

points(
  log2(dim_vec), 
  log(grhmc_final_ess_matrix$q2 / grhmc_final_ess_matrix$n_evals_ode * 1e6),
  col = "red",
  pch = 15
)

legend(
  "bottomleft",
  legend = c("Stan", "GRHMC"),
  pch = c(1, 15),
  col = c("black", "red"),
  cex = 0.75
)

par(mfrow = c(1, 1))

#######

# Stan vs GRHMC - ESS

par(mfrow = c(1, 2))

# q1

plot(
  log2(dim_vec),
  log(stan_final_ess_matrix$q1 / stan_final_ess_matrix$n_evals_ode * 1e6), 
  yaxt = "n",
  xaxt = "n",
  xlab = "Dimension",
  ylab = "ESS per 1 million gradient evaluations",
  ylim = c(2, 12),
  # main = expression(E(q[1]^2 + q[2]^2)),
  main = expression(q[1])
)


axis(1, at = log2(dim_vec), labels = dim_vec)
axis(2, at = log(10 ^ {0:5}), labels = 10 ^ {0:5})

points(
  log2(dim_vec), 
  log(grhmc_final_ess_matrix$q1 / grhmc_final_ess_matrix$n_evals_ode * 1e6),
  col = "red",
  pch = 15
)

legend(
  "bottomleft",
  legend = c("Stan", "GRHMC"),
  pch = c(1, 15),
  col = c("black", "red"),
  cex = 0.75
)

# q2 and similarly for the remaining q_j

plot(
  log2(dim_vec),
  # log(grhmc_final_ess_matrix$q2 / grhmc_final_ess_matrix$n_evals_ode * 1e6), 
  log(stan_final_ess_matrix$q2 / stan_final_ess_matrix$n_evals_ode * 1e6),
  yaxt = "n",
  xaxt = "n",
  xlab = "Dimension",
  ylab = "ESS per 1 million gradient evaluations",
  main = expression(q[2]),
  ylim = c(3, 12)
  # main = expression(E(q[1]^2 + q[2]^2)),
)


axis(1, at = log2(dim_vec), labels = dim_vec)
axis(2, at = log(10 ^ {0:5}), labels = 10 ^ {0:5})

points(
  log2(dim_vec), 
  log(grhmc_final_ess_matrix$q2 / grhmc_final_ess_matrix$n_evals_ode * 1e6),
  col = "red",
  pch = 15
)

legend(
  "bottomleft",
  legend = c("Stan", "GRHMC"),
  pch = c(1, 15),
  col = c("black", "red"),
  cex = 0.75
)

par(mfrow = c(1, 1))

# Stan vs GRHMC - Individual ESS - box plots

stan_final_individual_ess_matrix <- cbind(stan_final_individual_ess_matrix, Method = "Stan")
grhmc_final_individual_ess_matrix <- cbind(grhmc_final_individual_ess_matrix, Method = "GRHMC")
full_final_individual_ess_matrix <- rbind(stan_final_individual_ess_matrix, grhmc_final_individual_ess_matrix)
full_final_individual_ess_matrix <- data.frame(
  full_final_individual_ess_matrix,
  log_ess_per_mill_grad_eval_q1 = log(full_final_individual_ess_matrix$q1 / full_final_individual_ess_matrix$n_evals_ode * 1e6),
  log_ess_per_mill_grad_eval_q2 = log(full_final_individual_ess_matrix$q2 / full_final_individual_ess_matrix$n_evals_ode * 1e6)
)

individual_ess_plot_q1 <- ggplot(full_final_individual_ess_matrix, aes(x = as.factor(dim), y = (q1 / n_evals_ode * 1e6), fill = Method)) + 
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.5
  ) + 
  scale_y_continuous(
    trans = "log",
    breaks = c(1, 10, 100, 1000, 10000, 100000),
    labels = c("1", "10", "100", "1000", "10000", "100000")
  ) + 
  labs(
    title = expression(q[1]),
    x = "Dimension",
    y = "ESS per 1 million gradient evaluations"
  ) + 
  theme(
    plot.title = element_text(hjust = 0.5)
  )
individual_ess_plot_q1

individual_ess_plot_q2 <- ggplot(full_final_individual_ess_matrix, aes(x = as.factor(dim), y = (q2 / n_evals_ode * 1e6), fill = Method)) + 
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.5
  ) + 
  scale_y_continuous(
    trans = "log",
    breaks = c(1, 10, 100, 1000, 10000, 100000),
    labels = c("1", "10", "100", "1000", "10000", "100000")
  ) + 
  labs(
    title = expression(q[2]),
    x = "Dimension",
    y = "ESS per 1 million gradient evaluations"
  ) + 
  theme(
    plot.title = element_text(hjust = 0.5)
  )
individual_ess_plot_q2

gridExtra::grid.arrange(individual_ess_plot_q1, individual_ess_plot_q2, ncol = 2)
