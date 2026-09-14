rm(list = ls())

library(MASS)
library(brms)
library(ggplot2)
library(dplyr)
library(rstan)
library(foreach)
library(gridExtra)

n_iterations <- 10
n_samples_per_iteration <- 1000000

head(Boston)

x <- as.matrix(Boston[, -ncol(Boston)]) # 13 covariates
y <- Boston$medv
n_rows <- length(y)
length_beta <- ncol(x)

scale_x <- scale(x)
scale_y <- scale(y)
mean_y <- mean(y)

disc_grad_runs <- readRDS("disc_grad/example_scripts/regularized_linear_regression/disc_grad_run_identity_with_n_evals_ode_run_example.RDS")
sticky_runs <- readRDS("disc_grad/example_scripts/regularized_linear_regression/sticky_run_identity_with_n_evals_ode_run_example.RDS")

# Discontinuous gradient setup

disc_grad_samples_array <- array(dim = c(n_samples_per_iteration, n_iterations, length_beta + 1))

for (i in 1:n_iterations) {
  
  current_samples_matrix <- matrix(nrow = n_samples_per_iteration, ncol = length_beta + 1) # last column for variance parameter
  
  for (j in 1:length_beta) {
    
    current_beta_vec <- pmax(0, disc_grad_runs[[i]]$q_original_samples[, j]) - pmax(0, disc_grad_runs[[i]]$q_original_samples[, length_beta + j])
    current_samples_matrix[, j] <- current_beta_vec
  }
  
  current_samples_matrix[, length_beta + 1] <- disc_grad_runs[[i]]$q_original_samples[, 2 * length_beta + 1]
  
  current_samples_matrix[, 1:length_beta] <- t(t(current_samples_matrix[, 1:length_beta]) / apply(x, 2, sd)) * sd(y) # scale back to original scale
  
  current_samples_matrix[, length_beta + 1] <- exp(0.5 * current_samples_matrix[, length_beta + 1]) * sd(y) # scale back to original scale and sigma
  
  disc_grad_samples_array[, i, ] <- current_samples_matrix
}

# Sticky setup

sticky_samples_array <- array(dim = c(n_samples_per_iteration, n_iterations, length_beta + 1))

for (i in 1:n_iterations) {
  
  current_samples_matrix <- sticky_runs[[i]]$q_original_samples
  
  current_samples_matrix[, 1:length_beta] <- t(t(current_samples_matrix[, 1:length_beta]) / apply(x, 2, sd)) * sd(y) # scale back to original scale
  
  current_samples_matrix[, length_beta + 1] <- exp(0.5 * current_samples_matrix[, length_beta + 1]) * sd(y) # scale back to original scale and sigma
  
  sticky_samples_array[, i, ] <- current_samples_matrix
  
}

final_disc_grad_samples <- disc_grad_samples_array
dim(final_disc_grad_samples) <- c(n_iterations * n_samples_per_iteration, length_beta + 1)
final_disc_grad_samples <- cbind(final_disc_grad_samples, sort(rep(1:n_iterations, n_samples_per_iteration)))
final_disc_grad_samples <- final_disc_grad_samples[, -14] # sigma parameter
colnames(final_disc_grad_samples) <- c(colnames(x), "iteration")

final_sticky_samples <- sticky_samples_array
dim(final_sticky_samples) <- c(n_iterations * n_samples_per_iteration, length_beta + 1)
final_sticky_samples <- cbind(final_sticky_samples, sort(rep(1:n_iterations, n_samples_per_iteration)))
final_sticky_samples <- final_sticky_samples[, -14] # sigma parameter
colnames(final_sticky_samples) <- c(colnames(x), "iteration")

# Plots

# par(mfrow = c(3, 2))
# plot(final_disc_grad_samples[, "zn"], col = final_disc_grad_samples[, "iteration"], ylim = c(-0.02, 0.10), 
#      main = expression(beta %+-% ""), ylab = expression(beta[zn]))
# plot(final_sticky_samples[, "zn"], col = final_sticky_samples[, "iteration"], ylim = c(-0.02, 0.10), 
#      main = "Sticky", ylab = expression(beta[zn]))
# plot(final_disc_grad_samples[, "chas"], col = final_disc_grad_samples[, "iteration"], ylim = c(-1, 7),
#      main = expression(beta %+-% ""), ylab = expression(beta[chas]))
# plot(final_sticky_samples[, "chas"], col = final_sticky_samples[, "iteration"], ylim = c(-1, 7),
#      main = "Sticky", ylab = expression(beta[zn]))
# plot(final_disc_grad_samples[, "age"], col = final_disc_grad_samples[, "iteration"], ylim = c(-0.05, 0.03),
#      main = expression(beta %+-% ""), ylab = expression(beta[zn]))
# plot(final_sticky_samples[, "age"], col = final_sticky_samples[, "iteration"], ylim = c(-0.05, 0.03),
#      main = "Sticky", ylab = expression(beta[zn]))

par(mfrow = c(3, 2))
plot(final_disc_grad_samples[, "zn"], col = final_disc_grad_samples[, "iteration"], ylim = c(-0.02, 0.10), 
     main = expression(beta %+-% ""), ylab = expression(beta[zn]), pch = 19, cex = 0.1, xlab = "Sample number")
plot(final_sticky_samples[, "zn"], col = final_sticky_samples[, "iteration"], ylim = c(-0.02, 0.10), 
     main = "Sticky", ylab = expression(beta[zn]), pch = 19, cex = 0.1, xlab = "Sample number")
plot(final_disc_grad_samples[, "chas"], col = final_disc_grad_samples[, "iteration"], ylim = c(-1, 7),
     main = expression(beta %+-% ""), ylab = expression(beta[chas]), pch = 19, cex = 0.1, xlab = "Sample number")
plot(final_sticky_samples[, "chas"], col = final_sticky_samples[, "iteration"], ylim = c(-1, 7),
     main = "Sticky", ylab = expression(beta[zn]), pch = 19, cex = 0.1, xlab = "Sample number")
plot(final_disc_grad_samples[, "age"], col = final_disc_grad_samples[, "iteration"], ylim = c(-0.05, 0.03),
     main = expression(beta %+-% ""), ylab = expression(beta[zn]), pch = 19, cex = 0.1, xlab = "Sample number")
plot(final_sticky_samples[, "age"], col = final_sticky_samples[, "iteration"], ylim = c(-0.05, 0.03),
     main = "Sticky", ylab = expression(beta[zn]), pch = 19, cex = 0.1, xlab = "Sample number")

# Raster ggplot and downsampling

small_final_disc_grad_samples <- 
  final_disc_grad_samples[seq(from = 100, to = 1e7, by = 100), ]
small_final_disc_grad_samples <- cbind(small_final_disc_grad_samples, sample_number = 1:nrow(small_final_disc_grad_samples))

small_final_sticky_samples <- 
  final_sticky_samples[seq(from = 100, to = 1e7, by = 100), ]
small_final_sticky_samples <- cbind(small_final_sticky_samples, sample_number = 1:nrow(small_final_sticky_samples))


disc_grad_zn <- ggplot(
  small_final_disc_grad_samples,
  aes(x = sample_number, y = zn, colour = as.factor(iteration))
) + 
  geom_point_rast(
    raster.dpi = 150,
    alpha = 0.25,
    size = 0.1
  ) + 
  labs(
    title = expression(beta %+-% ""),
    x = "Sample number",
    y = expression(beta[zn])
  ) + 
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) + 
  ylim(
    c(-0.02, 0.10)
  )
disc_grad_zn

sticky_zn <- ggplot(
  small_final_sticky_samples,
  aes(x = sample_number, y = zn, colour = as.factor(iteration))
) + 
  geom_point_rast(
    raster.dpi = 150,
    alpha = 0.25,
    size = 0.1
  ) + 
  labs(
    title = "Sticky",
    x = "Sample number",
    y = expression(beta[zn])
  ) + 
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) + 
  ylim(
    c(-0.02, 0.10)
  )
sticky_zn

disc_grad_chas <- ggplot(
  small_final_disc_grad_samples,
  aes(x = sample_number, y = chas, colour = as.factor(iteration))
) + 
  geom_point_rast(
    raster.dpi = 150,
    alpha = 0.25,
    size = 0.1
  ) + 
  labs(
    title = expression(beta %+-% ""),
    x = "Sample number",
    y = expression(beta[chas])
  ) + 
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) + 
  ylim(
    c(-1, 7)
  )
disc_grad_chas

sticky_chas <- ggplot(
  small_final_sticky_samples,
  aes(x = sample_number, y = chas, colour = as.factor(iteration))
) + 
  geom_point_rast(
    raster.dpi = 150,
    alpha = 0.25,
    size = 0.1
  ) + 
  labs(
    title = "Sticky",
    x = "Sample number",
    y = expression(beta[chas])
  ) + 
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) + 
  ylim(
    c(-1, 7)
  )
sticky_chas

disc_grad_age <- ggplot(
  small_final_disc_grad_samples,
  aes(x = sample_number, y = age, colour = as.factor(iteration))
) + 
  geom_point_rast(
    raster.dpi = 150,
    alpha = 0.25,
    size = 0.1
  ) + 
  labs(
    title = expression(beta %+-% ""),
    x = "Sample number",
    y = expression(beta[age])
  ) + 
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) + 
  ylim(
    c(-0.05, 0.03)
  )
disc_grad_age

sticky_age <- ggplot(
  small_final_sticky_samples,
  aes(x = sample_number, y = age, colour = as.factor(iteration))
) + 
  geom_point_rast(
    raster.dpi = 150,
    alpha = 0.25,
    size = 0.1
  ) + 
  labs(
    title = "Sticky",
    x = "Sample number",
    y = expression(beta[age])
  ) + 
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) + 
  ylim(
    c(-0.05, 0.03)
  )
sticky_age

ggsave(
  "disc_grad/example_scripts/regularized_linear_regression/sticky_vs_disc_grad_traceplots.png",
  patchwork::wrap_plots(
    list(
      disc_grad_zn, 
      sticky_zn, 
      disc_grad_chas,
      sticky_chas,
      disc_grad_age,
      sticky_age
    ), 
    ncol = 2
  ) + 
    patchwork::plot_layout(guides = "collect"),
  width = 12,
  height = 5,
  dpi = 300,
  device = ragg::agg_png
)
