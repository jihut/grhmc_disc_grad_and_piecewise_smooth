rm(list = ls())

library(MASS)
library(brms)
library(ggplot2)
library(dplyr)
library(rstan)
library(foreach)
library(gridExtra)

n_iterations <- 10
n_samples_per_iteration <- 10000

head(Boston)

x <- as.matrix(Boston[, -ncol(Boston)]) # 13 covariates
y <- Boston$medv
n_rows <- length(y)
length_beta <- ncol(x)

scale_x <- scale(x)
scale_y <- scale(y)
mean_y <- mean(y)

disc_grad_samples <- readRDS("disc_grad/example_scripts/section_6/regularized_linear_regression/samples_disc_grad_run_example.RDS")
sticky_samples <- readRDS("disc_grad/example_scripts/section_6/regularized_linear_regression/samples_sticky_run_example.RDS")

# Discontinuous gradient setup

final_disc_grad_samples <- matrix(0, nrow = n_samples_per_iteration * n_iterations, ncol = length_beta) # store the scaled beta samples 

for (j in 1:length_beta) {
  
  current_beta_vec <- pmax(0, disc_grad_samples[, j]) - pmax(0, disc_grad_samples[, length_beta + j])
  
  final_disc_grad_samples[, j] <- current_beta_vec
  
}

final_disc_grad_samples <- t(t(final_disc_grad_samples) / apply(x, 2, sd)) * sd(y) # convert back to original scale
final_disc_grad_samples <- cbind(final_disc_grad_samples, sort(rep(1:n_iterations, n_samples_per_iteration)))
colnames(final_disc_grad_samples) <- c(colnames(x), "iteration")

# Sticky setup

final_sticky_samples <- t(t(sticky_samples[, 1:length_beta]) / apply(x, 2, sd)) * sd(y) # convert back to original scale
final_sticky_samples <- cbind(final_sticky_samples, sort(rep(1:n_iterations, n_samples_per_iteration)))
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

