rm(list = ls())
source("disc_grad/illustrations/illustration_global_error/comparison_nr1.R")
source("disc_grad/illustrations/illustration_global_error/comparison_nr2.R")
source("disc_grad/illustrations/illustration_global_error/comparison_nr3.R")

pdf("disc_grad/illustrations/illustration_global_error/base_comparison_plot.pdf", width = 11.7, height = 8.3)

par(mfrow = c(1, 3))

plot(
  h_vec_c_0_1, 
  leapfrog_store_error_c_0_1, 
  log = "xy", 
  ylim = c(1e-11, 1), 
  ylab = "L2 error at final state", 
  xlab = "h (step size)", 
  main = "c = 0.1",
  cex.axis = 1.5,
  cex.lab = 1.5,
  cex.main = 1.5
)
points(h_vec_c_0_1, grhmc_nr1_store_error_c_0_1, col = "red")
points(h_vec_c_0_1, grhmc_nr2_store_error_c_0_1, col = "blue", pch = 3)
lines(h_vec_c_0_1, 0.068 * h_vec_c_0_1 ^ 3, lty = 2, col = "green")
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
  cex = 1.25
)

plot(
  h_vec_c_1, 
  leapfrog_store_error_c_1, 
  log = "xy", 
  ylim = c(1e-10, 1), 
  ylab = "L2 error at final state", 
  xlab = "h (step size)", 
  main = "c = 1",
  cex.axis = 1.5,
  cex.lab = 1.5,
  cex.main = 1.5
)
points(h_vec_c_1, grhmc_nr1_store_error_c_1, col = "red")
points(h_vec_c_1, grhmc_nr2_store_error_c_1, col = "blue", pch = 3)
lines(h_vec_c_1, 0.25 * h_vec_c_1 ^ 3, lty = 2, col = "green")
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
  cex = 1.25
)

plot(
  h_vec_c_10, 
  leapfrog_store_error_c_10, 
  log = "xy", 
  ylim = c(1e-10, 1e-1), 
  ylab = "L2 error at final state", 
  xlab = "h (step size)", 
  main = "c = 10",
  cex.axis = 1.5,
  cex.lab = 1.5,
  cex.main = 1.5
)
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
  col = c("black", "red", "blue", "green"),
  cex = 1.25
)

dev.off()
