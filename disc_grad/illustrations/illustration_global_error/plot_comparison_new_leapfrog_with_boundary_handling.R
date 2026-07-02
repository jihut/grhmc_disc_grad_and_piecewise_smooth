rm(list = ls())
source("disc_grad/illustrations/illustration_global_error/comparison_nr1_leapfrog_with_boundary_handling_trial_nr1.R")
source("disc_grad/illustrations/illustration_global_error/comparison_nr2_leapfrog_with_boundary_handling_trial_nr1.R")
source("disc_grad/illustrations/illustration_global_error/comparison_nr3_leapfrog_with_boundary_handling_trial_nr1.R")

pdf("disc_grad/illustrations/illustration_global_error/base_comparison_plot_new_leapfrog_with_boundary_handling_trial_nr1.pdf", width = 14, height = 7)

par(mfrow = c(1, 3))
par(mar = c(5, 6, 4, 2))
plot(
  h_vec_c_0_1, 
  leapfrog_nr1_store_error_c_0_1, 
  log = "xy", 
  ylim = c(1e-11, 1), 
  ylab = "L2 error at final state", 
  xlab = "h (step size)", 
  main = "c = 0.1",
  cex = 2,
  cex.axis = 2,
  cex.lab = 2,
  cex.main = 2
)
points(h_vec_c_0_1, leapfrog_nr2_store_error_c_0_1, col = adjustcolor("orange", alpha.f = 1), cex = 2, pch = 3)
points(h_vec_c_0_1, grhmc_nr1_store_error_c_0_1, col = "red", cex = 2)
points(h_vec_c_0_1, grhmc_nr2_store_error_c_0_1, col = "blue", pch = 3, cex = 2)
lines(h_vec_c_0_1, 0.068 * h_vec_c_0_1 ^ 3, lty = 2, col = "green", lwd = 3)
lines(h_vec_c_0_1, 0.21 * h_vec_c_0_1 ^ 2, lty = 2, col = "darkgreen", lwd = 3)
legend(
  "bottomright",
  legend = c(
    "Leapfrog without boundary handling",
    "Leapfrog with boundary handling",
    "GRHMC without boundary handling",
    "GRHMC with boundary handling",
    latex2exp::TeX("$\\propto h^3$"),
    latex2exp::TeX("$\\propto h^2$")
  ),
  pch = c(1, 1, 1, 3, NA, NA),
  lty = c(NA, NA, NA, NA, 2, 2),
  col = c("black", "orange", "red", "blue", "green", "darkgreen"),
  cex = 1.3
)

plot(
  h_vec_c_1, 
  leapfrog_nr1_store_error_c_1, 
  log = "xy", 
  ylim = c(1e-10, 1), 
  ylab = "L2 error at final state", 
  xlab = "h (step size)", 
  main = "c = 1",
  cex = 2,
  cex.axis = 2,
  cex.lab = 2,
  cex.main = 2
)
points(h_vec_c_1, leapfrog_nr2_store_error_c_1, col = adjustcolor("orange", alpha.f = 1), cex = 2, pch = 3)
points(h_vec_c_1, grhmc_nr1_store_error_c_1, col = "red", cex = 2)
points(h_vec_c_1, grhmc_nr2_store_error_c_1, col = "blue", pch = 3, cex = 2)
lines(h_vec_c_1, 0.25 * h_vec_c_1 ^ 3, lty = 2, col = "green", lwd = 3)
lines(h_vec_c_1, 0.45 * h_vec_c_1 ^ 2, lty = 2, col = "darkgreen", lwd = 3)
legend(
  "bottomright",
  legend = c(
    "Leapfrog without boundary handling",
    "Leapfrog with boundary handling",
    "GRHMC without boundary handling",
    "GRHMC with boundary handling",
    latex2exp::TeX("$\\propto h^3$"),
    latex2exp::TeX("$\\propto h^2$")
  ),
  pch = c(1, 1, 1, 3, NA, NA),
  lty = c(NA, NA, NA, NA, 2, 2),
  col = c("black", "orange", "red", "blue", "green", "darkgreen"),
  cex = 1.3
)

plot(
  h_vec_c_10, 
  leapfrog_nr1_store_error_c_10, 
  log = "xy", 
  ylim = c(1e-10, 1e-1), 
  ylab = "L2 error at final state", 
  xlab = "h (step size)", 
  main = "c = 10",
  cex = 2,
  cex.axis = 2,
  cex.lab = 2,
  cex.main = 2
)
points(h_vec_c_10, leapfrog_nr2_store_error_c_10, col = adjustcolor("orange", alpha.f = 1), cex = 2, pch = 3)
points(h_vec_c_10, grhmc_nr1_store_error_c_10, col = "red", cex = 2)
points(h_vec_c_10, grhmc_nr2_store_error_c_10, col = "blue", pch = 3, cex = 2)
lines(h_vec_c_10, 105 * h_vec_c_10 ^ 3, lty = 2, col = "green", lwd = 3)
lines(h_vec_c_10, 16 * h_vec_c_10 ^ 2, lty = 2, col = "darkgreen", lwd = 3)
legend(
  "bottomright",
  legend = c(
    "Leapfrog without boundary handling",
    "Leapfrog with boundary handling",
    "GRHMC without boundary handling",
    "GRHMC with boundary handling",
    latex2exp::TeX("$\\propto h^3$"),
    latex2exp::TeX("$\\propto h^2$")
  ),
  pch = c(1, 1, 1, 3, NA, NA),
  lty = c(NA, NA, NA, NA, 2, 2),
  col = c("black", "orange", "red", "blue", "green", "darkgreen"),
  cex = 1.3
)

dev.off()
