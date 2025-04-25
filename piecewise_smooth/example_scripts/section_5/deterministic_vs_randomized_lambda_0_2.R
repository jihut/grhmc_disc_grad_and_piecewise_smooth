non_randomized_store_matrix <- readRDS("piecewise_smooth/example_scripts/section_5/normal_2d_unit_circle_non_randomized_reflection_non_adaptive_lambda_0_2.RDS")
randomized_store_matrix <- readRDS("piecewise_smooth/example_scripts/section_5/normal_2d_unit_circle_randomized_reflection_non_adaptive_lambda_0_2.RDS")

par(mfrow = c(1, 2))
hist(non_randomized_store_matrix[, 1], probability = T, breaks = 100, xlab = expression(q[1]), main = "Deterministic")
lines(seq(from = -10, to = 10, by = 0.01), marginal_q1_pdf(seq(from = -10, to = 10, by = 0.01)), col = "red")

hist(randomized_store_matrix[, 1], probability = T, breaks = 100, xlab = expression(q[1]), main = "Randomized")
lines(seq(from = -10, to = 10, by = 0.01), marginal_q1_pdf(seq(from = -10, to = 10, by = 0.01)), col = "red")
