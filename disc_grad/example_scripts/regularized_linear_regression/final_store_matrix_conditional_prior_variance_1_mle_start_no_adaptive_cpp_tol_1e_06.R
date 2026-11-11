rm(list = ls())
library(doParallel)
library(MASS)
library(brms)
library(ggplot2)
library(dplyr)
library(rstan)
library(foreach)
library(gridExtra)

source("disc_grad/implementation_scripts/general_scripts/cpp/grhmc_discontinuous_gradient_with_relevant_indices_function_cpp.R")

grhmc_wrapper_file_path <- normalizePath(
  "disc_grad/implementation_scripts/general_scripts/cpp/grhmc_discontinuous_gradient_with_relevant_indices_function_cpp.R"
)
model_file_path <- normalizePath(
  "disc_grad/example_scripts/regularized_linear_regression/regularized_linear_model.cpp"
)
model_cache <- normalizePath(
  "disc_grad/example_scripts/regularized_linear_regression/regularized_linear_model_cache",
  mustWork = F
)
dir.create(model_cache, showWarnings = FALSE)
Rcpp::sourceCpp(model_file_path, cacheDir = model_cache, rebuild = T)

trajectory_length <- 10000
n_samples_per_iteration <- n_samples <- 100000
n_trajectories <- n_iterations <- 10
n_cores <- 10
head(Boston)

x <- as.matrix(Boston[, -ncol(Boston)]) # 13 covariates
y <- Boston$medv
n_rows <- length(y)
length_beta <- ncol(x)

scale_x <- scale(x)
scale_y <- scale(y)
mean_y <- mean(y)

Boston_scale <- data.frame(scale_x, medv = scale_y)

# Fix prior Var(beta | beta != 0) to 1, tweak P(beta==0) by tweaking mu and rho.  
cond_prior_var <- 1
# prior_prob_vec <- seq(from = 0.005, to = 0.995, by = 0.005)
prior_prob_vec <- exp(seq(from = log(0.01), to = log(0.95), length.out = 50))

solve_for_mu_rho_given_p <- function(p) { # p = prior probability of beta == 0
  
  # Given the prior probability of beta and conditional variance --> want to solve for mu and rho
  
  mu_given_rho_and_p <- function(rho, p) { # use this equation to express mu w.r.t rho and p to have a one-dimensional root finding problem
    
    obj <- function(mu) {
      (-1 + 2 * pnorm(sqrt(2) / rho * mu * sqrt(2) / 2) - 1) ^ 2 / 4 - p
    }
    
    return(uniroot(obj, lower = -100, upper = 100, tol = 1e-12)$root)
    
  }
  
  cond_var_given_mu_and_rho <- function(mu, rho) { # the conditional prior variance based on mu and rho
    
    mu_star_1 <- rho * dnorm(mu / rho) + mu * pnorm(mu / rho)
    mu_star_2 <- rho * mu * dnorm(mu / rho) + (mu ^ 2 + rho ^ 2) * pnorm(mu / rho)
    uncond_var_beta <- 2 * (mu_star_2 - mu_star_1 ^ 2)
    cond_var_beta <- uncond_var_beta / (1 - pnorm(-mu / rho) ^ 2)
    return(cond_var_beta)
    
  }
  
  obj_for_rho <- function(rho) { # the root function for rho using the desired conditional prior variance
    cond_var_given_mu_and_rho(mu_given_rho_and_p(rho, p), rho) - cond_prior_var
  }
  
  rho_sol <- uniroot(obj_for_rho, lower = 0.01, upper = 10, tol = 1e-12)$root # solve for rho based on desired conditional prior variance
  mu_sol <- mu_given_rho_and_p(rho_sol, p) # use the result above to obtain the corresponding mu based on the desired prior probability
  
  mu_sol <- ifelse(abs(mu_sol) < 1e-13, 0, mu_sol)
  
  return(list(mu = mu_given_rho_and_p(rho_sol, p), rho = rho_sol))
  
}

system.time(
  
  for (i in 1:length(prior_prob_vec)) {
    progress_log <- normalizePath("disc_grad/example_scripts/regularized_linear_regression/progress.log", mustWork = FALSE)
    if (file.exists(progress_log)) file.remove(progress_log)
    
    print(paste0("i: ", i))
    print(paste0("Current prior probability of beta equal to zero: ", prior_prob_vec[i]))
    
    mu_rho_sols <- solve_for_mu_rho_given_p(prior_prob_vec[i])
    
    mu <- mu_rho_sols$mu
    rho <- mu_rho_sols$rho
    gamma_mu <- 0
    gamma_sigma <- 1
    
    init_cluster <- parallel::makeCluster(n_cores, outfile = progress_log)
    parallel::clusterExport(init_cluster, c("grhmc_wrapper_file_path", "model_file_path", "model_cache"))
    for (k in seq_along(init_cluster)) {
      parallel::clusterEvalQ(init_cluster[k], { # load cpp build one worker by one to avoid crash
        library(Rcpp)
        source(grhmc_wrapper_file_path)
        Rcpp::sourceCpp(model_file_path, cacheDir = model_cache, rebuild = FALSE)
        NULL
      })
    }
    doParallel::registerDoParallel(init_cluster)
    doRNG::registerDoRNG(seed = 42)
    
    samples_matrix <- foreach::foreach(
      i = 1:n_trajectories, 
      .combine = rbind,
      .noexport = c("grhmc_linear_model_cpp",
                    "grhmc_discontinuous_gradient_function_cpp")
    ) %dopar% {
      og_regmod <- lm(medv ~ ., data = Boston_scale)
      
      beta_plus_initial <- ifelse(og_regmod$coefficients[-1] < 0, mu, abs(og_regmod$coefficients[-1])) # start at MLE
      beta_minus_initial <- ifelse(og_regmod$coefficients[-1] > 0, mu, abs(og_regmod$coefficients[-1]))
      q_initial <- c(beta_plus_initial, beta_minus_initial, sum((og_regmod$residuals)^2) / og_regmod$df.residual)
      p_initial <- rnorm(length_beta * 2 + 1)
      u_initial <- rexp(1)
      
      test_run_warmup <- grhmc_discontinuous_gradient_function_cpp(
        model_list = function(settings) grhmc_linear_model_cpp(settings, mu = mu, rho = rho, gamma_mu = gamma_mu, gamma_sigma = gamma_sigma, x_matrix = scale_x, y_vector = scale_y),
        lambda = 0.2,
        T = trajectory_length,
        n_samples = n_samples,
        q_initial,
        p_initial,
        u_initial, 
        random_state = NULL,
        rtol = 1e-6,
        atol = 1e-6,
        h_max = 1.0,
        last_root_offset_lin_root_finder = 1e-10,
        precision_real_root_lin_root_finder = 1.0e-13,
        coeff_tol_real_root_lin_root_finder = 1.0e-13,
        push_eps_tol = 1e-4,
        relevant_indices = NULL,
        return_output_from_ode = FALSE,
        trajectory_info_string = paste0("warmup ", i),
        progress_ratio = 0.05
      )
      
      new_p_initial <- rnorm(length_beta * 2 + 1)
      new_u_initial <- rexp(1)
      
      test_run_sampling <- grhmc_discontinuous_gradient_function_cpp(
        model_list = function(settings) grhmc_linear_model_cpp(settings, mu = mu, rho = rho, gamma_mu = gamma_mu, gamma_sigma = gamma_sigma, x_matrix = scale_x, y_vector = scale_y),
        lambda = 0.2,
        T = trajectory_length,
        n_samples = n_samples,
        test_run_warmup$q_final,
        new_p_initial,
        new_u_initial, 
        random_state = NULL,
        rtol = 1e-6,
        atol = 1e-6,
        h_max = 1.0,
        last_root_offset_lin_root_finder = 1e-10,
        precision_real_root_lin_root_finder = 1.0e-13,
        coeff_tol_real_root_lin_root_finder = 1.0e-13,
        push_eps_tol = 1e-4,
        relevant_indices = NULL,
        return_output_from_ode = FALSE,
        trajectory_info_string = paste0("sampling ", i),
        progress_ratio = 0.05
      )
      test_run_sampling$q_original_samples
    }
    parallel::stopCluster(init_cluster)
    
    boston_scale_grhmc_beta_store_matrix <- matrix(0, nrow = n_samples_per_iteration * n_iterations, ncol = length_beta) # store the scaled beta samples 
    
    for (j in 1:length_beta) {
      
      current_beta_vec <- pmax(0, samples_matrix[, j]) - pmax(0, samples_matrix[, length_beta + j])
      
      boston_scale_grhmc_beta_store_matrix[, j] <- current_beta_vec
      
    }
    
    boston_grhmc_original_scale_beta_matrix <- t(t(boston_scale_grhmc_beta_store_matrix) / apply(x, 2, sd)) * sd(y) # convert back to original scale
    
    boston_grhmc_ratio_samples_zero <- sapply(1:ncol(boston_grhmc_original_scale_beta_matrix), function(i) mean(boston_grhmc_original_scale_beta_matrix[, i] == 0)) # posterior probability of beta equal to zero
    
    boston_grhmc_original_scale_beta_cred_int <- matrix(nrow = ncol(boston_grhmc_original_scale_beta_matrix), ncol = 2)
    colnames(boston_grhmc_original_scale_beta_cred_int) <- c("lower", "upper")
    
    for (j in 1:ncol(boston_grhmc_original_scale_beta_matrix)) {
      
      boston_grhmc_original_scale_beta_cred_int[j, ] <- quantile(boston_grhmc_original_scale_beta_matrix[, j], probs = c(0.025, 0.975))
      
    }
    
    final_boston_grhmc_original_scale_beta_cred_int <- data.frame(
      boston_grhmc_original_scale_beta_cred_int, 
      posterior_median = apply(boston_grhmc_original_scale_beta_matrix, 2, median),
      variable = colnames(x),
      mu = mu, 
      rho = rho,
      prior_prob = prior_prob_vec[i],
      ratio_samples_equal_zero = boston_grhmc_ratio_samples_zero
    )
    
    if (i == 1) {
      store_matrix <- final_boston_grhmc_original_scale_beta_cred_int
    } else {
      store_matrix <- rbind(store_matrix, final_boston_grhmc_original_scale_beta_cred_int)
    }
    
  }
  
)

saveRDS(store_matrix, "disc_grad/example_scripts/regularized_linear_regression/final_store_matrix_conditional_prior_variance_1_mle_start_no_adaptive_cpp_tol_1e_06.RDS")
# store_matrix <- readRDS("disc_grad/example_scripts/regularized_linear_regression/final_store_matrix_conditional_prior_variance_1_mle_start_no_adaptive_cpp_tol_1e_06.RDS")

store_matrix

final_store_matrix <- store_matrix

plot_nr1 <- final_store_matrix %>% 
  ggplot(aes(x = prior_prob, y = posterior_median, color = variable, group = variable, fill = variable)) + 
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = variable), alpha = 0.2, linetype = 0) + 
  labs(title = latex2exp::TeX(r"($Var (\beta_i | \beta_i \neq 0$) = 1)")) + 
  xlab(latex2exp::TeX(r"($Prior \, P(\beta_i=0)$)")) + 
  ylab("Posterior median (original scale)")
plot_nr1

plot_nr2 <- final_store_matrix %>% 
  filter(variable != "nox") %>%
  ggplot(aes(x = prior_prob, y = posterior_median, color = variable, group = variable, fill = variable)) + 
  geom_line() + 
  labs(title = latex2exp::TeX(r"($Var (\beta_i | \beta_i \neq 0$) = 1)")) + 
  xlab(latex2exp::TeX(r"($Prior \, P(\beta_i=0)$)")) + 
  ylab("Posterior median (original scale)")
plot_nr2

plot_nr3 <- final_store_matrix %>% 
  ggplot(aes(x = prior_prob, y = ratio_samples_equal_zero, color = variable, group = variable, fill = variable)) + 
  geom_line() +
  labs(title = latex2exp::TeX(r"($Var (\beta_i | \beta_i \neq 0$) = 1)")) + 
  xlab(latex2exp::TeX(r"($Prior \, P(\beta_i=0)$)")) + 
  ylab(latex2exp::TeX(r"($Posterior \, P(\beta_i=0)$)"))
plot_nr3

grid.arrange(plot_nr1, plot_nr2, plot_nr3, ncol = 3, top = "Boston dataset")

# If we want to plot P(\beta_i != 0)

final_store_matrix <- store_matrix %>%
  mutate(
    prior_prob_non_zero = 1 - prior_prob
  )

plot_nr1 <- final_store_matrix %>% 
  ggplot(aes(x = prior_prob_non_zero, y = posterior_median, color = variable, group = variable, fill = variable)) + 
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = variable), alpha = 0.2, linetype = 0) + 
  labs(title = latex2exp::TeX(r"($Var (\beta_i | \beta_i \neq 0$) = 1)")) + 
  xlab(latex2exp::TeX(r"($Prior \, P(\beta_i \neq 0)$)")) + 
  ylab("Posterior median (original scale)")
plot_nr1

plot_nr2 <- final_store_matrix %>% 
  filter(variable != "nox") %>%
  ggplot(aes(x = prior_prob_non_zero, y = posterior_median, color = variable, group = variable, fill = variable)) + 
  geom_line() + 
  labs(title = latex2exp::TeX(r"($Var (\beta_i | \beta_i \neq 0$) = 1)")) + 
  xlab(latex2exp::TeX(r"($Prior \, P(\beta_i \neq 0)$)")) + 
  ylab("Posterior median (original scale)")
plot_nr2

plot_nr3 <- final_store_matrix %>% 
  ggplot(aes(x = prior_prob_non_zero, y = ratio_samples_equal_zero, color = variable, group = variable, fill = variable)) + 
  geom_line() +
  labs(title = latex2exp::TeX(r"($Var (\beta_i | \beta_i \neq 0$) = 1)")) + 
  xlab(latex2exp::TeX(r"($Prior \, P(\beta_i \neq 0)$)")) + 
  ylab(latex2exp::TeX(r"($Posterior \, P(\beta_i=0)$)"))
plot_nr3

grid.arrange(plot_nr1, plot_nr2, plot_nr3, ncol = 3, top = "Boston dataset")

# If we want to plot log P(\beta_i != 0)

final_store_matrix <- store_matrix %>%
  mutate(
    prior_prob_non_zero = 1 - prior_prob
  ) %>%
  mutate(
    log_prior_prob_non_zero = log(prior_prob_non_zero)
  )

plot_nr1 <- final_store_matrix %>% 
  ggplot(aes(x = log_prior_prob_non_zero, y = posterior_median, color = variable, group = variable, fill = variable)) + 
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = variable), alpha = 0.2, linetype = 0) + 
  labs(title = latex2exp::TeX(r"($Var (\beta_i | \beta_i \neq 0$) = 1)")) + 
  xlab(latex2exp::TeX(r"($Prior \, log(P(\beta_i \neq 0))$)")) + 
  ylab("Posterior median (original scale)")
plot_nr1

plot_nr2 <- final_store_matrix %>% 
  filter(variable != "nox") %>%
  ggplot(aes(x = log_prior_prob_non_zero, y = posterior_median, color = variable, group = variable, fill = variable)) + 
  geom_line() + 
  labs(title = latex2exp::TeX(r"($Var (\beta_i | \beta_i \neq 0$) = 1)")) + 
  xlab(latex2exp::TeX(r"($Prior \, log(P(\beta_i \neq 0))$)")) + 
  ylab("Posterior median (original scale)")
plot_nr2

plot_nr3 <- final_store_matrix %>% 
  ggplot(aes(x = log_prior_prob_non_zero, y = ratio_samples_equal_zero, color = variable, group = variable, fill = variable)) + 
  geom_line() +
  labs(title = latex2exp::TeX(r"($Var (\beta_i | \beta_i \neq 0$) = 1)")) + 
  xlab(latex2exp::TeX(r"($Prior \, log(P(\beta_i \neq 0))$)")) + 
  ylab(latex2exp::TeX(r"($Posterior \, P(\beta_i=0)$)"))
plot_nr3

grid.arrange(plot_nr1, plot_nr2, plot_nr3, ncol = 3, top = "Boston dataset")
