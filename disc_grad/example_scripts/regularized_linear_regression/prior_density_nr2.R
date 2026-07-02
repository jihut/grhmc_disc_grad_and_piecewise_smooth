rm(list = ls())

# Want to find an expression for the slab (continuous part) of the prior

# First: Look how samples directly from prior looks like for a given set of mu and rho

# Fix prior Var(beta | beta != 0) to 1, tweak P(beta==0) by tweaking mu and rho.  
cond_prior_var <- 5
prior_prob_vec <- c(0.25, 0.5, 0.95)

solve_for_mu_rho_given_p <- function(p) { # p = prior probability of beta
  
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
  
  return(c(mu = mu_given_rho_and_p(rho_sol, p), rho = rho_sol))
  
}

mu_rho_sols <- as.matrix(sapply(prior_prob_vec, solve_for_mu_rho_given_p))
mu_rho_sols[1, ] <- ifelse(abs(mu_rho_sols[1, ]) < 1e-12, 0, mu_rho_sols[1, ])

# The Dirac delta spike part happens with probability pnorm(-mu / rho) ^ 2, i.e. the prior probability chosen above

# For the slab (continuous) part, there are three terms
# Two cross terms which essentially are some part of pdfs due to Dirac delta
# One term that actually requires computing the integral from convolution formula

prior_beta_pdf <- function(beta) {
  
  # pnorm(-mu / rho) * # the cross terms 
  #   ifelse(
  #     beta <= 0, 
  #     dnorm(-beta, mean = mu, sd = rho),
  #     dnorm(beta, mean = mu, sd = rho)
  #   ) + 
  #   exp(-beta ^ 2 / (4 * rho ^ 2)) / # the convolution term
  #   (4 * sqrt(pi) * rho) * 
  #   ifelse(
  #     beta <= 0,
  #     2 * pnorm((2 * mu + beta) / (sqrt(2) * rho)),
  #     2 * pnorm((2 * mu - beta) / (sqrt(2) * rho))
  #   )
  
  ifelse(
    beta <= 0,
    pnorm(-mu / rho) *
      dnorm(-beta, mean = mu, sd = rho) +
      4 * sqrt(pi) * rho * (dnorm(beta, mean = 0, sd = 2 * rho) ^ 2) *
      pnorm((2 * mu + beta) / (sqrt(2) * rho))
    ,
    pnorm(-mu / rho) *
      dnorm(beta, mean = mu, sd = rho) +
      4 * sqrt(pi) * rho * (dnorm(beta, mean = 0, sd = 2 * rho) ^ 2) *
      pnorm((2 * mu - beta) / (sqrt(2) * rho))
  )
  
}

grid <- seq(from = -15, to = 15, by = 0.01)

par(mfrow = c(1, 3))

# 1. Try first combination corresponding to P(beta == 0) = 0.25

set.seed(42)

mu <- mu_rho_sols[1, 1]
rho <- mu_rho_sols[2, 1]

beta_plus_nr1 <- rnorm(10000000, mean = mu, sd = rho)
beta_minus_nr1 <- rnorm(10000000, mean = mu, sd = rho)

beta_nr1 <- pmax(0, beta_plus_nr1) - pmax(0, beta_minus_nr1)
hist(beta_nr1, breaks = 100, probability = T, xlim = c(-10, 10), ylim = c(0, 0.3), xlab = expression(beta), main = expression(P(beta == 0) == 0.25))

lines(grid, prior_beta_pdf(grid), col = "red")

# # Compare empirical CDF to CDF obtained by numerical integration for the slab part
# 
# normalized_prior_beta_pdf <- function(beta) { # slab part
#   
#   prior_beta_pdf(beta) / (1 - pnorm(-mu / rho) ^ 2)
#   
# }
# 
# ecdf_slab_nr1 <- ecdf(beta_nr1[beta_nr1 != 0])
# num_cdf_slab <- sapply(1:length(grid), function(i) integrate(normalized_prior_beta_pdf, -Inf, grid[i])$value)
# 
# plot(ecdf_slab_nr1)
# lines(grid, num_cdf_slab, col = "red", lty = 2)

# 2. Try first combination corresponding to P(beta == 0) = 0.50

set.seed(42)

mu <- mu_rho_sols[1, 2]
rho <- mu_rho_sols[2, 2]

beta_plus_nr2 <- rnorm(10000000, mean = mu, sd = rho)
beta_minus_nr2 <- rnorm(10000000, mean = mu, sd = rho)

beta_nr2 <- pmax(0, beta_plus_nr2) - pmax(0, beta_minus_nr2)
hist(beta_nr2, breaks = 100, probability = T, xlim = c(-10, 10), ylim = c(0, 0.3), xlab = expression(beta), main = expression(P(beta == 0) == 0.50))

lines(grid, prior_beta_pdf(grid), col = "red")

# 3. Try first combination corresponding to P(beta == 0) = 0.95

set.seed(42)

mu <- mu_rho_sols[1, 3]
rho <- mu_rho_sols[2, 3]

beta_plus_nr3 <- rnorm(10000000, mean = mu, sd = rho)
beta_minus_nr3 <- rnorm(10000000, mean = mu, sd = rho)

beta_nr3 <- pmax(0, beta_plus_nr3) - pmax(0, beta_minus_nr3)
hist(beta_nr3, breaks = 100, probability = T, xlim = c(-10, 10), ylim = c(0, 0.3), xlab = expression(beta), main = expression(P(beta == 0) == 0.95))

lines(grid, prior_beta_pdf(grid), col = "red")
