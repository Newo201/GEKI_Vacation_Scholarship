pacman::p_load(pacman, MASS)

# Generate samples from the prior distributions
alpha_prior_sample <- function(alpha.sd) {
  return(rnorm(1, mean = 0, sd = alpha.sd))
}

# We need to have a prior sample that is from (-infinity, infinity)
# but we also need to restrict sigma2 to R+. Therefore we sample the log
# prior from a normal distribution
sigma2_prior_sample <- function(sigma2.sd) {
  
  log.sigma2 <- rnorm(1, mean = 0, sd = sigma2.sd)
  return(exp(log.sigma2))
}

# Generate samples from the normal distribution
likelihood_sample <- function(alpha, x, sigma2) {
  dim = length(x)
  return(mvrnorm(1, mu = alpha*x, Sigma = sigma2*diag(dim)))
}

likelihood_sample(1, c(1,0,1), 5)
alpha_prior_sample(5)
sigma2_prior_sample(2)

