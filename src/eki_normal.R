pacman::p_load(pacman, MASS)

# Generate samples from the prior distributions
alpha_prior_sample <- function(alpha.sd, iterations) {
  return(rnorm(iterations, mean = 0, sd = alpha.sd))
}

# We need to have a prior sample that is from (-infinity, infinity)
# but we also need to restrict sigma2 to R+. Therefore we sample the log
# prior from a normal distribution
sigma2_prior_sample <- function(sigma2.sd, iterations) {
  
  log.sigma2 <- rnorm(iterations, mean = 0, sd = sigma2.sd)
  return(exp(log.sigma2))
}

# Generate samples from the normal distribution
likelihood_sample <- function(alpha, x, sigma2) {
  dim = length(x)
  return(mvrnorm(1, mu = alpha*x, Sigma = sigma2*diag(dim)))
}

# Generate simulated data using the true parameters
generate_data <- function(iterations, parameters) {
  alpha.true <- parameters$alpha
  x.true <- parameters$x
  sigma2.true <- parameters$sigma2
  samples <- rep(NA, iterations)
  print(alpha.true, x.true, sigma2.true)
  for (i in 1:iterations) {
    samples[i] <- likelihood_sample(alpha.true, x.true, sigma2.true)
  }
  return(samples)
}

eki_normal <- function(iterations, parameters) {
  
  x.true <- parameters$x
  alpha.sd <- parameters$alpha.sd
  sigma2.sd <- parameters$sigma2.sd
  
  d_y <- length(x.true)
  
  # Simulate data using true parameters
  simulated_data <- generate_data(iterations, parameters)
  
  # Sample from the prior distribution
  prior_samples <- matrix(nrow = iterations, ncol = 2)
  prior_samples[, 1] <- alpha_prior_sample(alpha.sd, iterations)
  # Question: do I have to move particles in log(sigma2) space or in sigma2 space
  prior_samples[, 2] <- sigma2_prior_sample(sigma2.sd, iterations)
  
  # Until we reach a temperature of one do the following
  
    # Sample from the likelihood
  
  likelihood_samples <- matrix(nrow = iterations, ncol = d_y)
  for (i in 1:iterations) {
    likelihood_samples[i, ] <- likelihood_sample(prior_samples[i, 1], x.true, prior_samples[i, 2])
  }
  
    # Calculate the covariance matrices
  C_xx = cov(prior_samples)
  print(C_xx)
  C_yy = cov(likelihood_samples)
  C_xy = cov(prior_samples, likelihood_samples)
  C_yx = cov(likelihood_samples, prior_samples)
  
  C_y_given_x = C_yy - C_yx %*% solve(C_xx) %*% C_xy
  return(C_y_given_x)
  
    # Generate perturbations
  
    # Move the particles
  
    # Calculate the next temperature
}

parameters <- list(alpha = 2, sigma2 = 5, x = 5, alpha.sd = 5, sigma2.sd = 2)
eki_normal(100, parameters)

alpha_test <- alpha_prior_sample(5, 1000)
hist(alpha_test)

cov(alpha_test, alpha_test)


