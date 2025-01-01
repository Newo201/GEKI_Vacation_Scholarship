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

# Generate simulated data using the true parameters
generate_data <- function(iterations, parameters) {
  alpha.true <- parameters['alpha']
  x.true <- parameters['x']
  sigma2.true <- parameters['sigma2']
  samples <- rep(NA, iterations)
  for (i in 1:iterations) {
    samples[i] <- likelihood_sample(alpha.true, x.true, sigma2.true)
  }
  return(samples)
}

eki_normal <- function(iterations, parameters) {
  
  x.true <- parameters['x']
  
  # Simulate data using true parameters
  simulated_data <- generate_data(iterations, parameters)
  
  # Sample from the prior distribution
  
  # Until we reach a temperature of one do the following
  
    # Sample from the likelihood
  
    # Calculate the covariance matrices
  
    # Generate perturbations
  
    # Move the particles
  
    # Calculate the next temperature
}



