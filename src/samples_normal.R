###################### Sampling From Priors #################################

# Generate samples from the prior distributions
alpha_prior_sample <- function(alpha.sd, iterations) {
  return(rnorm(iterations, mean = 0, sd = alpha.sd))
}

# We need to have a prior sample that is from (-infinity, infinity)
# but we also need to restrict sigma2 to R+. Therefore we sample the log
# prior from a normal distribution
logsigma2_prior_sample <- function(sigma2.sd, iterations) {
  
  return(rnorm(iterations, mean = 0, sd = sigma2.sd))
}

###################### Sampling From Likelihood ###############################

# Generate samples from the normal distribution
likelihood_sample <- function(alpha, x, sigma2) {
  dim = length(x)
  return(rmvnorm(1, mean = alpha*x, sigma = sigma2*diag(dim)))
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