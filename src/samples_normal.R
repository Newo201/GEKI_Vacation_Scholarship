###################### Sampling From Priors #################################

# Generate samples from the prior distributions
alpha_prior_sample <- function(alpha.sd, num_samples) {
  return(rnorm(num_samples, mean = 0, sd = alpha.sd))
}

# We need to have a prior sample that is from (-infinity, infinity)
# but we also need to restrict sigma2 to R+. Therefore we sample the log
# prior from a normal distribution
logsigma2_prior_sample <- function(sigma2.sd, num_samples) {
  
  return(rnorm(num_samples, mean = 0, sd = sigma2.sd))
}

###################### Sampling From Likelihood ###############################

# Generate samples from the normal distribution
likelihood_sample <- function(parameters, num_samples) {
  alpha = parameters$alpha
  x = parameters$x
  sigma2 = parameters$sigma**2
  dim = length(x)
  print(sigma2)
  # print(sigma2*diag(dim))
  return(rmvnorm(num_samples, mean = alpha*x, sigma = sigma2*diag(dim)))
}