###################### Sampling From Priors #################################

# Generate samples from the prior distributions
alpha_prior_sample <- function(parameters, num_samples) {
  mean = parameters$alpha.mean
  sd = parameters$alpha.sd
  return(rnorm(num_samples, mean = mean, sd = sd))
}

# We need to have a prior sample that is from (-infinity, infinity)
# but we also need to restrict sigma2 to R+. Therefore we sample the log
# prior from a normal distribution
logsigma2_prior_sample <- function(parameters, num_samples) {
  
  mean = parameters$sigma2.mean
  sd = parameters$sigma2.sd
  
  return(rnorm(num_samples, mean = mean, sd = sd))
}

###################### Sampling From Likelihood ###############################

# Generate samples from the normal distribution
likelihood_normal <- function(likelihood_mean, parameters) {
  dim <- length(likelihood_mean)
  sigma2 = parameters$sigma**2
  return(rmvnorm(1, mean = likelihood_mean, sigma = sigma2*diag(dim)))
}