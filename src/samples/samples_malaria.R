# Generate samples from the prior distributions
# I think this might be stored in log space, but I have to check
din_prior_sample <- function(parameters, num_samples) {
  mean = parameters$din.mean # 0
  sd = parameters$din.sd # 2
  return(rnorm(num_samples, mean = mean, sd = sd))
}

logit_phi_prior_sample <- function(parameters, num_samples) {
  mean = parameters$phi.mean # 0
  sd = parameters$din.sd # 1
  return(rnorm(num_samples, mean = mean, sd = sd))
}

logit_eta0_prior_sample <- function(parameters, num_samples) {
  mean = parameters$eta0.mean # 0
  sd = parameters$eta0.sd # 1
  return(rnorm(num_samples, mean = mean, sd = sd))
}

log_sigma_prior_sample <- function(parameters, num_samples) {
  mean = parameters$sigma.mean # 10
  sd = parameters$sigma.sd # 4
  
  return(rnorm(num_samples, mean = mean, sd = sd))
}