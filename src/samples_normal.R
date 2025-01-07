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
  sigma = parameters$sigma
  dim = length(x)
  return(rmvnorm(num_samples, mean = alpha*x, sigma = (sigma^2)*diag(dim)))
}

# # Generate simulated data using the true parameters
# generate_data <- function(iterations, parameters) {
#   alpha.true <- parameters$alpha
#   x.true <- parameters$x
#   sigma2.true <- parameters$sigma2
#   samples <- matrix(nrow = iterations, ncol = length(x.true))
#   print(alpha.true, x.true, sigma2.true)
#   for (i in 1:iterations) {
#     samples[i, ] <- likelihood_sample(alpha.true, x.true, sigma2.true)
#   }
#   return(samples)
# }


# parameters = list(alpha = 2, sigma2 = 2, x = c(1, 1))
# test <- generate_data(100, parameters)
# likelihood_sample()

test <- rmvnorm(10, mean = c(1, 0, 1), sigma = diag(3))