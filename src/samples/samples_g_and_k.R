pacman::p_load(pacman, gk)

####################### Priors ################################

a_prior_sample <- function(parameters, num_samples) {
  samples <- runif(num_samples, 0, 10)
  # Need to convert them to (-infty, infty) scale
  return(qnorm(samples/10))
}

b_prior_sample <- function(parameters, num_samples) {
  samples <- runif(num_samples, 0, 10)
  # Need to convert them to (-infty, infty) scale
  return(qnorm(samples/10))
}

g_prior_sample <- function(parameters, num_samples) {
  samples <- runif(num_samples, 0, 10)
  # Need to convert them to (-infty, infty) scale
  return(qnorm(samples/10))
}

k_prior_sample <- function(parameters, num_samples) {
  samples <- runif(num_samples, 0, 10)
  # Need to convert them to (-infty, infty) scale
  return(qnorm(samples/10))
}

####################### Likelihood ###########################

likelihood_g_and_k <- function(parameters) {
  # Need to convert parameters back to the scale [0, 10]
  a <- pnorm(parameters$a)*10
  b <- pnorm(parameters$b)*10
  g <- pnorm(parameters$g)*10
  k <- pnorm(parameters$k)*10
  
  rgk(1, a, b, g, k)
}

pnorm(0)