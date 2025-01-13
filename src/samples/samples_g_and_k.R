pacman::p_load(pacman, gk, spatstat.utils)

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
  
  return(rgk(100, a, b, g, k))
}

likelihood_g_and_k_summary <- function(parameters) {
  # Need to convert parameters back to the scale [0, 10]
  a <- pnorm(parameters$a)*10
  b <- pnorm(parameters$b)*10
  g <- pnorm(parameters$g)*10
  k <- pnorm(parameters$k)*10
  
  observations <- rgk(1000, a, b, g, k)

  order_stats_sequence <- seq(1, 1000, length.out = 100)
  return(orderstats(observations, order_stats_sequence))
}

params <- list(a = 0, b = 0, g = 0, k = 0)
likelihood_g_and_k_summary(params)

x <- runif(10)
orderstats(x, 2)