loglike_malaria <- function(y, likelihood_means, sigma) {
  return(sum(dnorm(y, mean = likelihood_means, sd = sigma, log = T)))
}