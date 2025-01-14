loglike_malaria <- function(y, parameters) {
  likelihood_means = parameters$likelihood_means
  sigma = parameters$sigma
  return(sum(dnorm(y, mean = likelihood_means, sd = sigma, log = T)))
}