loglike_malaria <- function(y, likelihood_means, parameters) {
  # print(y)
  # print(y - likelihood_means)
  sigma = parameters$sigma
  return(sum(dnorm(y, mean = likelihood_means, sd = sigma, log = T)))
}