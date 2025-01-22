loglike_malaria <- function(y, likelihood_means, parameters, log_obs = T) {
  # print(y)
  # print(y - likelihood_means)
  sigma = parameters$sigma
  if (log_obs) {
    return(sum(dnorm(y, mean = likelihood_means, sd = sigma, log = T)))
  } else {
    return(sum(dlnorm(y, mean = likelihood_means, sd = sigma, log = T)))
  }

}