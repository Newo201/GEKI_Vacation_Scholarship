loglike_malaria <- function(y, likelihood_means, parameters, log_obs = T) {
  sigma = parameters$sigma
  # A check to determine if data is on the log scale
  if (log_obs) {
    return(sum(dnorm(y, mean = likelihood_means, sd = sigma, log = T)))
  } else {
    return(sum(dlnorm(y, mean = likelihood_means, sd = sigma, log = T)))
  }

}