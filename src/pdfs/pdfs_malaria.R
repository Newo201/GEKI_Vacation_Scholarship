loglike_malaria <- function(y, likelihood_means, sigma, log_obs = T) {
  if (log_obs) {
    return(sum(dnorm(y, mean = likelihood_means, sd = sigma, log = T)))
  }
  else {
    return(sum(dlnorm(y, mean = likelihood_means, sd = sigma, log = T)))
  }
}