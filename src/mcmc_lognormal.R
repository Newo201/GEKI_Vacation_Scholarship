lupost_lognormal <- function(y, theta.init, x, prior_params) {
  
  alpha.init <- theta.init[1]
  logsigma2.init <- theta.init[2]
  
  like_params <- list(alpha = alpha.init, x = x, sigma = sqrt(exp(logsigma2.init)))
  
  return(alpha_logprior_pdf(alpha.init, prior_params) + 
           logsigma2_logprior_pdf(logsigma2.init, prior_params) + 
           loglike_lognormal(y, like_params))
}

lognormal_mcmc <- function(true_data, true_params, prior_params, iterations = 1e4) {
  
  
  # 1. Create a partial function fixing data and parameters
  lupost_lognormal_mcmc <- partial(lupost_lognormal, y = true_data, 
                                x = true_params$x, prior_params = prior_params)
  
  # 2. Draw parameters from prior distribution
  alpha.init <- alpha_prior_sample(prior_params, 1)
  logsigma2.init <- logsigma2_prior_sample(prior_params, 1)
  theta.init <- c(alpha.init, logsigma2.init)
  
  # 3. Call on MCMC function
  chain <- metrop(lupost_lognormal_mcmc, theta.init, iterations)
  return(chain)
  
}