pacman::p_load(pacman, mvtnorm, purrr)
source('src/pdfs_normal.R')
source('src/samples_normal.R')

lupost_normal <- function(y, theta.init, x, params) {
  
  alpha.sd <- params$alpha.sd
  sigma2.sd <- params$sigma2.sd
  
  alpha.init <- theta.init[1]
  logsigma2.init <- theta.init[2]
  
  return(alpha_logprior_pdf(alpha.init, alpha.sd) + 
           logsigma2_logprior_pdf(logsigma2.init, sigma2.sd) + 
           loglike_pdf(y, alpha.init, x, exp(logsigma2.init)))
}

prior_params = list(alpha.sd = 2, sigma2.sd = 5)

normal_mcmc <- function(iterations, true_params, prior_params) {
  
  # 1. Generate data given true parameters
  simulated_data <- generate_data(iterations, true_params)
  
  # 2. Create a partial function fixing data and parameters
  lupost_normal_mcmc <- partial(lupost_normal, y = simulated_data, 
                                x = 1, params = prior_params)
  
  # 3. Draw parameters from prior distribution
  alpha.init <- alpha.prior_sample(prior_params$alpha.sd, 1)
  logsigma2.init <- logsigma2_prior_sample(prior_params$sigma2.sd, 1)
  theta.init <- c(alpha.init, logsigma2.init)
  
  # 4. Call on MCMC function
  
}