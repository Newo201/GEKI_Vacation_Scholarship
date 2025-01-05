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

lupost_normal(1, c(0, 1), 1, prior_params)

lupost_normal_mcmc <- partial(lupost_normal, x = 1, params = prior_params)
lupost_normal_mcmc(1, c(0, 1))

normal_mcmc <- function(true_params, prior_params) {
  
  y <- 
  
}