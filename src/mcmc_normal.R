pacman::p_load(pacman, mvtnorm, purrr, mcmc)
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/pdfs_normal.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/samples_normal.R')

lupost_normal <- function(y, theta.init, x, params) {
  
  alpha.sd <- params$alpha.sd
  sigma2.sd <- params$sigma2.sd
  
  alpha.init <- theta.init[1]
  logsigma2.init <- theta.init[2]
  
  like_params <- list(alpha = alpha.init, x = x, sigma = sqrt(exp(logsigma2.init)))
  
  return(alpha_logprior_pdf(alpha.init, alpha.sd) + 
           logsigma2_logprior_pdf(logsigma2.init, sigma2.sd) + 
           sum(loglike_pdf(y, like_params)))
}

normal_mcmc <- function(num_samples, true_params, prior_params) {
  
  # 1. Generate data given true parameters
  simulated_data <- likelihood_sample(true_params, num_samples)
  
  # 2. Create a partial function fixing data and parameters
  lupost_normal_mcmc <- partial(lupost_normal, y = simulated_data, 
                                x = true_params$x, params = prior_params)
  
  # 3. Draw parameters from prior distribution
  alpha.init <- alpha_prior_sample(prior_params$alpha.sd, 1)
  logsigma2.init <- logsigma2_prior_sample(prior_params$sigma2.sd, 1)
  theta.init <- c(alpha.init, logsigma2.init)
  
  # 4. Call on MCMC function
  chain <- metrop(lupost_normal_mcmc, theta.init, 1e3)
  return(chain)
  
}