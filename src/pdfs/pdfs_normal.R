# ToDo: make sure log likelihood can handle y containing more than one observation in the correct way
loglike_pdf <- function(y, parameters) {
  alpha = parameters$alpha
  sigma2 = parameters$sigma**2
  x = parameters$x
  # print(dim(y))
  d_y = length(x)
  return(dmvnorm(y, mean = alpha*x, sigma = sigma2*diag(d_y), log = T))
}

alpha_logprior_pdf <- function(alpha, prior_params) {
  alpha.mean = prior_params$alpha.mean
  alpha.sd = prior_params$alpha.sd
  
  return(dnorm(alpha, mean = alpha.mean, sd = alpha.sd, log = T))
}

logsigma2_logprior_pdf <- function(logsigma2, prior_params) {
  sigma2.mean = prior_params$sigma2.mean
  sigma2.sd  = prior_params$sigma2.sd
  
  return(dnorm(logsigma2, mean = sigma2.mean, sd = sigma2.sd))
}

# test <- loglike_pdf(matrix(data = c(1, 1, 1, 1, 1, 1, 1, 1, 1), nrow = 3, ncol = 3), 1, c(1,0,1), 1)