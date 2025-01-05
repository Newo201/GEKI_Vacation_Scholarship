# ToDo: make sure log likelihood can handle y containing more than one observation in the correct way
loglike_pdf <- function(y, alpha, x, sigma2) {
  # print(dim(y))
  d_y = length(x)
  return(dmvnorm(y, mean = alpha*x, sigma2*diag(d_y), log = T))
}

alpha_logprior_pdf <- function(alpha, alpha.sd) {
  return(dnorm(alpha, mean = 0, sd = alpha.sd, log = T))
}

logsigma2_logprior_pdf <- function(logsigma2, sigma.sd) {
  return(dnorm(logsigma2, mean = 0, sd = sigma.sd))
}

test <- loglike_pdf(matrix(data = c(1, 1, 1, 1, 1, 1, 1, 1, 1), nrow = 3, ncol = 3), 1, c(1,0,1), 1)