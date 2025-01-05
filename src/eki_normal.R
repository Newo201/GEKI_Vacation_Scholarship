pacman::p_load(pacman, mvtnorm)

source('src/pdfs_normal.R')
source('src/samples_normal.R')


########################## EKI Algorithm ###################################

eki_normal <- function(iterations, parameters) {
  
  x.true <- parameters$x
  alpha.sd <- parameters$alpha.sd
  sigma2.sd <- parameters$sigma2.sd
  
  d_y <- length(x.true)
  
  # Simulate data using true parameters
  simulated_data <- generate_data(iterations, parameters)
  
  # Sample from the prior distribution
  prior_samples <- matrix(nrow = iterations, ncol = 2)
  prior_samples[, 1] <- alpha_prior_sample(alpha.sd, iterations)
  # Question: do I have to move particles in log(sigma2) space or in sigma2 space
  # currently in log(sigma2) space
  prior_samples[, 2] <- logsigma2_prior_sample(sigma2.sd, iterations)
  
  particles <- prior_samples
  
  # Until we reach a temperature of one do the following
  for (temp in 1:10) {
    
    # Calculate the current temperature
    
    # Sample from the likelihood
    
    likelihood_samples <- matrix(nrow = iterations, ncol = d_y)
    for (i in 1:iterations) {
      likelihood_samples[i, ] <- likelihood_sample(particles[i, 1], x.true, exp(particles[i, 2]))
    }
    
    # Calculate the covariance matrices
    C_xx = cov(particles)
    # print(C_xx)
    C_yy = cov(likelihood_samples)
    C_xy = cov(particles, likelihood_samples)
    C_yx = cov(likelihood_samples, particles)
    
    C_y_given_x = C_yy - C_yx %*% solve(C_xx) %*% C_xy
    
    # Generate perturbations
    eta <- matrix(data = rmvnorm(n = iterations, mean = 0, sigma = 9*C_y_given_x), nrow = iterations, ncol = d_y)
    
    # print(C_xy)
    # print(solve((C_yy + 10*C_y_given_x)))
    # 
    # print(simulated_data - likelihood_samples - eta)
    # print(C_xy %*% solve((C_yy + 10*C_y_given_x)) %*% t((simulated_data - likelihood_samples - eta)))
    
    # Move the particles
    # Todo: fix up and check the matrix dimensions
    particles <- particles + t(C_xy %*% solve((C_yy + 9*C_y_given_x)) %*% t((simulated_data - likelihood_samples - eta)))
    

  }
  
  return(particles)
}

parameters <- list(alpha = 2, sigma2 = 5, x = 5, alpha.sd = 5, sigma2.sd = 2)
result = eki_normal(100, parameters)

hist(result[, 1])
hist(exp(result[, 2]))

alpha_test <- alpha_prior_sample(5, 1000)
hist(alpha_test)

cov(alpha_test, alpha_test)


