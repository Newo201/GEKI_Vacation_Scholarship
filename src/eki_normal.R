pacman::p_load(pacman, mvtnorm)

source('src/pdfs_normal.R')
source('src/samples_normal.R')

initialise_particles(num_particles, parameters) {
  
  x.true <- parameters$x
  d_y <- length(x.true)
  alpha.sd <- parameters$alpha.sd
  sigma2.sd <- parameters$sigma2.sd
  
  # We make a single draw from the likelihood using the true (unknown parameters)
  # I'm replicating this data for the number of particles to make the dimensions easier to work with
  simulated_data <- matrix(likelihood_sample(parameters, 1), nrow = num_particles, ncol = d_y, byrow = T)
  
  # Sample from the prior distribution
  prior_samples <- matrix(nrow = num_particles, ncol = 2)
  prior_samples[, 1] <- alpha_prior_sample(alpha.sd, num_particles)
  prior_samples[, 2] <- logsigma2_prior_sample(sigma2.sd, num_particles)
  
  # Initialise the particles and likelihood draws
  particles <- prior_samples
  
  return(particles)
}

update_particles(temp_difference, num_particles) {
  
  likelihood_samples <- matrix(nrow = num_particles, ncol = d_y)
  
  # For each particle, draw one observation from the likelihood
  for (particle in 1:num_particles) {
    current_params = list(alpha = particles[particle, 1], x = x.true, sigma2 = exp(particles[particle, 2]))
    likelihood_samples[particle, ] <- likelihood_sample(current_params, 1)
  }
  
  # Calculate the covariance matrices
  C_xx = cov(particles)
  # print(C_xx)
  C_yy = cov(likelihood_samples)
  C_xy = cov(particles, likelihood_samples)
  C_yx = cov(likelihood_samples, particles)
  
  C_y_given_x = C_yy - C_yx %*% solve(C_xx) %*% C_xy
  
  # Generate perturbations
  eta <- matrix(data = rmvnorm(n = num_particles, mean = rep(0, d_y), sigma = 9*C_y_given_x), nrow = num_particles, ncol = d_y)
  
  # Move the particles
  particles <- particles + t(C_xy %*% solve((C_yy + 9*C_y_given_x)) %*% t((simulated_data - likelihood_samples - eta)))
  return(particles)
}

########################## EKI Algorithm ####################################

eki_normal <- function(num_particles, parameters) {
  
  x.true <- parameters$x
  d_y <- length(x.true)
  
  # We make a single draw from the likelihood using the true (unknown parameters)
  # I'm replicating this data for the number of particles to make the dimensions easier to work with
  simulated_data <- matrix(likelihood_sample(parameters, 1), nrow = num_particles, ncol = d_y, byrow = T)
  
  # Initialise the particles and likelihood draws
  particles <- initialise_particles(num_particles, parameters)
  
  # Until we reach a temperature of one do the following
  for (temp in 1:10) {
    
    # ToDo:  Calculate the current temperature
    temp_difference = 1/10
    particles <- update_particles(temp_difference)
  }
  
  return(particles)
}

################## EKI Algorithm with Adaptive Temperature ##############################

eki_normal_adaptive <- function(num_particles, parameters) {
  
  x.true <- parameters$x
  d_y <- length(x.true)
  
  # We make a single draw from the likelihood using the true (unknown parameters)
  # I'm replicating this data for the number of particles to make the dimensions easier to work with
  simulated_data <- matrix(likelihood_sample(parameters, 1), nrow = num_particles, ncol = d_y, byrow = T)
  
  # Initialise the particles and likelihood draws
  particles <- initialise_particles(num_particles, parameters)
  
  # Until we reach a temperature of one do the following
  current_temp <- 0
  
  while (current_temp < 1) {
    
    # Find the next temperature
    next_temp <- find_next_temp(current_temp, summary_stat, likelihood_samples, num_particles, num_particles*0.5)
    temp_difference <- next_temp - current_temp
    particles <- update_particles(temp_difference)
  
  }
  
  return(particles)
}


num_particles <- 50
parameters <- list(alpha = 2, sigma = 5, x = c(1,1), alpha.sd = 5, sigma2.sd = 2)
result = eki_normal(100, parameters)
 


