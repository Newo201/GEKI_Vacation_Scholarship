pacman::p_load(pacman, mvtnorm)

source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/pdfs_normal.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/samples_normal.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/utils/tempering.R')

generate_likelihood_samples <- function(num_particles, particles, parameters) {
  
  x.true <- parameters$x
  sigma.true <- parameters$sigma
  d_y <- length(x.true)
  
  likelihood_samples <- matrix(nrow = num_particles, ncol = d_y)
  
  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    current_params = list(alpha = particles[particle, 1], x = x.true, sigma = sigma.true)
    likelihood_samples[particle, ] <- likelihood_sample(current_params, 1)
  }
  
  return(likelihood_samples)
  
}

initialise_particles <- function(num_particles, parameters) {
  
  x.true <- parameters$x
  d_y <- length(x.true)
  alpha.sd <- parameters$alpha.sd
  sigma2.sd <- parameters$sigma2.sd
  
  # We make a single draw from the likelihood using the true (unknown parameters)
  # I'm replicating this data for the number of particles to make the dimensions easier to work with
  simulated_data <- matrix(likelihood_sample(parameters, 1), nrow = num_particles, ncol = d_y, byrow = T)
  
  # Sample from the prior distribution
  prior_samples <- matrix(nrow = num_particles, ncol = 1)
  prior_samples[, 1] <- alpha_prior_sample(alpha.sd, num_particles)
  
  # Initialise the particles and likelihood draws
  particles <- prior_samples
  
  return(particles)
}

update_particles <- function(temp_difference, particles, simulated_data, likelihood_samples, num_particles) {
  
  # Calculate the covariance matrices
  C_xx = cov(particles)
  # print(C_xx)
  C_yy = cov(likelihood_samples)
  C_xy = cov(particles, likelihood_samples)
  C_yx = cov(likelihood_samples, particles)
  
  d_y <- dim(likelihood_samples)[2]
  
  C_y_given_x = C_yy - C_yx %*% solve(C_xx) %*% C_xy
  
  # Generate perturbations
  eta <- rmvnorm(n = num_particles, mean = rep(0, d_y), sigma = (1/temp_difference - 1)*C_y_given_x)
  
  # Move the particles
  # particles <- particles + t(C_xy %*% solve((C_yy + (1/temp_difference - 1)*C_y_given_x)) %*% t((simulated_data - likelihood_samples - eta)))
  # ToDo: make sure that adjusting the dimensions produces the same update
  particles <- particles + (simulated_data - likelihood_samples - eta) %*% solve((C_yy + (1/temp_difference - 1)*C_y_given_x)) %*% C_yx
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
  current_temp <- 0
  temp_sequence <- c(current_temp)
  
  # Until we reach a temperature of one do the following
  for (temp in 1:10) {
    
    likelihood_samples <- generate_likelihood_samples(num_particles, particles, parameters)
    
    # ToDo:  Calculate the current temperature
    temp_difference = 1/10
    particles <- update_particles(temp_difference, particles, simulated_data, likelihood_samples, num_particles)
    current_temp <- current_temp + 1/10
    temp_sequence <- c(temp_sequence, current_temp)
  }
  
  return(list(particles = particles, temp_sequence = temp_sequence))
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
  temp_sequence <- c(current_temp)
  
  while (current_temp < 1) {
    
    likelihood_samples <- generate_likelihood_samples(num_particles, particles, parameters)
    
    # Find the next temperature
    next_temp <- find_next_temp(current_temp, simulated_data, likelihood_samples, num_particles, num_particles*0.5)
    # print(next_temp)
    temp_difference <- next_temp - current_temp
    particles <- update_particles(temp_difference, particles, simulated_data, likelihood_samples, num_particles)
    current_temp <- next_temp
    temp_sequence <- c(temp_sequence, current_temp)
    
  }
  
  return(list(particles = particles, temp_sequence = temp_sequence))
}


