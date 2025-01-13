pacman::p_load(pacman, mvtnorm, purrr)

source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/eki.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/samples/samples_g_and_k.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/pdfs/pdfs_g_and_k.R')

densities_g_and_k <- function(true_data, num_particles, particles, parameters) {

  x.true <- parameters$x
  d_y <- length(x.true)

  likelihood_densities <- rep(0, num_particles)

  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    current_params = list(alpha = particles[particle, 1], x = x.true, sigma = sqrt(exp(particles[particle, 2])))
    likelihood_densities[particle] <- loglike_g_and_k(true_data, current_params)
  }

  return(likelihood_densities)

}

synthetic_g_and_k <- function(num_particles, particles, parameters) {
  
  x.true <- parameters$x
  d_y <- length(x.true)
  
  likelihood_samples <- matrix(nrow = num_particles, ncol = d_y)
  
  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    current_params = list(a = particles[particle, 1], b = particles[particle, 2], 
                          g = particles[particle, 3], k = particles[particle, 4])
    likelihood_samples[particle, ] <- likelihood_g_and_k(current_params)
  }
  
  return(likelihood_samples)
  
}

initialise_g_and_k_particles <- function(num_particles, parameters) {
  
  x.true <- parameters$x
  d_y <- length(x.true)
  
  # Sample from the prior distribution
  prior_samples <- matrix(nrow = num_particles, ncol = 2)
  prior_samples[, 1] <- a_prior_sample(parameters, num_particles)
  prior_samples[, 2] <- b_prior_sample(parameters, num_particles)
  prior_samples[, 3] <- g_prior_sample(parameters, num_particles)
  prior_samples[, 4] <- k_prior_sample(parameters, num_particles)
  
  # Initialise the particles
  particles <- prior_samples
  
  return(particles)
}

eki_normal <- function(num_particles, true_params, prior_params, adaptive = F) {
  
  initial_particles <- initialise_normal_particles(num_particles, prior_params)
  
  if (adaptive) {
    return(eki_adaptive(num_particles, initial_particles, true_params, 
                        likelihood_g_and_k, synthetic_g_and_k, densities_g_and_k))
  }
  else {
    return(eki(num_particles, initial_particles, true_params, 
               likelihood_g_and_k, synthetic_g_and_k, densities_g_and_k))
  }
  
}