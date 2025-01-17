pacman::p_load(pacman, mvtnorm, purrr)

source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/eki.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/samples/samples_g_and_k.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/pdfs/pdfs_g_and_k.R')

densities_g_and_k <- function(true_data, num_particles, particles, parameters) {

  likelihood_densities <- rep(0, num_particles)

  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    current_params = list(a = particles[particle, 1], b = particles[particle, 2], 
                          g = particles[particle, 3], k = particles[particle, 4])
    likelihood_densities[particle] <- loglike_g_and_k(true_data, current_params)
  }

  return(likelihood_densities)

}

synthetic_g_and_k <- function(num_particles, particles, parameters) {
  
  likelihood_samples <- matrix(nrow = num_particles, ncol = 100)
  
  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    current_params = list(a = particles[particle, 1], b = particles[particle, 2], 
                          g = particles[particle, 3], k = particles[particle, 4])
    likelihood_samples[particle, ] <- likelihood_g_and_k(current_params)
  }
  
  return(likelihood_samples)
  
}

synthetic_g_and_k_summary <- function(num_particles, particles, parameters) {
  
  likelihood_samples <- matrix(nrow = num_particles, ncol = 100)
  
  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    current_params = list(a = particles[particle, 1], b = particles[particle, 2], 
                          g = particles[particle, 3], k = particles[particle, 4])
    likelihood_samples[particle, ] <- likelihood_g_and_k_summary(current_params)
  }
  
  return(likelihood_samples)
  
}

initialise_g_and_k_particles <- function(num_particles, parameters) {
  
  # Sample from the prior distribution
  prior_samples <- matrix(nrow = num_particles, ncol = 4)
  prior_samples[, 1] <- a_prior_sample(parameters, num_particles)
  prior_samples[, 2] <- b_prior_sample(parameters, num_particles)
  prior_samples[, 3] <- g_prior_sample(parameters, num_particles)
  prior_samples[, 4] <- k_prior_sample(parameters, num_particles)
  
  # Initialise the particles
  particles <- prior_samples
  
  return(particles)
}

eki_g_and_k <- function(num_particles, true_params, prior_params, adaptive = F) {
  
  set.seed(2025)
  true_data <- likelihood_g_and_k_summary(true_params)
  set.seed(NULL)
  
  initial_particles <- initialise_g_and_k_particles(num_particles, prior_params)
  
  if (adaptive) {
    return(eki_adaptive(num_particles, initial_particles, true_data, synthetic_g_and_k_summary, 
                        densities_g_and_k))
  }
  else {
    return(eki(num_particles, initial_particles, true_data, synthetic_g_and_k_summary))
  }
  
}