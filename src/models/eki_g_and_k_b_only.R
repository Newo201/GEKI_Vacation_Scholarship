pacman::p_load(pacman, mvtnorm, purrr)

source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/eki.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/samples/samples_g_and_k.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/pdfs/pdfs_g_and_k.R')

densities_g_and_k_b_only <- function(true_data, num_particles, particles, parameters) {

  likelihood_densities <- rep(0, num_particles)

  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    current_params = list(a = parameters$a, b = particles[particle, 1], 
                          g = parameters$g, k = parameters$k)
    likelihood_densities[particle] <- loglike_g_and_k(true_data, current_params)
  }

  return(likelihood_densities)

}

synthetic_g_and_k_b_only <- function(num_particles, particles, parameters) {
  
  likelihood_samples <- matrix(nrow = num_particles, ncol = 100)
  
  
  
  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    current_params = list(a = parameters$a, b = particles[particle, 1], 
                          g = parameters$g, k = parameters$k)
    likelihood_samples[particle, ] <- likelihood_g_and_k(current_params)
  }
  
  return(likelihood_samples)
  
}

synthetic_g_and_k_b_only_summary <- function(num_particles, particles, parameters) {
  
  likelihood_samples <- matrix(nrow = num_particles, ncol = 100)
  
  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    current_params = list(a = parameters$a, b = particles[particle, 1], 
                          g = parameters$g, k = parameters$k)
    likelihood_samples[particle, ] <- likelihood_g_and_k_summary(current_params)
  }
  
  return(likelihood_samples)
  
}

initialise_g_and_k_particles_b_only <- function(num_particles, parameters) {
  
  # Sample from the prior distribution
  prior_samples <- matrix(nrow = num_particles, ncol = 1)
  prior_samples[, 1] <- b_prior_sample(parameters, num_particles)
  
  # Initialise the particles
  particles <- prior_samples
  
  return(particles)
}

eki_g_and_k_b_only <- function(num_particles, true_params, prior_params, adaptive = F) {
  
  set.seed(2025)
  true_data <- likelihood_g_and_k_summary(true_params)
  set.seed(NULL)
  
  initial_particles <- initialise_g_and_k_particles_b_only(num_particles, prior_params)
  
  if (adaptive) {
    return(eki_adaptive(num_particles, initial_particles, true_data, synthetic_g_and_k_b_only_summary, 
                        densities_g_and_k_b_only))
  }
  else {
    return(eki(num_particles, initial_particles, true_data, synthetic_g_and_k_b_only_summary))
  }
  
}