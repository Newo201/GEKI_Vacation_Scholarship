pacman::p_load(pacman, mvtnorm, purrr)

source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/eki.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/samples/samples_malaria.R')
# source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/pdfs_normal.R')

# densities_normal <- function(true_data, num_particles, particles, parameters) {
#   
#   x.true <- parameters$x
#   d_y <- length(x.true)
#   
#   likelihood_densities <- rep(0, num_particles)
#   
#   # For each particle, draw one observation from the likelihood
#   # ToDo: vectorise this operation
#   for (particle in 1:num_particles) {
#     current_params = list(alpha = particles[particle, 1], x = x.true, sigma = sqrt(exp(particles[particle, 2])))
#     likelihood_densities[particle] <- loglike_pdf(true_data, current_params)
#   }
#   
#   return(likelihood_densities)
#   
# }

synthetic_malaria <- function(num_particles, particles, parameters) {
  
  x.true <- parameters$x
  d_y <- length(x.true)
  
  likelihood_samples <- matrix(nrow = num_particles, ncol = d_y)
  
  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    
    # These parameters are in unconstrained form
    current_params <- c(d_in = particles[particle, 1],
                        phi = particles[particle, 2],
                        eta0 = particles[particle, 3],
                        sigma = particles[particle, 4])
    
    likelihood_samples[particle, ] <- likelihood_malaria(current_params)
    
  }
  
  return(likelihood_samples)
  
}

initialise_malaria_particles <- function(num_particles, parameters) {
  
  # Sample from the prior distribution
  prior_samples <- matrix(nrow = num_particles, ncol = 2)
  prior_samples[, 1] <- log_din_prior_sample(parameters, num_particles)
  prior_samples[, 2] <- logit_phi_prior_sample(parameters, num_particles)
  prior_samples[, 3] <- logit_eta0_prior_sample(parameters, num_particles)
  prior_samples[, 4] <- log_sigma_prior_sample(parameters, num_particles)
  
  # Initialise the particles and likelihood draws
  particles <- prior_samples
  
  return(particles)
}

eki_malaria <- function(num_particles, true_params, prior_params, adaptive = F) {
  
  initial_particles <- initialise_malaria_particles(num_particles, prior_params)
  
  if (adaptive) {
    return(eki_adaptive(num_particles, initial_particles, true_params, 
                        likelihood_malaria, synthetic_malaria, densities_malaria))
  }
  else {
    return(eki(num_particles, initial_particles, true_params, 
               likelihood_malaria, synthetic_malaria, densities_malaria))
  }
  
}