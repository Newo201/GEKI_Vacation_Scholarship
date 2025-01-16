source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/models/synthetic_malaria.R')

initialise_malaria_particles <- function(num_particles, parameters) {
  
  # Sample from the prior distribution
  prior_samples <- matrix(nrow = num_particles, ncol = 4)
  prior_samples[, 1] <- log_din_prior_sample(parameters, num_particles)
  prior_samples[, 2] <- logit_phi_prior_sample(parameters, num_particles)
  prior_samples[, 3] <- logit_eta0_prior_sample(parameters, num_particles)
  prior_samples[, 4] <- log_sigma_prior_sample(parameters, num_particles)
  
  # Initialise the particles and likelihood draws
  particles <- prior_samples
  
  return(particles)
}

eki_malaria <- function(num_particles, true_data, true_params, prior_params, adaptive = F) {
  
  initial_particles <- initialise_malaria_particles(num_particles, prior_params)
  
  # True parameters are specified noise
  if (adaptive) {
    return(eki_adaptive(num_particles, initial_particles, true_data, true_params, synthetic_malaria, densities_malaria))
  }
  else {
    return(eki(num_particles, initial_particles, true_data, true_params, synthetic_malaria))
  }
  
}