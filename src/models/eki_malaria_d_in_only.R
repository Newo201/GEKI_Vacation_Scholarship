source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/models/synthetic_malaria.R')

initialise_malaria_particles_d_in_only <- function(num_particles, parameters) {
  
  # Sample from the prior distribution
  prior_samples <- matrix(nrow = num_particles, ncol = 1)
  prior_samples[, 1] <- log_din_prior_sample(parameters, num_particles)
  
  # Initialise the particles and likelihood draws
  particles <- prior_samples
  
  return(particles)
}

eki_malaria_d_in_only <- function(num_particles, true_data, true_params, prior_params, adaptive = F, general = T) {
  
  initial_particles <- initialise_malaria_particles_d_in_only(num_particles, prior_params)
  
  # True parameters are specified noise
  if (adaptive && general) {
    return(eki_adaptive(num_particles, initial_particles, true_data, true_params, synthetic_malaria_d_in_only, densities_malaria_d_in_only))
  }
  else if (adaptive) {
    return(eki_adaptive_known_noise(num_particles, initial_particles, true_data, true_params, 
                                    synthetic_malaria_d_in_only, densities_malaria_d_in_only))
  } 
  else if (general) {
    return(eki(num_particles, initial_particles, true_data, true_params, synthetic_malaria_d_in_only))
  }
  else {
    return(eki_known_noise(num_particles, initial_particles, true_data, true_params, synthetic_malaria_d_in_only))
  }
  
}