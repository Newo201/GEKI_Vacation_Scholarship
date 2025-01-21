initialise_malaria_particles_known_var <- function(num_particles, parameters) {
  
  # Sample from the prior distribution
  prior_samples <- matrix(nrow = num_particles, ncol = 3)
  prior_samples[, 1] <- log_din_prior_sample(parameters, num_particles)
  prior_samples[, 2] <- logit_phi_prior_sample(parameters, num_particles)
  prior_samples[, 3] <- logit_eta0_prior_sample(parameters, num_particles)
  
  # Initialise the particles and likelihood draws
  particles <- prior_samples
  
  return(particles)
}

eki_malaria_known_var <- function(num_particles, true_data, true_params, prior_params, adaptive = F, general = T) {
  
  initial_particles <- initialise_malaria_particles_known_var(num_particles, prior_params)
  
  # True parameters are specified noise
  if (adaptive && general) {
    return(eki_adaptive(num_particles, initial_particles, true_data, true_params, synthetic_mean_malaria, 
                        synthetic_data_malaria_known_var, densities_malaria_known_var))
  }
  else if (adaptive) {
    return(eki_adaptive_known_noise(num_particles, initial_particles, true_data, true_params, 
                                    synthetic_mean_malaria, 
                                    densities_malaria_known_var))
  } 
  else if (general) {
    return(eki(num_particles, initial_particles, true_data, true_params, 
               synthetic_mean_malaria, synthetic_data_malaria_known_var))
  }
  else {
    return(eki_known_noise(num_particles, initial_particles, true_data, true_params, 
                           synthetic_mean_malaria))
  }
  
}