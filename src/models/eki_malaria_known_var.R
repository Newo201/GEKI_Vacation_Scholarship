eki_malaria_known_var <- function(num_particles, true_data, true_params, prior_params, adaptive = F, general = T) {
  
  initial_particles <- initialise_malaria_particles_known_var(num_particles, prior_params)
  
  # True parameters are specified noise
  if (adaptive && general) {
    return(eki_adaptive(num_particles, initial_particles, true_data, true_params, synthetic_malaria_known_var, densities_malaria_known_var))
  }
  else if (adaptive) {
    return(eki_adaptive_known_noise(num_particles, initial_particles, true_data, true_params, 
                                    synthetic_malaria_known_var, densities_malaria_known_var))
  } 
  else if (general) {
    return(eki(num_particles, initial_particles, true_data, true_params, synthetic_malaria_known_var))
  }
  else {
    return(eki_known_noise(num_particles, initial_particles, true_data, true_params, synthetic_malaria_known_var))
  }
  
}