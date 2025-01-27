initialise_malaria_particles_var_only <- function(num_particles, parameters) {
  
  # Sample from the prior distribution
  prior_samples <- matrix(nrow = num_particles, ncol = 1)
  prior_samples[, 1] <- log_sigma_prior_sample(parameters, num_particles)
  
  # Initialise the particles and likelihood draws
  particles <- prior_samples
  
  return(particles)
}

eki_malaria_var_only <- function(num_particles, true_data, true_params, prior_params, adaptive = F) {
  
  initial_particles <- initialise_malaria_particles_var_only(num_particles, prior_params)
  
  # True parameters are specified noise
  if (adaptive) {
    return(eki_adaptive(num_particles, initial_particles, true_data, true_params, synthetic_malaria_var_only, densities_malaria_var_only))
  }
  else {
    return(eki(num_particles, initial_particles, true_data, true_params, synthetic_malaria_var_only))
  }
  
}