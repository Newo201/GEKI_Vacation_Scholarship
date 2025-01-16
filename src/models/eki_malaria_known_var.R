densities_malaria_known_var <- function(true_data, num_particles, particles, parameters) {

  likelihood_densities <- rep(0, num_particles)

  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    current_params = list(sigma = log(parameters$sigma))
    likelihood_densities[particle] <- loglike_malaria(true_data, current_params)
  }

  return(likelihood_densities)

}

synthetic_malaria_known_var <- function(num_particles, particles, parameters, mean = F) {
  
  likelihood_samples <- matrix(nrow = num_particles, ncol = 129)
  
  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    
    # These parameters are in unconstrained form
    current_params <- list(d_in = particles[particle, 1],
                        phi = particles[particle, 2],
                        eta0 = particles[particle, 3],
                        sigma = log(parameters$sigma))
    
    # print(current_params)
    if (mean) {
      current_params <- constrain_malaria_params(current_params)
      sample <- log(diff(likelihood_malaria_mean(current_params)))
    }
    else {
      sample <- likelihood_malaria(current_params)
    }

    # print(sample)
    likelihood_samples[particle, ] <- sample
    
  }
  
  return(likelihood_samples)
  
}

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