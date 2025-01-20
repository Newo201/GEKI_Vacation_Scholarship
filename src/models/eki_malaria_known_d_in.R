densities_malaria_known_d_in <- function(true_data, num_particles, particles, parameters) {
  
  likelihood_densities <- rep(0, num_particles)
  
  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    # These parameters are in unconstrained form
    current_params <- list(d_in = parameters$d_in,
                           phi = particles[particle, 1],
                           eta0 = particles[particle, 2],
                           sigma = parameters$sigma)
    current_params <- constrain_malaria_params(current_params)
    # TODO: currently I am solving the diff equation twice: here and in sampling function
    # Need to figure out how to optimise the workflow to avoid this
    current_mean <- log(diff(likelihood_malaria_mean(current_params)))
    likelihood_densities[particle] <- loglike_malaria(true_data, current_mean, current_params)
    # print(likelihood_densities[particle])
  }
  
  return(likelihood_densities)
  
}

synthetic_malaria_known_d_in <- function(num_particles, particles, parameters, mean = F) {
  
  likelihood_samples <- matrix(nrow = num_particles, ncol = 129)
  
  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    
    # These parameters are in unconstrained form
    current_params <- list(d_in = parameters$d_in,
                           phi = particles[particle, 1],
                           eta0 = particles[particle, 2],
                           sigma = parameters$sigma)
    
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

initialise_malaria_particles_known_d_in <- function(num_particles, parameters) {
  
  # Sample from the prior distribution
  prior_samples <- matrix(nrow = num_particles, ncol = 3)
  prior_samples[, 1] <- logit_phi_prior_sample(parameters, num_particles)
  prior_samples[, 2] <- logit_eta0_prior_sample(parameters, num_particles)
  
  # Initialise the particles and likelihood draws
  particles <- prior_samples
  
  return(particles)
}

eki_malaria_known_d_in <- function(num_particles, true_data, true_params, prior_params, adaptive = F, general = T) {
  
  initial_particles <- initialise_malaria_particles_known_d_in(num_particles, prior_params)
  
  # True parameters are specified noise
  if (adaptive && general) {
    return(eki_adaptive(num_particles, initial_particles, true_data, true_params, synthetic_malaria_known_d_in, densities_malaria_known_d_in))
  }
  else if (adaptive) {
    return(eki_adaptive_known_noise(num_particles, initial_particles, true_data, true_params, 
                                    synthetic_malaria_known_d_in, densities_malaria_known_d_in))
  } 
  else if (general) {
    return(eki(num_particles, initial_particles, true_data, true_params, synthetic_malaria_known_d_in))
  }
  else {
    return(eki_known_noise(num_particles, initial_particles, true_data, true_params, synthetic_malaria_known_d_in))
  }
  
}