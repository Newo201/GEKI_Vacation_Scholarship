densities_normal_known_var <- function(true_data, num_particles, particles, parameters) {
  
  x.true <- parameters$x
  sigma.true <- parameters$sigma
  d_y <- length(x.true)
  
  likelihood_densities <- rep(0, num_particles)
  
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    current_params = list(alpha = particles[particle, 1], x = x.true, sigma = sigma.true)
    likelihood_densities[particle] <- loglike_pdf(true_data, current_params)
  }
  
  return(likelihood_densities)
}

synthetic_normal_known_var <- function(num_particles, particles, parameters, mean = F) {
  
  x.true <- parameters$x
  sigma.true <- parameters$sigma
  d_y <- length(x.true)
  
  likelihood_samples <- matrix(nrow = num_particles, ncol = d_y)
  
  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    current_params = list(alpha = particles[particle, 1], x = x.true, sigma = sigma.true)
    if (mean) {
      likelihood_samples[particle, ] <- current_params$alpha*current_params$x
    } else {
      likelihood_samples[particle, ] <- likelihood_normal(current_params)
    }

  }
  
  return(likelihood_samples)
  
}

initialise_normal_particles_known_var <- function(num_particles, parameters) {
  
  x.true <- parameters$x
  d_y <- length(x.true)
  
  # Sample from the prior distribution
  prior_samples <- matrix(nrow = num_particles, ncol = 1)
  prior_samples[, 1] <- alpha_prior_sample(parameters, num_particles)
  
  # Initialise the particles and likelihood draws
  particles <- prior_samples
  
  return(particles)
}

eki_normal_known_var <- function(num_particles, true_data, true_params, prior_params, adaptive = F, general = T) {
  
  initial_particles <- initialise_normal_particles_known_var(num_particles, prior_params)
  
  if (adaptive && general) {
    return(eki_adaptive(num_particles, initial_particles, true_data, true_params, 
                        synthetic_normal_known_var,
                        densities_normal_known_var))
  }
  else if (adaptive) {
    return(eki_adaptive_known_noise(num_particles, initial_particles, true_data, true_params,
                                    synthetic_normal_known_var, densities_normal_known_var))
  }
  else if (general) {
    return(eki(num_particles, initial_particles, true_data, 
               true_params, synthetic_normal_known_var))
  }
  else {
    return(eki_known_noise(num_particles, initial_particles, true_data, true_params,
                           synthetic_normal_known_var))
  }
}





