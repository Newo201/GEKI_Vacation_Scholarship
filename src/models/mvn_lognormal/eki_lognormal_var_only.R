densities_lognormal_var_only <- function(true_data, num_particles, particles, parameters) {
  
  x.true <- parameters$x
  d_y <- length(x.true)
  
  likelihood_densities <- rep(0, num_particles)
  
  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    current_params = list(alpha = parameters$alpha, x = x.true, sigma = sqrt(exp(particles[particle, 1])))
    likelihood_densities[particle] <- loglike_lognormal(true_data, current_params)
  }
  
  return(likelihood_densities)
  
}

synthetic_lognormal_var_only <- function(num_particles, particles, parameters, mean = F) {
  
  x.true <- parameters$x
  d_y <- length(x.true)
  
  likelihood_samples <- matrix(nrow = num_particles, ncol = d_y)
  
  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    current_params = list(alpha = parameters$alpha, x = x.true, sigma = sqrt(exp(particles[particle, 1])))
    if (mean) {
      likelihood_samples[particle, ] <- current_params$alpha*current_params$x
    } else {
      likelihood_samples[particle, ] <- likelihood_lognormal(current_params)
    }
  }
  
  return(likelihood_samples)
  
}

initialise_lognormal_particles_var_only <- function(num_particles, parameters) {
  
  x.true <- parameters$x
  d_y <- length(x.true)
  
  # Sample from the prior distribution
  prior_samples <- matrix(nrow = num_particles, ncol = 1)
  prior_samples[, 1] <- logsigma2_prior_sample(parameters, num_particles)
  
  # Initialise the particles and likelihood draws
  particles <- prior_samples
  
  return(particles)
}

eki_lognormal_var_only <- function(num_particles, true_data, true_params, prior_params, adaptive = F) {
  
  initial_particles <- initialise_lognormal_particles_var_only(num_particles, prior_params)
  
  if (adaptive) {
    return(eki_adaptive(num_particles, initial_particles, true_data, 
                        true_params, synthetic_lognormal_var_only, densities_lognormal_var_only))
  }
  else {
    return(eki(num_particles, initial_particles, true_data, 
               true_params, synthetic_lognormal_var_only))
  }
  
}