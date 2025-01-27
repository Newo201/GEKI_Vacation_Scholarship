############################# Malaria Model ######################################

densities_log_malaria <- function(true_data, num_particles, particles, parameters, unknown = "all") {
  
  likelihood_densities <- rep(0, num_particles)
  
  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    current_params <- get_current_malaria_params(particle, particles, parameters, unknown = unknown)

    current_params <- constrain_malaria_params(current_params)
    # TODO: currently I am solving the diff equation twice: here and in sampling function
    # Need to figure out how to optimise the workflow to avoid this
    current_mean <- log(diff(likelihood_malaria_mean(current_params)))
    likelihood_densities[particle] <- loglike_malaria(true_data, current_mean, current_params)
  }
  
  return(likelihood_densities)
  
}

densities_malaria <- function(true_data, num_particles, particles, parameters, unknown = "all") {
  
  likelihood_densities <- rep(0, num_particles)
  
  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    current_params <- get_current_malaria_params(particle, particles, parameters, unknown = unknown)
    
    current_params <- constrain_malaria_params(current_params)
    # TODO: currently I am solving the diff equation twice: here and in sampling function
    # Need to figure out how to optimise the workflow to avoid this
    current_mean <- log(diff(likelihood_malaria_mean(current_params)))
    likelihood_densities[particle] <- loglike_malaria(true_data, current_mean, current_params, log_obs = F)
  }
  
  return(likelihood_densities)
  
}

##################################### Normal Model ###############################################

densities_normal <- function(true_data, num_particles, particles, parameters) {
  
  likelihood_densities <- rep(0, num_particles)
  
  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    current_params = get_current_normal_params(particle, particles, parameters, unknown = unknown)
    likelihood_densities[particle] <- loglike_pdf(true_data, current_params)
  }
  
  return(likelihood_densities)
  
}

densities_normal <- function(true_data, num_particles, particles, parameters) {
  
  likelihood_densities <- rep(0, num_particles)
  
  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    current_params = get_current_normal_params(particle, particles, parameters, unknown = unknown)
    likelihood_densities[particle] <- loglike_lognormal(true_data, current_params)
  }
  
  return(likelihood_densities)
  
}

