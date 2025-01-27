############################################ Normal Model ##########################################

synthetic_normal <- function(num_particles, particles, parameters, mean = F, unknown = "all") {
  
  x.true <- parameters$x
  d_y <- length(x.true)
  
  likelihood_samples <- matrix(nrow = num_particles, ncol = d_y)
  
  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    current_params = get_current_normal_params(particle, particles, parameters, unknown = unknown)
    if (mean) {
      likelihood_samples[particle, ] <- current_params$alpha*current_params$x
    } else {
      likelihood_samples[particle, ] <- likelihood_normal(current_params)
    }
  }
  
  return(likelihood_samples)
  
}

synthetic_lognormal <- function(num_particles, particles, parameters, unknown = "all") {
  
  x.true <- parameters$x
  d_y <- length(x.true)
  
  likelihood_samples <- matrix(nrow = num_particles, ncol = d_y)
  
  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    current_params = get_current_normal_params(particle, particles, parameters, unknown = unknown)
    likelihood_samples[particle, ] <- likelihood_lognormal(current_params)
  }
  
  return(likelihood_samples)
  
}
############################################ Malaria Model #########################################

synthetic_log_malaria <- function(num_particles, particles, parameters, mean = F, unknown = "all") {
  
  likelihood_samples <- matrix(nrow = num_particles, ncol = 129)
  
  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    
    # These parameters are in unconstrained form
    current_params <- get_current_malaria_params(particle, particles, parameters, unknown = unknown)
    
    # print(current_params)
    if (mean) {
      current_params <- constrain_malaria_params(current_params)
      sample <- log(diff(likelihood_malaria_mean(current_params)))
    }
    else {
      sample <- likelihood_malaria(current_params)
    }
    
    likelihood_samples[particle, ] <- sample
    
  }
  
  return(likelihood_samples)
  
}

synthetic_malaria <- function(num_particles, particles, parameters, unknown = "all") {
  
  likelihood_samples <- matrix(nrow = num_particles, ncol = 129)
  
  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    
    # These parameters are in unconstrained form
    current_params <- get_current_malaria_params(particle, particles, parameters, unknown = unknown)
    sample <- likelihood_malaria(current_params, log_obs = F)
    likelihood_samples[particle, ] <- sample
    
  }
  
  return(likelihood_samples)
  
}