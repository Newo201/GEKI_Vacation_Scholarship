densities_malaria <- function(true_data, num_particles, particles, likelihood_means, 
                              parameters, known_var = F) {
  
  likelihood_densities <- rep(0, num_particles)
  # 
  # loglike_vectorised <- Vectorize(loglike_malaria, vectorize.args = c("likelihood_means", "sigma"))
  # 
  # if (known_var) {
  #   return(likelihood_densities <- loglike_vectorised(true_data, likelihood_means, exp(parameters$sigma)))
  # } else {
  #   return(likelihood_densities <- loglike_vectorised(true_data, likelihood_means, exp(particles[2, ])))
  # }
  # 
  
  for (particle in 1:num_particles) {
    
    current_mean = likelihood_means[particle, ]
    
    if (known_var) {
      current_sd = exp(parameters$sigma)
    } else {
      current_sd = exp(particles[particle, 2])
    }
    
    likelihood_densities[particle] <- loglike_malaria(true_data, current_mean, current_sd)
    
  }
  
  # # For each particle, draw one observation from the likelihood
  # # ToDo: vectorise this operation
  # for (particle in 1:num_particles) {
  #   if (d_in_only) {
  #     current_params <- list(d_in = particles[particle, 1],
  #                            phi = parameters$phi,
  #                            eta0 = parameters$eta0,
  #                            sigma = parameters$sigma)
  #   }
  #   else if (known_var) {
  #     current_params <- list(d_in = particles[particle, 1],
  #                            phi = particles[particle, 2],
  #                            eta0 = particles[particle, 3],
  #                            sigma = parameters$sigma)
  #   }
  #   else {
  #     # These parameters are in unconstrained form
  #     current_params <- list(d_in = particles[particle, 1],
  #                            phi = particles[particle, 2],
  #                            eta0 = particles[particle, 3],
  #                            sigma = particles[particle, 4])
  #   }
  #   current_params <- constrain_malaria_params(current_params)
  #   current_mean <- likelihood_means[particle, ]
  #   likelihood_densities[particle] <- loglike_malaria(true_data, current_mean, current_params)
  # }
  # 
  return(likelihood_densities)
  
}

densities_malaria_known_var <- partial(densities_malaria, known_var = T)

synthetic_mean_malaria <- function(num_particles, particles, parameters, d_in_only = F) {
  
  likelihood_means <- matrix(nrow = num_particles, ncol = 129)
  # mean_vectorised <- Vectorize(likelihood_malaria_mean, vectorize.args = c("d_in", "phi", "eta0"))
  # 
  # if (d_in_only) {
  #   
  #   d_in_particles <- particles[, 1]
  #   likelihood_means <- mean_vectorised(d_in_particles, parameters$phi, parameters$eta0)
  # }
  # else {
  #   d_in_particles <- particles[, 1]
  #   phi_particles <- particles[, 2]
  #   eta0_particles <- particles[, 3]
  #   likelihood_means <- mean_vectorised(d_in_particles, phi_particles, eta0_particles)
  # }
  
  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {

    if (d_in_only) {
      # These parameters are in unconstrained form
      current_params <- list(d_in = particles[particle, 1],
                             phi = parameters$phi,
                             eta0 = parameters$eta0,
                             sigma = 0)
    } else {
      # These parameters are in unconstrained form
      # Don't need sigma here just adding it so the code doesn't throw an error
      current_params <- list(d_in = particles[particle, 1],
                             phi = particles[particle, 2],
                             eta0 = particles[particle, 3],
                             sigma = 0)
    }


    # print(current_params)
    current_constrained_params <- constrain_malaria_params(current_params)
    means <- likelihood_malaria_mean(current_params$d_in, 
                                     current_params$phi,
                                     current_params$eta0)
    # print(means)
    # Convert to log difference
    likelihood_means[particle, ] <- log(diff(means))

  }

  return(likelihood_means)
}

synthetic_mean_malaria_d_in_only <- partial(synthetic_mean_malaria, d_in_only = T)

synthetic_data_malaria <- function(num_particles, particles, likelihood_means, 
                                   parameters, known_var = F, d_in_only = F) {
  
  likelihood_samples <- matrix(nrow = num_particles, ncol = 129)
  
  # print(num_particles)
  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    
    if (d_in_only) {
      current_params <- list(d_in = particles[particle, 1],
                             phi = parameters$phi,
                             eta0 = parameters$eta0,
                             sigma = parameters$sigma)
    }
    else if (known_var) {
      current_params <- list(d_in = particles[particle, 1],
                             phi = particles[particle, 2],
                             eta0 = particles[particle, 3],
                             sigma = parameters$sigma)
    }
    else {
      # These parameters are in unconstrained form
      current_params <- list(d_in = particles[particle, 1],
                             phi = particles[particle, 2],
                             eta0 = particles[particle, 3],
                             sigma = particles[particle, 4])
    }

    
    # print(current_params)
    current_constrained_params <- constrain_malaria_params(current_params)
    print(current_constrained_params)
    sample <- likelihood_malaria(likelihood_means[particle, ], current_constrained_params)
    likelihood_samples[particle, ] <- sample
    
  }
  
  return(likelihood_samples)
  
}

synthetic_data_malaria_known_var <- partial(synthetic_data_malaria, known_var = T)
synthetic_data_malaria_d_in_only <- partial(synthetic_data_malaria, d_in_only = T)