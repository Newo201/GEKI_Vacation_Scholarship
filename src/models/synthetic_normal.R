densities_normal <- function(true_data, num_particles, particles, likelihood_means, parameters, known_var = F) {
  
  likelihood_densities <- rep(0, num_particles)
  
  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    if (known_var) {
      current_params = list(sigma = parameters$sigma)
    } else {
      current_params = list(sigma = sqrt(exp(particles[particle, 2])))
    }
    current_mean = likelihood_means[particle, ]
    likelihood_densities[particle] <- loglike_pdf(true_data, current_mean, current_params)
  }
  
  return(likelihood_densities)

}

densities_normal_known_var <- partial(densities_normal, known_var = T)

synthetic_data_normal <- function(num_particles, particles, likelihood_means, parameters, known_var = F) {
  
  x.true <- parameters$x
  d_y <- length(x.true)
  
  likelihood_samples <- matrix(nrow = num_particles, ncol = d_y)
  
  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    if (known_var) {
      current_params = list(sigma = parameters$sigma)
    } else {
      current_params = list(sigma = sqrt(exp(particles[particle, 2])))
    }
    likelihood_samples[particle, ] <- likelihood_normal(likelihood_means, current_params)
  }
  
  return(likelihood_samples)
  
}

synthetic_data_normal_known_var <- partial(synthetic_data_normal, known_var = T)

synthetic_mean_normal <- function(num_particles, particles, parameters, known_mean = F) {
  
  x.true <- parameters$x
  if (known_mean) {
    return(matrix(parameters$alpha * x.true, nrow = num_particles, ncol = length(x.true)))
  } else {
    return(particles[, 1] %*% x.true)
  }

}

synthetic_mean_normal_known_mean <- partial(synthetic_mean_normal, known_mean = T)