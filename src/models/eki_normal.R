# densities_normal <- function(true_data, num_particles, particles, parameters) {
#   
#   x.true <- parameters$x
#   d_y <- length(x.true)
#   
#   likelihood_densities <- rep(0, num_particles)
#   
#   # For each particle, draw one observation from the likelihood
#   # ToDo: vectorise this operation
#   for (particle in 1:num_particles) {
#     current_params = list(alpha = particles[particle, 1], x = x.true, sigma = sqrt(exp(particles[particle, 2])))
#     likelihood_densities[particle] <- loglike_pdf(true_data, current_params)
#   }
#   
#   return(likelihood_densities)
#   
# }

synthetic_normal <- function(num_particles, particles, parameters, mean = F) {
  
  x.true <- parameters$x
  d_y <- length(x.true) + 1
  
  likelihood_samples <- matrix(nrow = num_particles, ncol = d_y)
  
  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    current_params = list(alpha = particles[particle, 1], x = x.true, sigma = sqrt(exp(particles[particle, 2])))
    if (mean) {
      likelihood_samples[particle, ] <- current_params$alpha*current_params$x
    } else {
      sample <- likelihood_normal(current_params)
      # Summarise into sufficient statistics
      likelihood_samples[particle, ] <- c(sample, sd(sample))
    }
  }
  
  return(likelihood_samples)

}

initialise_normal_particles <- function(num_particles, parameters) {
  
  x.true <- parameters$x
  d_y <- length(x.true)
  
  # Sample from the prior distribution
  prior_samples <- matrix(nrow = num_particles, ncol = 2)
  prior_samples[, 1] <- alpha_prior_sample(parameters, num_particles)
  prior_samples[, 2] <- logsigma2_prior_sample(parameters, num_particles)
  
  # Initialise the particles and likelihood draws
  particles <- prior_samples
  
  return(particles)
}

eki_normal <- function(num_particles, true_data, true_params, prior_params, adaptive = F) {
  
  initial_particles <- initialise_normal_particles(num_particles, prior_params)
  
  if (adaptive) {
    return(eki_adaptive(num_particles, initial_particles, true_data, 
                        true_params, synthetic_normal))
  }
  else {
    return(eki(num_particles, initial_particles, true_data, 
               true_params, synthetic_normal))
  }
  
}