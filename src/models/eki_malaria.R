# densities_malaria <- function(true_data, num_particles, particles, parameters) {
# 
#   likelihood_densities <- rep(0, num_particles)
# 
#   # For each particle, draw one observation from the likelihood
#   # ToDo: vectorise this operation
#   for (particle in 1:num_particles) {
#     current_params <- list(d_in = particles[particle, 1],
#                            phi = particles[particle, 2],
#                            eta0 = particles[particle, 3],
#                            sigma = particles[particle, 4])
#     current_params <- constrain_malaria_params(current_params)
#     # TODO: currently I am solving the diff equation twice: here and in sampling function
#     # Need to figure out how to optimise the workflow to avoid this
#     current_mean <- log(diff(likelihood_malaria_mean(current_params)))
#     likelihood_densities[particle] <- loglike_malaria(true_data, current_mean, current_params, log_obs = F)
#   }
# 
#   return(likelihood_densities)
# 
# }

synthetic_malaria <- function(num_particles, particles, parameters) {
  
  likelihood_samples <- matrix(nrow = num_particles, ncol = 130)
  
  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    
    # These parameters are in unconstrained form
    current_params <- list(d_in = particles[particle, 1],
                        phi = particles[particle, 2],
                        eta0 = particles[particle, 3],
                        sigma = particles[particle, 4])
    
    # print(current_params)
    sample <- likelihood_malaria(current_params, log_obs = F)
    # print(sample)
    likelihood_samples[particle, ] <- c(sample, sd(sample))
    
  }
  
  return(likelihood_samples)
  
}

initialise_malaria_particles <- function(num_particles, parameters) {
  
  # Sample from the prior distribution
  prior_samples <- matrix(nrow = num_particles, ncol = 4)
  prior_samples[, 1] <- log_din_prior_sample(parameters, num_particles)
  prior_samples[, 2] <- logit_phi_prior_sample(parameters, num_particles)
  prior_samples[, 3] <- logit_eta0_prior_sample(parameters, num_particles)
  prior_samples[, 4] <- log_sigma_prior_sample(parameters, num_particles)
  
  # Initialise the particles and likelihood draws
  particles <- prior_samples
  
  return(particles)
}

eki_malaria <- function(num_particles, true_data, true_params, prior_params, adaptive = F) {
  
  initial_particles <- initialise_malaria_particles(num_particles, prior_params)
  
  # True parameters are specified noise
  if (adaptive) {
    return(eki_adaptive(num_particles, initial_particles, true_data, true_params, synthetic_malaria))
  }
  else {
    return(eki(num_particles, initial_particles, true_data, true_params, synthetic_malaria))
  }
  
}