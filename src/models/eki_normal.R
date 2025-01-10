pacman::p_load(pacman, mvtnorm, purrr)

source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/eki.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/samples_normal.R')

synthetic_normal <- function(true_data, num_particles, particles, parameters) {
  
  x.true <- parameters$x
  d_y <- length(x.true)
  
  likelihood_samples <- matrix(nrow = num_particles, ncol = d_y)
  likelihood_densities <- rep(0, num_particles)
  
  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    current_params = list(alpha = particles[particle, 1], x = x.true, sigma = sqrt(exp(particles[particle, 2])))
    likelihood_samples[particle, ] <- likelihood_normal(current_params)
    likelihood_densities[particle] <- loglike_pdf(true_data, current_params)
  }
  
  return(list(samples = likelihood_samples, ll_densities = likelihood_densities))

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

eki_normal <- function(num_particles, true_params, prior_params, adaptive = F) {
  
  initial_particles <- initialise_normal_particles(num_particles, prior_params)
  
  if (adaptive) {
    return(eki_adaptive(num_particles, initial_particles, true_params, likelihood_normal, synthetic_normal))
  }
  else {
    return(eki(num_particles, initial_particles, true_params, likelihood_normal, synthetic_normal))
  }
  
}