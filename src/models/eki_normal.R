pacman::p_load(pacman, mvtnorm, purrr)

source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/eki.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/samples_normal.R')

synthetic_normal <- function(num_particles, particles, parameters) {
  
  x.true <- parameters$x
  d_y <- length(x.true)
  
  likelihood_samples <- matrix(nrow = num_particles, ncol = d_y)
  
  # For each particle, draw one observation from the likelihood
  # ToDo: vectorise this operation
  for (particle in 1:num_particles) {
    current_params = list(alpha = particles[particle, 1], x = x.true, sigma = sqrt(exp(particles[particle, 2])))
    likelihood_samples[particle, ] <- likelihood_normal(current_params)
  }
  
  return(likelihood_samples)

}

initialise_normal_particles <- function(num_particles, parameters) {
  
  x.true <- parameters$x
  d_y <- length(x.true)
  alpha.sd <- parameters$alpha.sd
  sigma2.sd <- parameters$sigma2.sd
  
  # We make a single draw from the likelihood using the true (unknown parameters)
  # I'm replicating this data for the number of particles to make the dimensions easier to work with
  simulated_data <- matrix(likelihood_normal(parameters), nrow = num_particles, ncol = d_y, byrow = T)
  
  # Sample from the prior distribution
  prior_samples <- matrix(nrow = num_particles, ncol = 2)
  prior_samples[, 1] <- alpha_prior_sample(alpha.sd, num_particles)
  prior_samples[, 2] <- logsigma2_prior_sample(sigma2.sd, num_particles)
  
  # Initialise the particles and likelihood draws
  particles <- prior_samples
  
  return(particles)
}

eki_normal <- function(num_particles, true_params, adaptive = F) {
  
  initial_particles <- initialise_normal_particles(num_particles, true_params)
  
  if (adaptive) {
    return(eki_adaptive(num_particles, initial_particles, true_params, likelihood_normal, synthetic_normal))
  }
  else {
    return(eki(num_particles, initial_particles, true_params, likelihood_normal, synthetic_normal))
  }
  
}