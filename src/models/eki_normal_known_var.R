source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/models/synthetic_normal.R')

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
                        synthetic_mean_normal, synthetic_data_normal_known_var,
                        densities_normal_known_var))
  }
  else if (adaptive) {
    return(eki_adaptive_known_noise(num_particles, initial_particles, true_data, true_params,
                                    synthetic_mean_normal, densities_normal_known_var))
  }
  else if (general) {
    return(eki(num_particles, initial_particles, true_data, 
               true_params, synthetic_mean_normal, synthetic_data_normal_known_var))
  }
  else {
    return(eki_known_noise(num_particles, initial_particles, true_data, true_params,
                           synthetic_mean_normal))
  }
}





