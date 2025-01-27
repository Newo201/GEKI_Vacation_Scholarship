########################## EKI Algorithm ####################################

eki_known_noise <- function(num_particles, initial_particles, true_data, true_params, synthetic_data_func) {

  # Synthetic_data_func -> a function which calculates likelihood means from a given set of parameters
  ## Takes num_particles, particles and true_params as arguments
  known_noise = exp(true_params$sigma)
  
  print(known_noise)
  
  d_y <- length(true_data)
  # I'm replicating this data for the number of particles to make the dimensions easier to work with
  simulated_data <- matrix(true_data, nrow = num_particles, ncol = d_y, byrow = T)
  
  # Initialise the particles and likelihood draws
  particles <- initial_particles
  
  current_temp <- 0
  temp_sequence <- c(current_temp)
  
  # Until we reach a temperature of one do the following
  for (temp in 1:10) {
    
    # Generate synthetic means
    likelihood_means <- synthetic_data_func(num_particles, particles, true_params, mean = T)
    # Calculate empirical covariances
    covariances <- calculate_covariances_known_noise(particles, likelihood_means)
    
    # Get the next temperature
    print(glue("Next temp is {temp/10}"))
    temp_difference = 1/10
    # Update the particles
    particles <- update_particles_known_noise(temp_difference, particles, simulated_data, 
                                              likelihood_means, covariances, num_particles,
                                              known_noise)
    current_temp <- current_temp + 1/10
    temp_sequence <- c(temp_sequence, current_temp)
  }
  
  return(list(particles = particles, temp_sequence = temp_sequence))
}

################## EKI Algorithm with Adaptive Temperature ##############################
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/utils/eki_helper.R')

eki_adaptive_known_noise <- function(num_particles, initial_particles, true_data, true_params, 
                                     synthetic_data_func, density_func) {
  
  # Assumes params are in constrained form
  known_noise = exp(true_params$sigma)
  print(known_noise)
  
  # Synthetic_data_func -> a function which draws samples from the likelihood using particle parameters
  ## Takes true_params, particles and number of particles as arguments
  
  d_y <- length(true_data)
  # I'm replicating this data for the number of particles to make the dimensions easier to work with
  simulated_data <- matrix(true_data, nrow = num_particles, ncol = d_y, byrow = T)
  
  # Initialise the particles and likelihood draws
  particles <- initial_particles
  
  # Until we reach a temperature of one do the following
  current_temp <- 0
  temp_sequence <- c(current_temp)
  
  while (current_temp < 1) {
    
    likelihood_means <- synthetic_data_func(num_particles, particles, true_params, mean = T)
    ll_densities <- density_func(true_data, num_particles, particles, true_params)
    
    covariances <- calculate_covariances_known_noise(particles, likelihood_means)
    
    # Find the next temperature
    next_temp <- find_next_temp(current_temp, ll_densities, num_particles*0.5)
    print(glue("Next temp is {next_temp}"))
    temp_difference <- next_temp - current_temp
    particles <- update_particles_known_noise(temp_difference, particles, simulated_data, 
                                              likelihood_means, covariances, num_particles,
                                              known_noise)
    current_temp <- next_temp
    temp_sequence <- c(temp_sequence, current_temp)
    
  }
  
  return(list(particles = particles, temp_sequence = temp_sequence))
}