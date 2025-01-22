########################## EKI Algorithm ####################################

eki <- function(num_particles, initial_particles, true_data, true_params, synthetic_mean_func, synthetic_data_func) {

  # Synthetic_data_func -> a function which draws samples from the likelihood using particle parameters
  ## Takes num_particles, particles and true_params as arguments
  
  d_y <- length(true_data)
  # I'm replicating this data for the number of particles to make the dimensions easier to work with
  simulated_data <- matrix(true_data, nrow = num_particles, ncol = d_y, byrow = T)
  
  # Initialise the particles and likelihood draws
  particles <- initial_particles
  
  # Until we reach a temperature of one do the following
  current_temp <- 0
  temp_sequence <- c(current_temp)
  
  # Until we reach a temperature of one do the following
  for (temp in 1:10) {
    
    likelihood_means <- synthetic_mean_func(num_particles, particles, true_params)
    likelihood_samples <- synthetic_data_func(num_particles, particles, likelihood_means, true_params)

    covariances <- calculate_covariances(particles, likelihood_samples)
    
    # ToDo:  Calculate the current temperature
    temp_difference = 1/10
    print(temp)
    particles <- update_particles(temp_difference, particles, simulated_data, likelihood_samples, covariances, num_particles)
    current_temp <- current_temp + 1/10
    temp_sequence <- c(temp_sequence, current_temp)
  }
  
  return(list(particles = particles, temp_sequence = temp_sequence))
}

################## EKI Algorithm with Adaptive Temperature ##############################

eki_adaptive <- function(num_particles, initial_particles, true_data, true_params, 
                         synthetic_mean_func, synthetic_data_func, density_func) {
  
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
    
    # print(particles)
    likelihood_means <- synthetic_mean_func(num_particles, particles, true_params)
    likelihood_samples <- synthetic_data_func(num_particles, particles, likelihood_means, true_params)
    print(likelihood_samples)
    ll_densities <- density_func(true_data, num_particles, particles, likelihood_means, true_params)
    # print(ll_densities)
    
    covariances <- calculate_covariances(particles, likelihood_samples)
    # print(covariances)
    
    # Find the next temperature
    next_temp <- find_next_temp(current_temp, ll_densities, num_particles*0.5)
    print(next_temp)
    temp_difference <- next_temp - current_temp
    particles <- update_particles(temp_difference, particles, simulated_data, likelihood_samples, covariances, num_particles)
    current_temp <- next_temp
    temp_sequence <- c(temp_sequence, current_temp)
    
  }
  
  return(list(particles = particles, temp_sequence = temp_sequence))
}