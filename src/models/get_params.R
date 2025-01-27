get_current_malaria_params <- function(particle, particles, parameters, unknown = "all") {
  
  if (unknown == "all") {
    current_params <- list(d_in = particles[particle, 1],
                           phi = particles[particle, 2],
                           eta0 = particles[particle, 3],
                           sigma = particles[particle, 4])
  } else if (unknown == "var") {
    current_params <- list(d_in = parameters$d_in,
                           phi = parameters$phi,
                           eta0 = parameters$eta0,
                           sigma = particles[particle, 1])
  } else if (unknown == "mean") {
    current_params <- list(d_in = particles[particle, 1],
                           phi = particles[particle, 2],
                           eta0 = particles[particle, 3],
                           sigma = parameters$sigma)
  } else if (unknown == "d_in") {
    current_params <- list(d_in = particles[particle, 1],
                           phi = parameters$phi,
                           eta0 = parameters$eta0,
                           sigma = parameters$sigma)
  } else {
    current_params <- list(d_in = parameters$d_in,
                           phi = particles[particle, 1],
                           eta0 = particles[particle, 2],
                           sigma = parameters$sigma)
  }
  return(current_params)
  
}

get_current_normal_params <- function(particle, particles, parameters, unknown = "all") {
  
  if (unknown == "all") {
    current_params = list(alpha = particles[particle, 1], 
                          x = parameters$x, 
                          sigma = sqrt(exp(particles[particle, 2])))
  } else if (unknown == "mean") {
    current_params = list(alpha = particles[particle, 1], 
                          x = parameters$x, 
                          sigma = parameters$sigma)
  } else {
    current_params = list(alpha = parameters$alpha, 
                          x = parameters$x, 
                          sigma = sqrt(exp(particles[particle, 1])))
  }
  
  
}