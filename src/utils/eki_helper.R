calculate_covariances <- function(particles, likelihood_samples) {
  
  # Calculate the covariance matrices
  C_xx = cov(particles)
  # print(C_xx)
  C_yy = cov(likelihood_samples)
  C_xy = cov(particles, likelihood_samples)
  C_yx = cov(likelihood_samples, particles)
  
  C_y_given_x = C_yy - C_yx %*% solve(C_xx) %*% C_xy
  # Presolving the inverse so we don't have to recalculate it every time
  C_y_given_x_inv <- solve(C_y_given_x)
  
  return(list(C_xx = C_xx, C_yy = C_yy, C_xy = C_xy, C_yx = C_yx, 
              C_y_given_x = C_y_given_x, C_y_given_x_inv = C_y_given_x_inv))
  
}

update_particles <- function(temp_difference, particles, simulated_data, likelihood_samples, covariances, num_particles) {
  
  C_yx <- covariances$C_yx
  C_yy <- covariances$C_yy
  C_y_given_x <- covariances$C_y_given_x
  
  d_y <- dim(likelihood_samples)[2]
  
  # Generate perturbations
  eta <- rmvnorm(n = num_particles, mean = rep(0, d_y), sigma = (1/temp_difference - 1)*C_y_given_x)
  
  # Move the particles
  # particles <- particles + t(C_xy %*% solve((C_yy + (1/temp_difference - 1)*C_y_given_x)) %*% t((simulated_data - likelihood_samples - eta)))
  # ToDo: make sure that adjusting the dimensions produces the same update
  particles <- particles + (simulated_data - likelihood_samples - eta) %*% solve((C_yy + (1/temp_difference - 1)*C_y_given_x)) %*% C_yx
  return(particles)
}