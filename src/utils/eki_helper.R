calculate_covariances <- function(particles, likelihood_samples, correlation = F) {
  
  # Calculate the covariance matrices
  C_xx = cov(particles)

  # Correlation is used for generating diagnostics to better understand the
  # performance of the EKI algorithm
  if (correlation) {
    C_yy = cor(likelihood_samples)
    C_xy = cor(particles, likelihood_samples)
    C_yx = cor(likelihood_samples, particles)
  } else {
    C_yy = cov(likelihood_samples)
    C_xy = cov(particles, likelihood_samples)
    C_yx = cov(likelihood_samples, particles)
  }
  
  C_y_given_x = C_yy - C_yx %*% ginv(C_xx) %*% C_xy
  # Presolving the inverse so we don't have to recalculate it every time
  C_y_given_x_inv <- ginv(C_y_given_x)
  
  return(list(C_xx = C_xx, C_yy = C_yy, C_xy = C_xy, C_yx = C_yx, 
              C_y_given_x = C_y_given_x, C_y_given_x_inv = C_y_given_x_inv))
  
}

calculate_covariances_known_noise <- function(particles, likelihood_means) {
  
  C_xx = cov(particles)
  C_hh = cov(likelihood_means)
  C_xh = cov(particles, likelihood_means)
  C_hx = cov(likelihood_means, particles)
  
  return(list(C_xx = C_xx, C_hh = C_hh, C_xh = C_xh, C_hx = C_hx))
  
}

update_particles <- function(temp_difference, particles, simulated_data, likelihood_samples, covariances, num_particles) {
  
  C_yx <- covariances$C_yx
  C_yy <- covariances$C_yy
  C_y_given_x <- covariances$C_y_given_x
  
  d_y <- dim(likelihood_samples)[2]
  
  # Generate perturbations
  eta <- rmvnorm(n = num_particles, mean = rep(0, d_y), sigma = (1/temp_difference - 1)*C_y_given_x)
  
  # Move the particles
  particles <- particles + (simulated_data - likelihood_samples - eta) %*% ginv(C_yy + (1/temp_difference - 1)*C_y_given_x) %*% C_yx
  if (sum(is.infinite(particles)) > 0) {
    print(det(C_yy + (1/temp_difference - 1)*C_y_given_x))
    stop("Some of the particle values are infinite")
  }
  return(particles)
}

update_particles_known_noise <- function(temp_difference, particles, simulated_data, 
                                         likelihood_means, covariances, num_particles,
                                         known_noise) {
  
  C_hx <- covariances$C_hx
  C_hh <- covariances$C_hh
  
  d_y <- dim(likelihood_means)[2]
  R <- known_noise**2 * diag(d_y)
  
  # Generate perturbations
  eta <- rmvnorm(n = num_particles, mean = rep(0, d_y), sigma = (1/temp_difference)*R)
  # Move the particles
  particles <- particles + (simulated_data - likelihood_means - eta) %*% ginv(C_hh + (1/temp_difference)*R) %*% C_hx

  return(particles)
}