pacman::p_load(pacman, purrr)

get_weights <- function(next_temp, current_temp, simulated_data, likelihood_samples, covariances, num_particles) {
  
  C_y_given_x_inv <- covariances$C_y_given_x_inv
  
  change_in_temp <- next_temp - current_temp
  weights <- rep(0, num_particles)
  for (particle in 1:num_particles) {
    weights[particle] <- exp(-1/2*change_in_temp*t((simulated_data[particle, ] - likelihood_samples[particle, ])) %*% C_y_given_x_inv %*% (simulated_data[particle, ] - likelihood_samples[particle, ]))
    # weights[particle] <- exp(-1/2*change_in_temp*t((simulated_data[particle, ] - likelihood_samples[particle, ])) %*% (simulated_data[particle, ] - likelihood_samples[particle, ]))
  }
  # print(length(weights))
  # print(sum(weights))
  if (sum(weights) == 0) {
    print(change_in_temp)
    # print(simulated_data)
    # print(likelihood_samples)
  }
  weights <- weights/sum(weights)
  return(weights)
}

estimate_ess <- function(next_temp, current_temp, simulated_data, likelihood_samples, covariances, num_particles) {
  
  weights <- get_weights(next_temp, current_temp, simulated_data, likelihood_samples, covariances, num_particles)
  return(1/sum(weights**2))
}

get_ess_diff <- function(next_temp, current_temp, simulated_data, likelihood_samples, covariances, num_particles, target_ess) {
  
  ess <- estimate_ess(next_temp, current_temp, simulated_data, likelihood_samples, covariances, num_particles)
  return(ess - target_ess)
}

find_next_temp <- function(current_temp, simulated_data, likelihood_samples, covariances, num_particles, target_ess) {
  
  # print(current_temp)
  # print(dim(simulated_data))
  # print(dim(likelihood_samples))
  # print(covariances)
  # print(num_particles)
  # print(target_ess)
  
  # Fix the arguments of the estimate ess function
  estimate_ess_partial <- partial(get_ess_diff, current_temp = current_temp,
                                  simulated_data = simulated_data, likelihood_samples = likelihood_samples, 
                                  covariances = covariances, num_particles = num_particles, 
                                  target_ess = target_ess)
  
  # print(estimate_ess_partial(1))
  
  # We want the next temperature to be chosen such that the ESS equals the target ESS
  if (estimate_ess_partial(1) >= 0) {
    return(1.0)
  }
  else {
    next_temp <- uniroot(estimate_ess_partial, interval = c(current_temp, 1), tol = 1e-5)
    return(next_temp$root)
  }
  
  # print(current_temp)
  # next_temp <- uniroot(estimate_ess_partial, interval = c(current_temp, current_temp + 0.1), tol = 1e-5)
  # print(next_temp$root)
  # return(min(next_temp$root, 1))
  
  return(next_temp$root, 1)
  
}