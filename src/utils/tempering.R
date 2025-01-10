pacman::p_load(pacman, purrr)

get_weights <- function(next_temp, current_temp, ll_densities) {
  
  change_in_temp <- next_temp - current_temp
  weights <- exp(change_in_temp*ll_densities)
  
  if (sum(weights) == 0) {
    print(change_in_temp)
  }
  weights <- weights/sum(weights)
  return(weights)
}

estimate_ess <- function(next_temp, current_temp, ll_densities) {
  
  weights <- get_weights(next_temp, current_temp, ll_densities)
  return(1/sum(weights**2))
}

get_ess_diff <- function(next_temp, current_temp, ll_densities) {
  
  ess <- estimate_ess(next_temp, current_temp, ll_densities)
  return(ess - target_ess)
}

find_next_temp <- function(current_temp, ll_densities) {
  
  # Fix the arguments of the estimate ess function
  estimate_ess_partial <- partial(get_ess_diff, current_temp = current_temp, 
                                  ll_densities = ll_densities)
  
  # We want the next temperature to be chosen such that the ESS equals the target ESS
  if (estimate_ess_partial(1) >= 0) {
    return(1.0)
  }
  else {
    next_temp <- uniroot(estimate_ess_partial, interval = c(current_temp, 1), tol = 1e-5)
    return(next_temp$root)
  }
  
  return(next_temp$root, 1)
  
}