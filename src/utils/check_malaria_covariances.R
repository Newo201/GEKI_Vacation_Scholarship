####################### Covariance Functions ##################################

check_c_yy <- function(num_particles, prior_params, true_params, column) {
  
  particles <- initialise_malaria_particles(num_particles, prior_params)
  likelihood_samples <- synthetic_malaria(num_particles, particles, true_params)
  covariances <- calculate_covariances(particles, likelihood_samples)
  
  return(mean(diag(covariances$C_yy)))
}

check_c_y_given_x <- function(num_particles, prior_params, true_params, column) {
  particles <- initialise_malaria_particles(num_particles, prior_params)
  likelihood_samples <- synthetic_malaria(num_particles, particles, true_params)
  
  covariances <- calculate_covariances(particles, likelihood_samples)
  
  return(mean(diag(covariances$C_y_given_x)))
}

check_c_yx <- function(num_particles, prior_params, true_params, column) {
  particles <- initialise_malaria_particles(num_particles, prior_params)
  likelihood_samples <- synthetic_malaria(num_particles, particles, true_params)
  
  covariances <- calculate_covariances(particles, likelihood_samples)
  
  mean(covariances$C_yx[, column])
}

check_stepsize_var <- function(num_particles, prior_params, true_params, column) {
  
  particles <- initialise_malaria_particles(num_particles, prior_params)
  likelihood_samples <- synthetic_malaria(num_particles, particles, true_params)
  covariances <- calculate_covariances(particles, likelihood_samples)
  
  # print(mean(diag(covariances$C_y_given_x)))
  
  # # print(covariances$C_yx)
  stepsize_var <- ginv(covariances$C_yy + covariances$C_y_given_x)
  return(mean(diag(stepsize_var)))
}

check_stepsize <- function(num_particles, prior_params, true_params, column) {
  particles <- initialise_malaria_particles(num_particles, prior_params)
  likelihood_samples <- synthetic_malaria(num_particles, particles, true_params)
  
  covariances <- calculate_covariances(particles, likelihood_samples)
  
  # print(mean(diag(covariances$C_y_given_x)))
  
  # # print(covariances$C_yx)
  stepsize_var <- ginv(covariances$C_yy)
  stepsize <- stepsize_var %*% covariances$C_yx
  return(mean(stepsize[, column]))
}

############################# Plots ###################################

plot_covariance_against_sigma <- function(covariance_func, covariance_type, column) {
  sigma_seq <- seq(-2, 0, length.out = 5)
  res <- c()
  for (sigma in sigma_seq) {
    prior_params <- list(din.sd = 2,
                         phi.mean = 0, phi.sd = 1,
                         eta0.mean = 0, eta0.sd = 1,
                         sigma.mean = sigma, sigma.sd = 1)
    true_params <- list(sigma = 0.5, phi = 0.25, eta0 = 0.11, d_in = 0.5)
    res <- c(res, covariance_func(50, prior_params, true_params, column))
  }
  print(res)
  plot(exp(sigma_seq), res, xlab = expression(sigma), ylab = 'Covariance', 
       main = glue('Empirical {covariance_type}'), cex.main = 2, cex.lab = 1.5)
  abline(h = 0, col = 'red')
}

# plot_covariance_against_sigma2_dispersion <- function(covariance_func) {
#   sigma2_seq <- seq(0.01, 3, length.out = 50)
#   res <- c()
#   for (sigma2 in sigma2_seq) {
#     prior_params = list(alpha.mean = 0, alpha.sd = 5, sigma2.mean = 0, sigma2.sd = sigma2)
#     true_params = list(alpha = 2, sigma = 2, x = rep(0, 50))
#     res <- c(res, covariance_func(1000, prior_params, true_params))
#   }
#   plot(sigma2_seq, res, xlab = 'sigma2_prior_variance', ylab ='')
# }
# 
# plot_covariance_against_alpha <- function(covariance_func) {
#   alpha_seq <- seq(-2, 10, length.out = 50)
#   res <- c()
#   for (alpha in alpha_seq) {
#     prior_params = list(alpha.mean = alpha, alpha.sd = 5, sigma2.mean = 0, sigma2.sd = 2)
#     true_params <- list(alpha = alpha, sigma = 2, x = rep(0, 50))
#     res <- c(res, covariance_func(1000, prior_params, true_params))
#   }
#   plot(alpha_seq, res, xlab = 'alpha', ylab = '')
# s
