############################## MCMC ###########################################

plot_mcmc_trace_plots <- function(chain.1, chain.2, burnin = 0, iterations = 1e4) {
  
  par(mfrow = c(1,2))
  
  plot(chain.1[(burnin + 1):iterations, 1], type = 'l', main = 'MCMC', 
       xlab = 'Iteration', ylab = expression(alpha))
  lines(chain.2[(burnin + 1):iterations, 1], col = 'red')
  
  plot(chain.1[(burnin + 1):iterations, 2], type = 'l', main = 'MCMC', 
       xlab = 'Iteration', ylab = expression(sigma^2))
  lines(chain.2[(burnin + 1):iterations, 2], col = 'red')
}

plot_mcmc_histogram <- function(chain.1, chain.2, true_params, burnin = 0, 
                           iterations = 1e4) {
  chain.1.burnin <- chain.1[(burnin + 1): iterations, ]
  chain.2.burnin <- chain.2[(burnin + 1): iterations, ]
  combined_chains <- rbind(chain.1.burnin, chain.2.burnin)
  
  alpha_particles <- combined_chains[, 1]
  alpha_sequence <- seq(min(alpha_particles), max(alpha_particles), 
                        length.out = 20)
  alpha_prior_density <- dnorm(alpha_sequence, mean = prior_params$alpha.mean, 
                               sd = prior_params$alpha.sd)
  true_alpha <- true_params$alpha
  sigma2_particles <- combined_chains[, 2]
  sigma2_sequence <- seq(min(sigma2_particles), max(sigma2_particles),
                         length.out = 20)
  sigma2_prior_density <- dnorm(sigma2_sequence, mean = prior_params$sigma2.mean, 
                                sd = prior_params$sigma2.sd)
  
  par(mfrow = c(1,2))
  hist(alpha_particles, freq = F, main = 'MCMC', xlab = expression(alpha), ylab = 'Density')
  lines(alpha_sequence, alpha_prior_density, col = 'blue')
  abline(v = true_params$alpha, col = 'red')
  hist(combined_chains[, 2], freq = F, main = 'MCMC', xlab = expression(log(sigma^2)), ylab = '')
  lines(sigma2_sequence, sigma2_prior_density, col = 'blue')
  abline(v = log(true_params$sigma**2), col = 'red')
  
}


############################### EKI #########################################

plot_eki_normal <- function(eki_result, true_params, prior_params) {
  
  alpha_particles <- eki_result$particles[, 1]
  alpha_sequence <- seq(min(alpha_particles), max(alpha_particles), 
                        length.out = 20)
  alpha_prior_density <- dnorm(alpha_sequence, mean = prior_params$alpha.mean, 
                               sd = prior_params$alpha.sd)
  true_alpha <- true_params$alpha
  sigma2_particles <- eki_result$particles[, 2]
  sigma2_sequence <- seq(min(sigma2_particles), max(sigma2_particles),
                         length.out = 20)
  sigma2_prior_density <- dnorm(sigma2_sequence, mean = prior_params$sigma2.mean, 
                                sd = prior_params$sigma2.sd)
  true_sigma2 <- log(true_params$sigma**2)
  
  par(mfrow = c(1, 2))
  
  hist(alpha_particles, freq = F, main = 'EKI', 
       xlab = expression(alpha), ylab = 'Density')
  lines(alpha_sequence, alpha_prior_density, col = 'blue')
  abline(v = true_alpha, col = 'red')
  hist(sigma2_particles, freq = F, main = 'EKI', 
       xlab = expression(log(sigma^2)), ylab = '')
  lines(sigma2_sequence, sigma2_prior_density, col = 'blue')
  abline(v = true_sigma2, col = 'red')
}

# ToDo: add plot of posterior density
plot_eki_normal_known_var <- function(eki_result, true_data, true_params, prior_params) {
  
  alpha_particles <- eki_result$particles[, 1]
  alpha_sequence <- seq(min(alpha_particles), max(alpha_particles), 
                        length.out = 20)
  alpha_prior_density <- dnorm(alpha_sequence, mean = prior_params$alpha.mean, 
                               sd = prior_params$alpha.sd)
  
  true_alpha <- true_params$alpha
  
  par(mfrow = c(1, 1))
  
  hist(alpha_particles, freq = F, main = 'EKI Known Variance', 
       xlab = expression(alpha), ylab = 'Density')
  lines(alpha_sequence, alpha_prior_density, col = 'blue')
  abline(v = true_alpha, col = 'red')
  
  Q <- (prior_params$alpha.sd**2) * diag(1)
  m <- prior_params$alpha.mean * diag(1)
  H <- t(t(true_params$x))
  R <- (true_params$sigma**2)*diag(length(H))
  
  post_mean <- m + Q %*% t(H) %*% ginv(H %*% Q %*% t(H) + R) %*% (t(true_data) - H %*% m)
  post_var <- Q - Q %*% t(H) %*% ginv(H %*% Q %*% t(H) + R) %*% H %*% Q
  
  alpha_post_density <- dnorm(alpha_sequence, mean = post_mean, sd = sqrt(post_var))
  lines(alpha_sequence, alpha_post_density, col = 'green')
}

plot_eki_normal_known_mean <- function(eki_result, true_params, prior_params) {
  
  sigma2_particles <- eki_result$particles[, 1]
  sigma2_sequence <- seq(min(sigma2_particles), max(sigma2_particles),
                         length.out = 20)
  sigma2_prior_density <- dnorm(sigma2_sequence, mean = prior_params$sigma2.mean, 
                                sd = prior_params$sigma2.sd)
  true_sigma2 <- log(true_params$sigma**2)
  
  par(mfrow = c(1, 1))
  
  hist(sigma2_particles, freq = F, main = 'EKI Known Mean', 
       xlab = expression(log(sigma^2)), ylab = 'Density')
  lines(sigma2_sequence, sigma2_prior_density, col = 'blue')
  abline(v = true_sigma2, col = 'red')
}