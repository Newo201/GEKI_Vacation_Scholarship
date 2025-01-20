########################### Helper Functions ##################################

plot_alpha_particles <- function(alpha_particles, true_params, prior_params, algorithm,
                                 kde = T) {
  
  alpha_sequence <- seq(min(alpha_particles), max(alpha_particles), 
                        length.out = 20)
  alpha_prior_density <- dnorm(alpha_sequence, mean = prior_params$alpha.mean, 
                               sd = prior_params$alpha.sd)
  true_alpha <- true_params$alpha
  
  if (kde) {
    alpha_post_density <- density(alpha_particles)
    plot(alpha_post_density, main = algorithm, xlab = expression(alpha), ylab = 'Density')
  } else {
    hist(alpha_particles, freq = F, main = algorithm, xlab = expression(alpha), ylab = 'Density')
  }
  
  lines(alpha_sequence, alpha_prior_density, col = 'blue')
  abline(v = true_params$alpha, col = 'red')
  
  legend("topleft", 
         legend = c("Estimated Posterior", "Prior", "True Value"),
         col = c("black", "blue", "red"),
         pch = 16)
  
}

plot_sigma2_particles <- function(sigma2_particles, true_params, prior_params, algorithm, kde = T) {
  
  sigma2_sequence <- seq(min(sigma2_particles), max(sigma2_particles),
                         length.out = 20)
  sigma2_prior_density <- dnorm(sigma2_sequence, mean = prior_params$sigma2.mean, 
                                sd = prior_params$sigma2.sd)
  true_sigma2 <- true_params$sigma**2
  
  if (kde) {
    sigma2_post_density <- density(sigma2_particles)
    plot(sigma2_post_density, main = algorithm, xlab = expression(sigma^2), ylab = 'Density')
  } else {
    hist(sigma2_particles, freq = F, main = algorithm, 
         xlab = expression(sigma^2), ylab = 'Density')
  }

  lines(sigma2_sequence, sigma2_prior_density, col = 'blue')
  abline(v = true_sigma2, col = 'red')
  legend("topleft", 
         legend = c("Estimated Posterior", "Prior", "True Value"),
         col = c("black", "blue", "red"),
         pch = 16)
  
}

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

plot_mcmc_histogram <- function(chain.1, chain.2, true_params, prior_params, 
                                burnin = 0, iterations = 1e4, kde = T) {
  
  chain.1.burnin <- chain.1[(burnin + 1): iterations, ]
  chain.2.burnin <- chain.2[(burnin + 1): iterations, ]
  combined_chains <- rbind(chain.1.burnin, chain.2.burnin)
  
  alpha_particles <- combined_chains[, 1]
  sigma2_particles <- exp(combined_chains[, 2])
  
  par(mfrow = c(1,2))
  plot_alpha_particles(alpha_particles, true_params, prior_params, 'MCMC', kde = kde)
  plot_sigma2_particles(sigma2_particles, true_params, prior_params, 'MCMC', kde = kde)
  
}


############################### EKI #########################################

plot_eki_normal <- function(eki_result, true_params, prior_params, kde = T) {
  
  alpha_particles <- eki_result$particles[, 1]
  sigma2_particles <- eki_result$particles[, 2]
  
  par(mfrow = c(1, 2))
  plot_alpha_particles(alpha_particles, true_params, prior_params, 'EKI Normal', kde = kde)
  
  plot_sigma2_particles(sigma2_particles, true_params, prior_params, 'EKI Normal', kde = kde)
  
}

# ToDo: add plot of posterior density
plot_eki_normal_known_var <- function(eki_result, true_data, true_params, 
                                      prior_params, kde = T) {
  
  alpha_particles <- eki_result$particles[, 1]
  alpha_sequence <- seq(min(alpha_particles), max(alpha_particles), 
                        length.out = 20)
  
  par(mfrow = c(1, 1))
  
  plot_alpha_particles(alpha_particles, true_params, prior_params, 
                       'EKI Known Var', kde = kde)
  Q <- (prior_params$alpha.sd**2) * diag(1)
  m <- prior_params$alpha.mean * diag(1)
  H <- t(t(true_params$x))
  R <- (true_params$sigma**2)*diag(length(H))
  
  post_mean <- m + Q %*% t(H) %*% ginv(H %*% Q %*% t(H) + R) %*% (t(true_data) - H %*% m)
  post_var <- Q - Q %*% t(H) %*% ginv(H %*% Q %*% t(H) + R) %*% H %*% Q
  
  alpha_post_density <- dnorm(alpha_sequence, mean = post_mean, sd = sqrt(post_var))
  lines(alpha_sequence, alpha_post_density, col = 'green')
  legend("topleft", 
         legend = c("Estimated Posterior", "Prior", "True Value", "Analytical Posterior"),
         col = c("black", "blue", "red", "green"),
         pch = 16)
}

plot_eki_normal_known_mean <- function(eki_result, true_params, prior_params, kde = T) {
  
  sigma2_particles <- eki_result$particles[, 1]
  
  par(mfrow = c(1, 1))
  
  plot_sigma2_particles(sigma2_particles, true_params, prior_params, 'EKI Known Mean', kde = kde)
  
}