########################### Helper Functions ##################################

plot_alpha_overlay <- function(alpha_eki_particles, alpha_mcmc_particles, true_params, prior_params) {
  
  true_alpha <- true_params$alpha
  alpha_sequence <- seq(min(c(alpha_eki_particles, alpha_mcmc_particles, true_alpha)), 
                        max(c(alpha_eki_particles, alpha_mcmc_particles, true_alpha)), 
                        length.out = 50)
  
  alpha_prior_density <- dnorm(alpha_sequence, mean = prior_params$alpha.mean, 
                               sd = prior_params$alpha.sd)
  
  alpha_eki_density <- density(alpha_eki_particles) 
  alpha_mcmc_density <- density(alpha_mcmc_particles) 
  
  max_density <- max(c(alpha_eki_density$y, alpha_mcmc_density$y))
  
  plot(alpha_eki_density, main = 'EKI and MCMC', xlab = expression(alpha), ylab = 'Density',
       cex.main = 2, cex.lab = 1.5, xlim = c(min(alpha_sequence), max(alpha_sequence)), 
       ylim = c(0, max_density))
  lines(alpha_mcmc_density, col = 'purple')
  lines(alpha_sequence, alpha_prior_density, col = 'blue')
  abline(v = true_alpha, col = 'red')
  
  legend("topleft", 
         legend = c("EKI Posterior", "MCMC Posterior", "Prior", "True Value"),
         col = c("black", "purple", "blue", "red"),
         lty = 1, cex = 0.8)
  
}

plot_sigma2_overlay <- function(sigma2_eki_particles, sigma2_mcmc_particles, true_params, prior_params) {
  
  true_sigma2 <- log(true_params$sigma**2)
  sigma2_sequence <- seq(min(c(sigma2_eki_particles, sigma2_mcmc_particles, true_sigma2)), 
                        max(c(sigma2_eki_particles, sigma2_mcmc_particles, true_sigma2)), 
                        length.out = 50)
  
  sigma2_prior_density <- dnorm(sigma2_sequence, mean = prior_params$sigma2.mean, 
                               sd = prior_params$sigma2.sd)
  
  sigma2_eki_density <- density(sigma2_eki_particles) 
  sigma2_mcmc_density <- density(sigma2_mcmc_particles) 
  
  max_density <- max(c(sigma2_eki_density$y, sigma2_mcmc_density$y))
  
  plot(sigma2_eki_density, main = 'EKI and MCMC', xlab = expression(log(sigma^2)), ylab = 'Density',
       cex.main = 2, cex.lab = 1.5, xlim = c(min(sigma2_sequence), max(sigma2_sequence)),
       ylim = c(0, max_density))
  lines(sigma2_mcmc_density, col = 'purple')
  lines(sigma2_sequence, sigma2_prior_density, col = 'blue')
  abline(v = true_sigma2, col = 'red')
  
  legend("topleft", 
         legend = c("EKI Posterior", "MCMC Posterior", "Prior", "True Value"),
         col = c("black", "purple", "blue", "red"),
         lty = 1, cex = 0.8)
  
}

plot_alpha_particles <- function(alpha_particles, true_params, prior_params, algorithm,
                                 kde = T) {
  
  true_alpha <- true_params$alpha
  alpha_sequence <- seq(min(c(alpha_particles, true_alpha)), max(c(alpha_particles, true_alpha)), 
                        length.out = 50)
  alpha_prior_density <- dnorm(alpha_sequence, mean = prior_params$alpha.mean, 
                               sd = prior_params$alpha.sd)

  
  if (kde) {
    alpha_post_density <- density(alpha_particles)
    plot(alpha_post_density, main = algorithm, xlab = expression(alpha), ylab = 'Density',
         cex.main = 2, cex.lab = 1.5, xlim = c(min(alpha_sequence), max(alpha_sequence)))
  } else {
    hist(alpha_particles, freq = F, main = algorithm, xlab = expression(alpha), ylab = 'Density',
         cex.main = 2, cex.lab = 1.5)
  }
  
  lines(alpha_sequence, alpha_prior_density, col = 'blue')
  abline(v = true_params$alpha, col = 'red')
  
  legend("topleft", 
         legend = c("Estimated Posterior", "Prior", "True Value"),
         col = c("black", "blue", "red"),
         lty = 1, cex = 0.8)
  
}

plot_sigma2_particles <- function(sigma2_particles, true_params, prior_params, algorithm, kde = T) {
  
  true_sigma2 <- log(true_params$sigma**2)
  sigma2_sequence <- seq(min(c(sigma2_particles, true_sigma2)), max(c(sigma2_particles, true_sigma2)),
                         length.out = 50)
  sigma2_prior_density <- dnorm(sigma2_sequence, mean = prior_params$sigma2.mean, 
                                sd = prior_params$sigma2.sd)

  
  if (kde) {
    sigma2_post_density <- density(sigma2_particles)
    plot(sigma2_post_density, main = algorithm, xlab = expression(log(sigma^2)), ylab = 'Density',
         cex.main = 2, cex.lab = 1.5, xlim = c(min(sigma2_sequence), max(sigma2_sequence)))
  } else {
    hist(sigma2_particles, freq = F, main = algorithm, 
         xlab = expression(log(sigma^2)), ylab = 'Density',
         cex.main = 2, cex.lab = 1.5)
  }

  lines(sigma2_sequence, sigma2_prior_density, col = 'blue')
  abline(v = true_sigma2, col = 'red')
  legend("topleft", 
         legend = c("Estimated Posterior", "Prior", "True Value"),
         col = c("black", "blue", "red"),
         lty = 1, cex = 0.8)
  
}

############################## MCMC ###########################################

plot_mcmc_trace_plots <- function(chain.1, chain.2, burnin = 500, iterations = 1e4) {
  
  # par(mfrow = c(1,2))
  
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
  sigma2_particles <- combined_chains[, 2]
  
  # par(mfrow = c(1,2))
  plot_alpha_particles(alpha_particles, true_params, prior_params, 'MCMC', kde = kde)
  plot_sigma2_particles(sigma2_particles, true_params, prior_params, 'MCMC', kde = kde)
  
}


############################### EKI #########################################

plot_eki_normal <- function(eki_result, true_params, prior_params, kde = T) {
  
  alpha_particles <- eki_result$particles[, 1]
  sigma2_particles <- eki_result$particles[, 2]
  
  # par(mfrow = c(1, 2))
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
         lty = 1, cex = 0.8)
}

plot_eki_normal_known_mean <- function(eki_result, true_params, prior_params, kde = T) {
  
  sigma2_particles <- eki_result$particles[, 1]
  
  par(mfrow = c(1, 1))
  
  plot_sigma2_particles(sigma2_particles, true_params, prior_params, 'EKI Known Mean', kde = kde)
  
}

########################## Combined ###############################

plot_overlay <- function(eki_result, chain.1, chain.2, true_params, 
                         prior_params, burnin = 500, iterations = 1e4) {
  
  alpha_eki_particles <- eki_result$particles[, 1]
  sigma2_eki_particles <- eki_result$particles[, 2]
  
  chain.1.burnin <- chain.1[(burnin + 1): iterations, ]
  chain.2.burnin <- chain.2[(burnin + 1): iterations, ]
  combined_chains <- rbind(chain.1.burnin, chain.2.burnin)
  
  alpha_mcmc_particles <- combined_chains[, 1]
  sigma2_mcmc_particles <- combined_chains[, 2]
  
  plot_alpha_overlay(alpha_eki_particles, alpha_mcmc_particles, true_params, prior_params)
  plot_sigma2_overlay(sigma2_eki_particles, sigma2_mcmc_particles, true_params, prior_params)
  
  
}