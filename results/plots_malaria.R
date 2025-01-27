pacman::p_load(pacman, matrixStats, glue)

############################# Helper Functions #####################################

plot_d_in_particles <- function(d_in_particles, true_params, prior_params, kde = T) {
  
  d_in_seq <- seq(min(d_in_particles), max(d_in_particles), length.out = 50)
  # ToDo: add parameters for prior
  d_in_prior <- dhnorm(d_in_seq - 0.16, sigma = prior_params$din.sd)
  
  if (kde) {
    d_in_post_density <- density(d_in_particles)
    plot(d_in_post_density)
  } else {
    hist(d_in_particles, breaks = d_in_seq, freq = F, xlab = expression(d[inf]))
  }

  lines(d_in_seq, d_in_prior, col = 'blue')
  abline(v = true_params$d_in, col = 'red')
  
  legend("topleft", 
         legend = c("Estimated Posterior", "Prior", "True Value"),
         col = c("black", "blue", "red"),
         pch = 16)
  
}

plot_phi_particles <- function(phi_particles, true_params, prior_params, kde = T) {
  
  phi_seq <- seq(min(phi_particles), max(phi_particles), length.out = 50)
  phi_prior <- dlogitnorm(phi_seq, mu = prior_params$phi.mean, sigma = prior_params$phi.sd)
  
  if (kde) {
    phi_post_density <- density(phi_particles)
    plot(phi_post_density)
  } else {
    hist(phi_particles, breaks = phi_seq, freq = F, xlab = expression(phi))
  }
  
  lines(phi_seq, phi_prior, col = 'blue')
  abline(v = true_params$phi, col = 'red')
  
  legend("topleft", 
         legend = c("Estimated Posterior", "Prior", "True Value"),
         col = c("black", "blue", "red"),
         pch = 16)
}

plot_eta0_particles <- function(eta0_particles, true_params, prior_params, kde = T) {
  
  eta0_seq <- seq(min(eta0_particles), max(eta0_particles), length.out = 50)
  # Todo: fix up prior
  eta0_prior <- dlogitnorm(eta0_seq, mu = prior_params$eta0.mean, sigma = prior_params$eta0.sd)
  
  if (kde) {
    eta0_post_density <- density(eta0_particles)
    plot(eta0_post_density)
  } else {
    hist(eta0_particles, breaks = eta0_seq, freq = F, xlab = expression(eta[0]))
  }
  
  lines(eta0_seq, eta0_prior, col = 'blue')
  abline(v = true_params$eta0, col = 'red')
  
  legend("topleft", 
         legend = c("Estimated Posterior", "Prior", "True Value"),
         col = c("black", "blue", "red"),
         pch = 16)
}

plot_sigma_particles <- function(sigma_particles, true_params, prior_params, kde = kde) {
  
  sigma_seq <- seq(min(sigma_particles), max(sigma_particles), length.out = 50)
  # Add parameters
  sigma_prior <- dlnorm(sigma_seq, mean = prior_params$sigma.mean, sd = prior_params$sigma.sd)
  
  if (kde) {
    sigma_post_density <- density(sigma_particles)
    plot(sigma_post_density)
  } else {
    hist(sigma_particles, freq = F, xlab = expression(sigma))
  }
  
  lines(sigma_seq, sigma_prior, col = 'blue')
  abline(v = true_params$sigma, col = 'red')
  
  legend("topleft", 
         legend = c("Estimated Posterior", "Prior", "True Value"),
         col = c("black", "blue", "red"),
         pch = 16)
  
}

############################# Posterior Predictive ##################################

plot_eki_posterior_predictive_known_var <- function(eki_result, true_data, true_params) {
  
  final_particles <- eki_result$particles
  
  num_particles <- dim(final_particles)[1]
  
  # Generate samples using the latest particles
  likelihood_prediction <- synthetic_malaria_known_var(num_particles, final_particles, true_params)
  
  time_seq <- seq(1/12, 10.75, by = 1/12)
  
  # Find the 2.5% and 97.5% quantiles
  pred.lower <- colQuantiles(likelihood_prediction, probs = 0.025)
  pred.mean <- colMeans(likelihood_prediction)
  pred.upper <- colQuantiles(likelihood_prediction, probs = 0.975)
  
  plot(time_seq, pred.mean, type = 'l', ylim = c(min(pred.lower), max(pred.upper)), 
       xlab = 'Time', ylab = 'New cases (log scale)', main = '95% Posterior Predictive Distribution',
       cex.main = 1.5, cex.lab = 1.2)
  polygon(c(time_seq, rev(time_seq)), c(pred.lower, rev(pred.upper)), col = 'lightblue', border = F)
  points(time_seq, true_data)
  lines(time_seq, pred.mean, lwd = 2, col = 'darkblue')
  # lines(time_seq, pred.lower, col = 'blue')
  # lines(time_seq, pred.upper, col = 'blue')
  
}

plot_eki_posterior_predictive_known_d_in <- function(eki_result, true_data, true_params) {
  
  final_particles <- eki_result$particles
  
  num_particles <- dim(final_particles)[1]
  
  # Generate samples using the latest particles
  likelihood_prediction <- synthetic_malaria_known_d_in(num_particles, final_particles, true_params)
  
  time_seq <- seq(1/12, 10.75, by = 1/12)
  
  # Find the 2.5% and 97.5% quantiles
  pred.lower <- colQuantiles(likelihood_prediction, probs = 0.025)
  pred.mean <- colMeans(likelihood_prediction)
  pred.upper <- colQuantiles(likelihood_prediction, probs = 0.975)
  
  plot(time_seq, pred.mean, type = 'l', ylim = c(min(pred.lower), max(pred.upper)), 
       xlab = 'Time', ylab = 'New cases (log scale)', main = '95% Posterior Predictive Distribution',
       cex.main = 1.5, cex.lab = 1.2)
  polygon(c(time_seq, rev(time_seq)), c(pred.lower, rev(pred.upper)), col = 'lightblue', border = F)
  points(time_seq, true_data)
  lines(time_seq, pred.mean, lwd = 2, col = 'darkblue')
  # lines(time_seq, pred.lower, col = 'blue')
  # lines(time_seq, pred.upper, col = 'blue')
  
}

plot_eki_posterior_predictive_d_in_only <- function(eki_result, true_data, true_params) {
  
  final_particles <- eki_result$particles
  
  num_particles <- dim(final_particles)[1]
  
  # Generate samples using the latest particles
  likelihood_prediction <- synthetic_malaria_d_in_only(num_particles, final_particles, true_params)
  
  time_seq <- seq(1/12, 10.75, by = 1/12)
  
  # Find the 2.5% and 97.5% quantiles
  pred.lower <- colQuantiles(likelihood_prediction, probs = 0.025)
  pred.mean <- colMeans(likelihood_prediction)
  pred.upper <- colQuantiles(likelihood_prediction, probs = 0.975)
  
  plot(time_seq, pred.mean, type = 'l', ylim = c(min(pred.lower), max(pred.upper)), 
       xlab = 'Time', ylab = 'New cases (log scale)', main = '95% Posterior Predictive Distribution',
       cex.main = 1.5, cex.lab = 1.2)
  polygon(c(time_seq, rev(time_seq)), c(pred.lower, rev(pred.upper)), col = 'lightblue', border = F)
  points(time_seq, true_data)
  lines(time_seq, pred.mean, lwd = 2, col = 'darkblue')
  # lines(time_seq, pred.lower, col = 'blue')
  # lines(time_seq, pred.upper, col = 'blue')
  
}

################################### Prior Predictive ##############################

plot_eki_prior_predictive_d_in_only <- function(prior_particles, true_data, true_params) {
  
  num_particles <- dim(prior_particles)[1]
  
  # Generate samples using the latest particles
  likelihood_prediction <- synthetic_malaria_d_in_only(num_particles, prior_particles, true_params)
  
  time_seq <- seq(1/12, 10.75, by = 1/12)
  
  # Find the 2.5% and 97.5% quantiles
  pred.lower <- colQuantiles(likelihood_prediction, probs = 0.025)
  pred.mean <- colMeans(likelihood_prediction)
  pred.upper <- colQuantiles(likelihood_prediction, probs = 0.975)
  
  plot(time_seq, pred.mean, type = 'l', ylim = c(min(pred.lower), max(pred.upper)))
  polygon(c(time_seq, rev(time_seq)), c(pred.lower, rev(pred.upper)), col = 'lightblue', border = F)
  points(time_seq, true_data)
  lines(time_seq, pred.mean, lwd = 2, col = 'darkblue')
  # lines(time_seq, pred.lower, col = 'blue')
  # lines(time_seq, pred.upper, col = 'blue')
  
}

############################## Marginal Posteriors ################################

plot_eki_malaria <- function(eki_result, true_params, prior_params, kde = T) {
  
  d_in_particles <- exp(eki_result$particles[, 1]) + 0.16
  phi_particles <- plogis(eki_result$particles[, 2])
  eta0_particles <- plogis(eki_result$particles[, 3])
  sigma_particles <- exp(eki_result$particles[, 4])
  
  plot_d_in_particles(d_in_particles, true_params, prior_params, kde = kde)
  plot_phi_particles(phi_particles, true_params, prior_params, kde = kde)
  plot_eta0_particles(eta0_particles, true_params, prior_params, kde = kde)
  plot_sigma_particles(sigma_particles, true_params, prior_params, kde = kde)
  
}

plot_eki_malaria_known_var <- function(eki_result, true_params, prior_params, kde = T) {
  
  d_in_particles <- exp(eki_result$particles[, 1]) + 0.16
  phi_particles <- plogis(eki_result$particles[, 2])
  eta0_particles <- plogis(eki_result$particles[, 3])
  
  plot_d_in_particles(d_in_particles, true_params, prior_params, kde = kde)
  plot_phi_particles(phi_particles, true_params, prior_params, kde = kde)
  plot_eta0_particles(eta0_particles, true_params, prior_params, kde = kde)
  
}

plot_eki_malaria_known_d_in <- function(eki_result, true_params, prior_params, kde = T) {
  
  phi_particles <- plogis(eki_result$particles[, 1])
  eta0_particles <- plogis(eki_result$particles[, 2])
  
  plot_phi_particles(phi_particles, true_params, prior_params, kde = kde)
  plot_eta0_particles(eta0_particles, true_params, prior_params, kde = kde)
  
}


plot_eki_malaria_d_in_only <- function(eki_result, true_params, prior_params, kde = T) {
  
  d_in_particles <- exp(eki_result$particles[, 1]) + 0.16
  plot_d_in_particles(d_in_particles, true_params, prior_params, kde = kde)
  
}
