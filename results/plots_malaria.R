pacman::p_load(pacman, matrixStats, glue)

plot_eki_posterior_predictive <- function(eki_result, true_data, true_params) {
  
  final_particles <- eki_result$particles
  
  num_particles <- dim(final_particles)[1]
  
  # Generate samples using the latest particles
  likelihood_prediction <- synthetic_malaria_known_var(num_particles, final_particles, true_params)
  
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
  
  plot(time_seq, pred.mean, type = 'l', ylim = c(min(pred.lower), max(pred.upper)))
  polygon(c(time_seq, rev(time_seq)), c(pred.lower, rev(pred.upper)), col = 'lightblue', border = F)
  points(time_seq, true_data)
  lines(time_seq, pred.mean, lwd = 2, col = 'darkblue')
  # lines(time_seq, pred.lower, col = 'blue')
  # lines(time_seq, pred.upper, col = 'blue')
  
}

plot_eki_malaria <- function(eki_result) {
  
  d_in_particles <- exp(eki_result$particles[, 1]) + 0.16
  d_in_seq <- seq(min(d_in_particles), max(d_in_particles), length.out = 50)
  phi_particles <- plogis(eki_result$particles[, 2])
  eta0_particles <- plogis(eki_result$particles[, 3]) * 0.96 + 0.04
  sigma_particles <- exp(eki_result$particles[, 4])
  
  hist(d_in_particles, breaks = d_in_seq, freq = F, xlab = expression(d[inf]))
  abline(v = 0.5, col = 'red')
  hist(phi_particles, freq = F, xlab = expression(phi))
  abline(v = 0.25, col = 'red')
  hist(eta0_particles, freq = F, xlab = expression(eta[0]))
  abline(v = 0.11, col = 'red')
  hist(sigma_particles, freq = F, xlab = expression(sigma))
  abline(v = 0.5, col = 'red')
  
}

plot_eki_malaria_known_var <- function(eki_result) {
  
  d_in_particles <- exp(eki_result$particles[, 1]) + 0.16
  d_in_seq <- seq(min(d_in_particles), max(d_in_particles), length.out = 50)
  phi_particles <- plogis(eki_result$particles[, 2])
  eta0_particles <- plogis(eki_result$particles[, 3]) * 0.96 + 0.04
  
  hist(d_in_particles, breaks = d_in_seq, freq = F, xlab = expression(d[inf]))
  abline(v = 0.5, col = 'red')
  hist(phi_particles, freq = F, xlab = expression(phi))
  abline(v = 0.25, col = 'red')
  hist(eta0_particles, freq = F, xlab = expression(eta[0]))
  abline(v = 0.11, col = 'red')
  
}

plot_eki_malaria_d_in_only <- function(eki_result) {
  
  d_in_particles <- exp(eki_result$particles[, 1]) + 0.16
  d_in_seq <- seq(min(d_in_particles), max(d_in_particles), length.out = 50)
  
  hist(d_in_particles, breaks = d_in_seq, freq = F, xlab = expression(d[inf]))
  abline(v = 0.5, col = 'red')
  
}

m <- matrix(rnorm(10), nrow = 10, ncol = 10)
colQuantiles(m, probs = 0.25)

next_temp = 2
print(glue("Next temp is {next_temp}"))