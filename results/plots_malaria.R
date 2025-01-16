plot_eki_malaria <- function(eki_result) {
  
  d_in_particles <- exp(eki_result$particles[, 1]) + 0.16
  d_in_seq <- seq(min(d_in_particles), max(d_in_particles), length.out = 50)
  phi_particles <- plogis(eki_result$particles[, 2])
  eta0_particles <- plogis(eki_result$particles[, 3]) * 0.96 + 0.04
  sigma_particles <- exp(eki_result$particles[, 4])
  
  hist(d_in_particles, breaks = d_in_seq, freq = F, xlab = expression(d[in]))
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
  
  hist(d_in_particles, breaks = d_in_seq, freq = F, xlab = expression(d[in]))
  abline(v = 0.5, col = 'red')
  hist(phi_particles, freq = F, xlab = expression(phi))
  abline(v = 0.25, col = 'red')
  hist(eta0_particles, freq = F, xlab = expression(eta[0]))
  abline(v = 0.11, col = 'red')
  
}