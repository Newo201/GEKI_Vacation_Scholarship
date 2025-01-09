source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/eki_normal.R')

num_particles <- 400
num_dimensions <- 100
true_params_nd = list(alpha = 3, sigma = 2, x = rep(1, num_dimensions), alpha.sd = 5, sigma2.sd = 2)
eki_result <- eki_normal_adaptive(num_particles, true_params_nd)
# eki_result <- eki_normal(num_particles, true_params_nd)

par(mfrow = c(1, 1))
initial_particles <- initialise_particles(num_particles, true_params_nd)
sigma2_prior <- exp(initial_particles[, 2])
particle_sequence <- seq(min(sigma2_prior), max(sigma2_prior), length.out = 50)

density <- dlnorm(particle_sequence, meanlog = 0, sdlog = 2)
hist(sigma2_prior, freq = F, breaks = particle_sequence)
lines(particle_sequence, density, col = 'red')

plot_eki_normal <- function(eki_result) {
  
  alpha_particles <- eki_result$particles[, 1]
  sigma2_particles <- exp(eki_result$particles[, 2])
  
  par(mfrow = c(1, 2))
  
  hist(alpha_particles, freq = F, main = 'Alpha', xlab = '')
  hist(sigma2_particles, freq = F, main = 'Sigma2', xlab = '', ylab = '')
}

plot_eki_normal(eki_result)

initial_particles <- initialise_particles(num_particles, true_params_nd)
cov(initial_particles)

plot_temp_sequence <- function(eki_result) {
  temp_sequence = eki_result$temp_sequence
  plot(1:length(temp_sequence), temp_sequence)
}

plot_temp_sequence(eki_result)