pacman::p_load(pacman, MASS)

source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/models/eki_g_and_k_k_only.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/utils/eki_helper.R')

num_particles <- 400
initial_particles <- initialise_g_and_k_particles_k_only(num_particles, c())

true_params <- list(a = qnorm(3/10), b = qnorm(1/10), g = qnorm(2/10), k = qnorm(1/20))
# We make a single draw from the likelihood using the true (unknown parameters)
true_data <- likelihood_g_and_k_summary(true_params)
d_y <- length(true_data)
# I'm replicating this data for the number of particles to make the dimensions easier to work with
simulated_data <- matrix(true_data, nrow = num_particles, ncol = d_y, byrow = T)

likelihood_samples <- synthetic_g_and_k_k_only_summary(num_particles, initial_particles, true_params)

covariances <- calculate_covariances(initial_particles, likelihood_samples)

covariances$C_xx
covariances$C_yy
covariances$C_yx
covariances$C_y_given_x

ginv((covariances$C_yy + 9*covariances$C_y_given_x)) %*% covariances$C_yx

(simulated_data - likelihood_samples) %*% ginv((covariances$C_yy + 9*covariances$C_y_given_x)) %*% covariances$C_yx
(simulated_data - likelihood_samples) %*% ginv((covariances$C_yy)) %*% covariances$C_yx

gk_seq <- seq(-1000, 1000, length.out = 100)

plot(gk_seq, dgk(gk_seq, 3, 1, 2, 1/2), type = 'l')
lines(gk_seq, dgk(gk_seq, 3, 1, 2, 8), col = 'blue')

hist(rgk(1000, 3, 1, 2, 2), freq = F)
lines(gk_seq, dgk(gk_seq, 3, 1, 2, 2), col = 'red')

cor(likelihood_samples, initial_particles)
covariances$C_yx