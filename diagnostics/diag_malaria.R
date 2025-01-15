pacman::p_load(pacman, testthat, deSolve, rootSolve)

source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/models/eki_malaria.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/samples/samples_malaria.R')

# Observations decrease as d_in increases
d_in_seq <- seq(0.2, 50, by = 0.5)
res <- c()
for (d_in in d_in_seq) {
  constrained_params = list(d_in = d_in, phi = 0.25, eta0 = 0.11)
  res <- c(res, mean(log(diff(likelihood_malaria_mean(constrained_params)))))
}
plot(d_in_seq, res)

# We get a wave like pattern as the value of phi changes
phi_seq <- seq(0, 1, by = 0.05)
res <- c()
for (phi in phi_seq) {
  constrained_params = list(d_in = 0.5, phi = phi, eta0 = 0.11)
  res <- c(res, mean(log(diff(likelihood_malaria_mean(constrained_params)))))
}
plot(phi_seq, res)

# Observations increase as eta0 increases
eta0_seq <- seq(0.05, 1, by = 0.05)
res <- c()
for (eta0 in eta0_seq) {
  constrained_params = list(d_in = 0.5, phi = 0.25, eta0 = eta0)
  res <- c(res, mean(log(diff(likelihood_malaria_mean(constrained_params)))))
}
plot(eta0_seq, res)

