pacman::p_load(pacman, testthat, deSolve, rootSolve)

source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/models/eki_malaria.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/samples/samples_malaria.R')

data_path = "C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/data/Malariah_data.rds"
true_data = log(readRDS(data_path))

######################## Sensitivity Analays ##################################

# Observations decrease as d_in increases
d_in_seq <- seq(0.2, 10, by = 0.1)
res <- c()
for (d_in in d_in_seq) {
  constrained_params = list(d_in = d_in, phi = 0.75, eta0 = 0.04)
  res <- c(res, mean(log(diff(likelihood_malaria_mean(constrained_params)))))
}
plot(d_in_seq, res, ylim = c(8, 12))
abline(h = mean(true_data))

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

######################### Plotting Data #################################

constrained_params = list(d_in = 5, phi = 0.5, eta0 = 0.04, sigma = 0.5)
unconstrained_params <- unconstrain_malaria_params(constrained_params)
time = seq(1/12, 10.75, by = 1/12)
means <- log(diff(likelihood_malaria_mean(constrained_params)))
means
# samples <- likelihood_malaria(unconstrained_params)
# 
# print(mean(res))
# plot(time, means, type = 'l', col = 'red', ylim = c(min(min(true_data), min(res)), max(max(true_data, max(res)))))
# lines(time, samples, col = 'blue')
# lines(time, true_data)
# 
# res
# 
# f <- function(x, y, z) {
#   x + y + z
# }
# vf <- Vectorize(f, vectorize.args = c("x", "y"), SIMPLIFY = F)
# f(1:3, 1:3, 2)
# vf(1:3, 1:3)
# vf(y = 1:3)
