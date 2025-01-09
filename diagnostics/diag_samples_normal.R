source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/samples_normal.R')
getwd()

########################### Prior Samples #####################################

## alpha

alpha.sd <- 2
alpha_samples <- alpha_prior_sample(alpha.sd, 1000)
alpha_sequence <- seq(min(alpha_samples), max(alpha_samples), length = 40)
density <- dnorm(alpha_sequence, mean = 0, sd = alpha.sd)
hist(alpha_samples, freq = F)
lines(alpha_sequence, density, col = 'red')

## sigma2

sigma2.sd <- 5
logsigma2_samples <- logsigma2_prior_sample(sigma2.sd, 1000)
logsigma2_sequence <- seq(min(logsigma2_samples), max(logsigma2_samples), length = 40)
density <- dnorm(logsigma2_sequence, mean = 0, sd = sigma2.sd)
hist(logsigma2_samples, freq = F)
lines(logsigma2_sequence, density, col = 'red')

######################### Likelihood Samples ################################

par(mfrow = c(1,1))
num_dimensions <- 100
parameters <- list(alpha = 2, sigma = 5, x = rep(1, num_dimensions), alpha.sd = 5, sigma2.sd = 2)
simulated_data <- likelihood_sample(parameters, 1)

hist(simulated_data, freq = F)
data_sequence <- seq(min(simulated_data), max(simulated_data), length = 40)
density <- dnorm(data_sequence, mean = parameters$alpha*parameters$x[1], sd = parameters$sigma)
lines(data_sequence, density, col = 'red')




