source('src/samples_normal.R')

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

parameters <- list(alpha = 2, sigma = 5, x = 5, alpha.sd = 5, sigma2.sd = 2)
num_observations <- 100
simulated_data <- likelihood_sample(parameters, num_observations)

hist(simulated_data, freq = F)
data_sequence <- seq(min(simulated_data), max(simulated_data), length = 40)
density <- dnorm(data_sequence, mean = parameters$alpha*parameters$x, sd = parameters$sigma)
lines(data_sequence, density, col = 'red')




