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




