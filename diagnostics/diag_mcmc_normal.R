source('src/mcmc_normal.R')
prior_params = list(alpha.sd = 2, sigma2.sd = 5)

######################## 1 Dimensional Normal #################################

true_params = list(alpha = 2, sigma = 5, x = rep(1, 100))


run.1 <- normal_mcmc(true_params, prior_params)
chain.1 <- run.1$batch
run.2 <- normal_mcmc(true_params, prior_params)
chain.2 <- run.2$batch

## Alpha

# No burn-in
plot(chain.1[, 1], type = 'l')
lines(chain.2[, 1], col = 'red')

# A burn-in of 500
plot(chain.1[500:1000, 1], type = 'l')
lines(chain.2[500:1000, 1], col = 'red')

## Sigma2

# No burn-in
plot(exp(chain.1[, 2]), type = 'l')
lines(exp(chain.2[, 2]), col = 'red')

# A burn-in of 500
plot(exp(chain.1[500:1000, 2]), type = 'l')
lines(exp(chain.2[500:1000, 2]), col = 'red')
