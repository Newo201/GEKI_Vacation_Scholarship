pacman::p_load(pacman, testthat)

source('samples_normal.R')

# Testing dimensions for prior samples
test_that('Testing dimensions of prior sampling functions', 
         {
           expect_length(alpha_prior_sample(2, 100), 100)
           expect_length(logsigma2_prior_sample(5, 100), 100)
         })

# Testing means of prior samples
test_mean <- function(sampler, std, iterations) {
  samples <- sampler(std, iterations)
  return(abs(mean(samples)))
}

test_that('Testing that means of prior sampling functions are close to true mean', 
          {
            expect_lt(test_mean(alpha_prior_sample, 2, 100), 0.1)
            expect_lt(test_mean(logsigma2_prior_sample, 5, 100), 0.1)
          })

# Testing standard deviations of prior samples
test_sd <- function(sampler, std, iterations) {
  samples <- sampler(std, iterations)
  return(abs(sd(samples) - std))
}

test_that('Testing that sd of prior sampling functions are close to true sd', 
          {
            expect_lt(test_sd(alpha_prior_sample, 2, 10000), 0.1)
            expect_lt(test_sd(logsigma2_prior_sample, 5, 10000), 0.1)
          })

# Test sampling consistency
test_mean_difference <- function(sampler, std, iterations) {
  samples1 <- sampler(std, iterations)
  samples2 <- sampler(std, iterations)
  return(abs(mean(samples1) - mean(samples2)))
}

test_that('Testing sampling consistency',
          {
            expect_lt(test_mean_difference(alpha_prior_sample, 2, 100), 0.1)
            expect_lt(test_mean_difference(logsigma2_prior_sample, 5, 100), 0.1)
          })

# Test that random seed is working correctly