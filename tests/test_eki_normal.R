source('src/eki_normal.R')
source('src/utils.R')

pacman::p_load(pacman, testthat, purrr, matrixcalc)

num_particles <- 50

######################### Fixtures #######################################

# ToDo: see if there is a proper way of implementing fixtures like in Python

true_params_1d = list(alpha = 2, sigma = 2, x = 1, alpha.sd = 5, sigma2.sd = 2)
true_params_2d = list(alpha = 2, sigma = 2, x = c(1, 2), alpha.sd = 5, sigma2.sd = 2)

likelihood_samples_1d <- likelihood_sample(true_params_1d, num_particles)
simulated_data_1d <- matrix(likelihood_sample(true_params_1d, 1), nrow = num_particles, ncol = 1)

likelihood_samples_2d <- likelihood_sample(true_params_2d, num_particles)
simulated_data_2d <- matrix(likelihood_sample(true_params_2d, 1), nrow = num_particles, ncol = 1)

particles <- initialise_particles(num_particles, true_params_1d)

######################### Testing For Dimensions ##############################

test_conditional_cov <- function(particles, likelihood_samples) {
  
  # Calculate the covariance matrices
  C_xx = cov(particles)
  # print(C_xx)
  C_yy = cov(likelihood_samples)
  C_xy = cov(particles, likelihood_samples)
  C_yx = cov(likelihood_samples, particles)
  
  C_y_given_x = C_yy - C_yx %*% solve(C_xx) %*% C_xy
  return(C_y_given_x)
}

generate_pertubations <- function(num_particles, particles, likelihood_samples) {
  d_y <- dim(likelihood_samples)[2]
  C_y_given_x <- test_conditional_cov(particles, likelihood_samples)
  # Generate perturbations
  eta <- matrix(data = rmvnorm(n = num_particles, mean = rep(0, d_y), sigma = C_y_given_x), nrow = num_particles, ncol = d_y)
  return(eta)
}

test_that('Dimensions of particles are correct', {
  expect_equal(dim(particles), c(50, 2))
  expect_equal(dim(particles), c(50, 2))
})

test_that('Dimensions of weights are correct', {
  expect_equal(length(get_weights(0.1, 0, likelihood_samples_1d, simulated_data_1d)), num_particles)
})

test_that('Dimensions of pertubations are correct', {
  expect_equal(dim(generate_pertubations(num_particles, particles, likelihood_samples_1d)), c(num_particles, 1))
  expect_equal(dim(generate_pertubations(num_particles, particles, likelihood_samples_2d)), c(num_particles, 2))
})

test_that('Dimensions of particles are correct after updating', {
  expect_equal(dim(update_particles(0.1, particles, simulated_data_1d, likelihood_samples_1d, num_particles)), dim(particles))
})

######################## Testing For Covariance Matrices #####################

test_that('Dimensions of covariance matrices are correct', {
  expect_equal(dim(cov(likelihood_samples_1d)), c(1, 1))
  expect_equal(dim(cov(likelihood_samples_2d)), c(2, 2))
  expect_equal(dim(cov(particles)), c(2,2))
  expect_equal(dim(cov(particles, likelihood_samples_1d)), c(2,1))
  expect_equal(dim(cov(particles, likelihood_samples_2d)), c(2,2))
  expect_equal(dim(test_conditional_cov(particles, likelihood_samples_1d)), c(1, 1))
  expect_equal(dim(test_conditional_cov(particles, likelihood_samples_2d)), c(2, 2))
})

test_that('Covariance matrices are positive definite', {
  expect_equal(is.positive.definite(cov(likelihood_samples_1d)), TRUE)
  expect_equal(is.positive.definite(cov(likelihood_samples_2d)), TRUE)
  expect_equal(is.positive.definite(cov(particles)), TRUE)
  # Add test for conditional covariance matrix
})

######################## Testing For Adaptive Temperature ####################


test_that('ESS is between 1 and N', {
  expect_lte(estimate_ess(0.1, 0, simulated_data_1d, likelihood_samples_1d),  num_particles)
  expect_gte(estimate_ess(0.1, 0, simulated_data_1d, likelihood_samples_1d),  1)
})

test_that('Next selected temperature is always >= current_temp and <= 1', {
  expect_lte(find_next_temp(0.1, simulated_data_1d, likelihood_samples_1d, num_particles*0.5), 1)
  expect_gt(find_next_temp(0.1, simulated_data_1d, likelihood_samples_1d, num_particles*0.5), 0.1)
})

test_next_temp_valid <- function(current_temp, simulated_data, likelihood_samples, target_ess) {
  
  next_temp <- find_next_temp(current_temp, simulated_data, likelihood_samples, target_ess)
  return(estimate_ess(next_temp, current_temp, simulated_data, likelihood_samples))
}

# ToDo: generate variable likelihood samples so this test works better
test_that('Next selected temperature gives the correct ESS', {
  expect_equal(test_next_temp_valid(0.1, simulated_data_1d, likelihood_samples_1d, num_particles*0.5), num_particles*0.5)
})