source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/models/eki_normal.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/utils/tempering.R')

pacman::p_load(pacman, testthat, purrr, matrixcalc)

# I believe we need to make sure the number of particles > number of dimensions which makes logical sense
num_particles <- 400
num_dimensions <- 5

######################### Fixtures #######################################

# ToDo: see if there is a proper way of implementing fixtures like in Python
prior_params = list(alpha.mean = 0, alpha.sd = 5, sigma2.mean = 0, sigma2.sd = 2)
true_params_1d = list(alpha = 2, sigma = 2, x = 1)
true_params_2d = list(alpha = 2, sigma = 2, x = c(1, 2))
true_params_nd = list(alpha = 2, sigma = 2, x = rep(1, num_dimensions))

particles <- initialise_normal_particles(num_particles, prior_params)

likelihood_samples_1d <- synthetic_normal(num_particles, particles, true_params_1d)

true_data_1d <- likelihood_normal(true_params_1d)
simulated_data_1d <- matrix(true_data_1d, nrow = num_particles, ncol = 1, byrow = T)

likelihood_samples_2d <- synthetic_normal(num_particles, particles, true_params_2d)
true_data_2d <- likelihood_normal(true_params_2d)
simulated_data_2d <- matrix(true_data_2d, nrow = num_particles, ncol = 2, byrow = T)


likelihood_samples_nd <- synthetic_normal(num_particles, particles, true_params_nd)
true_data_nd <-likelihood_normal(true_params_nd)
simulated_data_nd <- matrix(true_data_nd, nrow = num_particles, ncol = num_dimensions, byrow = T)

covariances_1d <- calculate_covariances(particles, likelihood_samples_1d)
covariances_2d <- calculate_covariances(particles, likelihood_samples_2d)
covariances_nd <- calculate_covariances(particles, likelihood_samples_nd)

# covariances <- calculate_covariances(particles, likelihood_samples_1d)
# covariances$C_yy

######################### Testing For Dimensions ##############################

test_conditional_cov <- function(particles, likelihood_samples) {
  
  # Calculate the covariance matrices
  C_xx = cov(particles)
  # print(C_xx)
  C_yy = cov(likelihood_samples)
  C_xy = cov(particles, likelihood_samples)
  C_yx = cov(likelihood_samples, particles)
  
  C_y_given_x = C_yy - C_yx %*% solve(C_xx) %*% C_xy
  return(round(C_y_given_x, 5))
}

generate_pertubations <- function(num_particles, particles, likelihood_samples) {
  d_y <- dim(likelihood_samples)[2]
  C_y_given_x <- test_conditional_cov(particles, likelihood_samples)
  # Generate perturbations
  eta <- matrix(data = rmvnorm(n = num_particles, mean = rep(0, d_y), sigma = C_y_given_x), nrow = num_particles, ncol = d_y)
  return(eta)
}

test_that('Dimensions of particles are correct', {
  expect_equal(dim(particles), c(num_particles, 2))
  expect_equal(dim(particles), c(num_particles, 2))
})

test_that('Dimensions of weights are correct', {
  sum_of_sq <- get_sum_of_sq(simulated_data_1d, likelihood_samples_1d, covariances_1d, num_particles)
  expect_equal(length(get_weights(0.1, 0, sum_of_sq)), num_particles)
})

test_that('Dimensions of pertubations are correct', {
  expect_equal(dim(generate_pertubations(num_particles, particles, likelihood_samples_1d)), c(num_particles, 1))
  expect_equal(dim(generate_pertubations(num_particles, particles, likelihood_samples_2d)), c(num_particles, 2))
  expect_equal(dim(generate_pertubations(num_particles, particles, likelihood_samples_nd)), c(num_particles, num_dimensions))
})

test_that('Dimensions of particles are correct after updating', {
  expect_equal(dim(update_particles(0.1, particles, simulated_data_1d, likelihood_samples_1d, covariances_1d, num_particles)), dim(particles))
})



######################## Testing For Covariance Matrices #####################

test_that('Dimensions of covariance matrices are correct', {
  expect_equal(dim(cov(likelihood_samples_1d)), c(1, 1))
  expect_equal(dim(cov(likelihood_samples_2d)), c(2, 2))
  expect_equal(dim(cov(likelihood_samples_nd)), c(num_dimensions, num_dimensions))
  expect_equal(dim(cov(particles)), c(2,2))
  expect_equal(dim(cov(particles, likelihood_samples_1d)), c(2,1))
  expect_equal(dim(cov(particles, likelihood_samples_2d)), c(2,2))
  expect_equal(dim(cov(particles, likelihood_samples_nd)), c(2, num_dimensions))
  expect_equal(dim(test_conditional_cov(particles, likelihood_samples_1d)), c(1, 1))
  expect_equal(dim(test_conditional_cov(particles, likelihood_samples_2d)), c(2, 2))
  expect_equal(dim(test_conditional_cov(particles, likelihood_samples_nd)), c(num_dimensions, num_dimensions))
})

test_that('Covariance matrices are symmetric', {
  expect_equal(isSymmetric(cov(likelihood_samples_1d)), TRUE)
  expect_equal(isSymmetric(cov(likelihood_samples_2d)), TRUE)
  expect_equal(isSymmetric(cov(likelihood_samples_nd)), TRUE)
  expect_equal(isSymmetric(cov(particles)), TRUE)
  # Add test for conditional covariance matrix
  expect_equal(isSymmetric(test_conditional_cov(particles, likelihood_samples_nd)), TRUE)
})

test_that('Covariance matrices are positive definite', {
  expect_equal(is.positive.definite(cov(likelihood_samples_1d)), TRUE)
  expect_equal(is.positive.definite(cov(likelihood_samples_2d)), TRUE)
  expect_equal(is.positive.definite(cov(likelihood_samples_nd)), TRUE)
  expect_equal(is.positive.definite(cov(particles)), TRUE)
  expect_equal(is.positive.definite(test_conditional_cov(particles, likelihood_samples_nd)), TRUE)
})

test_that('Covariance matrices are positive semi definite', {
  expect_equal(is.positive.semi.definite(cov(likelihood_samples_1d)), TRUE)
  expect_equal(is.positive.semi.definite(cov(likelihood_samples_2d)), TRUE)
  expect_equal(is.positive.semi.definite(cov(likelihood_samples_nd)), TRUE)
  expect_equal(is.positive.semi.definite(cov(particles)), TRUE)
  expect_equal(is.positive.semi.definite(test_conditional_cov(particles, likelihood_samples_nd), tol = 1e-3), TRUE)
})

######################## Testing For Adaptive Temperature ####################

test_that('ESS is between 1 and N', {

  ll_densities = densities_normal(true_data_1d, num_particles, particles, true_params_1d)
  expect_lte(estimate_ess(0.1, 0, ll_densities),  num_particles)
  expect_gte(estimate_ess(0.1, 0, ll_densities),  1)

})

test_that('Next selected temperature is always >= current_temp and <= 1', {
  ll_densities = densities_normal(true_data_1d, num_particles, particles, true_params_1d)
  expect_lte(find_next_temp(0.1, ll_densities, num_particles*0.5), 1)
  expect_gt(find_next_temp(0.1, ll_densities, num_particles*0.5), 0.1)
})

test_next_temp_valid <- function(current_temp, ll_densities, target_ess) {
  
  next_temp <- find_next_temp(current_temp, ll_densities, target_ess)
  return(abs(get_ess_diff(next_temp, current_temp, ll_densities, target_ess)))

}

# ToDo: account for when next temperature is 1
test_that('Next selected temperature gives the correct ESS', {
  # expect_lt(test_next_temp_valid(0.1, simulated_data_1d, likelihood_samples_1d, covariances_1d, num_particles, num_particles*0.5), 0.5)
  # expect_lt(test_next_temp_valid(0.1, simulated_data_2d, likelihood_samples_2d, covariances_2d, num_particles, num_particles*0.5), 0.5)
  ll_densities = densities_normal(true_data_nd, num_particles, particles, true_params_nd)
  expect_lt(test_next_temp_valid(0.1, ll_densities, num_particles*0.5), 0.5)
})

test_that('Dimensions of weights are correct', {
  ll_densities = densities_normal(true_data_1d, num_particles, particles, true_params_1d)
  expect_equal(length(get_weights(0.1, 0, ll_densities)), num_particles)
})

