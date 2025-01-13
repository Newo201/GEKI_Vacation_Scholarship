source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/models/eki_normal.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/models/eki_normal_known_var.R')

pacman::p_load(pacman, testthat, purrr, matrixcalc)

# I believe we need to make sure the number of particles > number of dimensions which makes logical sense
num_particles <- 1000
num_dimensions <- 5

######################### Fixtures #######################################

# ToDo: see if there is a proper way of implementing fixtures like in Python
prior_params = list(alpha.mean = 0, alpha.sd = 5, sigma2.mean = 0, sigma2.sd = 2)
true_params_1d = list(alpha = 2, sigma = 2, x = 1)
true_params_2d = list(alpha = 2, sigma = 2, x = c(1, 2))
true_params_nd = list(alpha = 2, sigma = 2, x = rep(1, num_dimensions))

# normal_particles <- initialise_normal_particles(num_particles, prior_params)
normal_particles_known_var <- initialise_normal_particles_known_var(num_particles, prior_params)

likelihood_samples_1d <- synthetic_normal_known_var(num_particles, normal_particles_known_var, true_params_1d)
simulated_data_1d <- matrix(likelihood_normal(true_params_1d), nrow = num_particles, ncol = 1, byrow = T)

likelihood_samples_2d <- synthetic_normal_known_var(num_particles, normal_particles_known_var, true_params_2d)
simulated_data_2d <- matrix(likelihood_normal(true_params_2d), nrow = num_particles, ncol = 2, byrow = T)

likelihood_samples_nd <- synthetic_normal_known_var(num_particles, normal_particles_known_var, true_params_nd)
simulated_data_nd <- matrix(likelihood_normal(true_params_nd), nrow = num_particles, ncol = num_dimensions, byrow = T)

covariances_1d <- calculate_covariances(normal_particles_known_var, likelihood_samples_1d)
covariances_2d <- calculate_covariances(normal_particles_known_var, likelihood_samples_2d)
covariances_nd <- calculate_covariances(normal_particles_known_var, likelihood_samples_nd)

######################## Testing For Covariance Matrices #####################

test_that('Dimensions of covariance matrices are correct', {
  expect_equal(dim(cov(likelihood_samples_1d)), c(1, 1))
  expect_equal(dim(cov(likelihood_samples_2d)), c(2, 2))
  expect_equal(dim(cov(likelihood_samples_nd)), c(num_dimensions, num_dimensions))
  expect_equal(dim(cov(normal_particles_known_var)), c(2,2))
  expect_equal(dim(cov(normal_particles_known_var, likelihood_samples_1d)), c(2,1))
  expect_equal(dim(cov(normal_particles_known_var, likelihood_samples_2d)), c(2,2))
  expect_equal(dim(cov(normal_particles_known_var, likelihood_samples_nd)), c(2, num_dimensions))
  expect_equal(dim(test_conditional_cov(normal_particles_known_var, likelihood_samples_1d)), c(1, 1))
  expect_equal(dim(test_conditional_cov(normal_particles_known_var, likelihood_samples_2d)), c(2, 2))
  expect_equal(dim(test_conditional_cov(normal_particles_known_var, likelihood_samples_nd)), c(num_dimensions, num_dimensions))
})

test_that('Covariance matrices are symmetric', {
  expect_equal(isSymmetric(cov(likelihood_samples_1d)), TRUE)
  expect_equal(isSymmetric(cov(likelihood_samples_2d)), TRUE)
  expect_equal(isSymmetric(cov(likelihood_samples_nd)), TRUE)
  expect_equal(isSymmetric(cov(normal_particles_known_var)), TRUE)
  # Add test for conditional covariance matrix
  expect_equal(isSymmetric(test_conditional_cov(normal_particles_known_var, likelihood_samples_nd)), TRUE)
})

test_that('Covariance matrices are positive definite', {
  expect_equal(is.positive.definite(cov(likelihood_samples_1d)), TRUE)
  expect_equal(is.positive.definite(cov(likelihood_samples_2d)), TRUE)
  expect_equal(is.positive.definite(cov(likelihood_samples_nd)), TRUE)
  expect_equal(is.positive.definite(cov(normal_particles_known_var)), TRUE)
  expect_equal(is.positive.definite(test_conditional_cov(normal_particles_known_var, likelihood_samples_nd)), TRUE)
})

test_that('Covariance matrices are positive semi definite', {
  expect_equal(is.positive.semi.definite(cov(likelihood_samples_1d)), TRUE)
  expect_equal(is.positive.semi.definite(cov(likelihood_samples_2d)), TRUE)
  expect_equal(is.positive.semi.definite(cov(likelihood_samples_nd)), TRUE)
  expect_equal(is.positive.semi.definite(cov(normal_particles_known_var)), TRUE)
  expect_equal(is.positive.semi.definite(test_conditional_cov(normal_particles_known_var, likelihood_samples_nd), tol = 1e-3), TRUE)
})

########################## Large Sample Limit #############################
test_that('Large sample limits of covariance matrices are as expected', {
  covariances <- calculate_covariances(normal_particles_known_var, likelihood_samples_nd)
  Q <- prior_params$alpha.sd**2 * diag(1)
  H <- true_params_nd$x
  R <- true_params_nd$sigma**2 * diag(num_dimensions)
  
  expect_lt(abs(covariances$C_xx - Q), 1)
  expect_lt(sum(abs(covariances$C_xy - Q %*% t(H))), num_dimensions * 1)
  expect_lt(sum(abs(covariances$C_yy - (H %*% Q %*% t(H) + R))), num_dimensions**2 * 1)
  expect_lt(sum(abs(covariances$C_y_given_x - R)), num_dimensions*1)
})

sum(abs(covariances_nd$C_yy - (H %*% Q %*% t(H) + R)))

cov(normal_particles_known_var)
cov(normal_particles_known_var, likelihood_samples_nd)
cov(likelihood_samples_nd)

true_params_nd$x %*% t(true_params_nd$x)

Q <- prior_params$alpha.sd**2 * diag(1)
H <- true_params_nd$x
R <- true_params_nd$sigma**2 * diag(num_dimensions)

Q %*% t(H)

covariances_nd$C_xy

t(H)

t(H) %*% Q