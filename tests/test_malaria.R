pacman::p_load(pacman, testthat, deSolve, rootSolve)

source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/models/eki_malaria_d_in_only.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/samples/samples_malaria.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/utils/eki_helper.R')

##################### Fixtures #########################
num_particles <- 400

prior_params <- list(din.sd = 2,
                     phi.mean = 0, phi.sd = 1,
                     eta0.mean = 0, eta0.sd = 1,
                     sigma.mean = 10, sigma.sd = 4)

######################## Likelihood Mean #####################################

# Not really sure how we test for the correctness of the differential equation solution
## We can look at more informal intuition - e.g. when eta0 is lower dW/dt is lower so
## we expect the number of cases to be lower also

# test_that('Solution to first equation gives correct parameters')

test_that('Sum of variables equals population', {
  parameters <- c(N = 29203486,
                  L = 66.67,
                  dimm = 0.93,
                  d_in = 0.25, 
                  d_treat = 3/52,
                  p1 = 0.87,
                  p2 = 0.08, 
                  amp = 0.67,
                  R_m = 1.23,
                  phi = 0.11,
                  eta0 = 0.11)
  
  #The initial conditions for solving the ODE 
  ICs <- solve_steady_state(parameters)
  expect_equal(Reduce("+", ICs), parameters['N'][[1]])
})

test_that('Length of likelihood mean is correct', {
  
  constrained_parameters <- list(d_in = 0.3, phi = 0.25, eta0 = 0.05, sigma = 10)
  l_mean <- likelihood_malaria_mean(constrained_parameters)
  expect_equal(length(l_mean), 130)
})

test_that('Likelihood mean is monotonically increasing', {
  constrained_parameters <- list(d_in = 0.3, phi = 0.25, eta0 = 0.05, sigma = 10)
  l_mean <- likelihood_malaria_mean(constrained_parameters)
  # print(l_mean)
  expect_equal(all(diff(l_mean) >= 0), TRUE)
})

test_that('Length of likelihood samples is correct', {
  
  constrained_parameters <- list(d_in = 26/52, phi = 0.25, eta0 = 0.11, sigma = 10)
  unconstrained_parameters <- unconstrain_malaria_params(constrained_parameters)
  likelihood_sample <- likelihood_malaria(unconstrained_parameters)
  expect_equal(length(likelihood_sample), 129)
  
})

######################### Initialising Particles ############################

test_that('Dimensions of particles are as expected', {

  particles <- initialise_malaria_particles(num_particles, prior_params)
  expect_equal(dim(particles), c(num_particles, 4))
})

test_that('Dimensions of synthetic data is as expected', {

  particles <- initialise_malaria_particles(num_particles, prior_params)
  likelihood_samples <- synthetic_malaria(num_particles, particles, c())
  expect_equal(dim(likelihood_samples), c(num_particles, 129))
  expect_equal(sum(is.na(likelihood_samples)), 0)
})

######################### Other #############################################

test_that('Parameters are being constrained', {
  
  unconstrained_parameters <- list(d_in = -0.1, phi = -0.5, eta0 = 3, sigma = -1)
  constrained_parameters <- constrain_malaria_params(unconstrained_parameters)
  expect_gt(constrained_parameters$d_in, 0.16)
  expect_gt(constrained_parameters$phi, 0)
  expect_lt(constrained_parameters$phi, 1)
  expect_gt(constrained_parameters$eta0, 0)
  expect_lt(constrained_parameters$eta0, 1)
  expect_gt(constrained_parameters$sigma, 0)
  
})

test_that('Constraints are working at extreme values', {
  unconstrained_parameters <- list(d_in = -100, phi = -100, eta0 = -100, sigma = -100)
  constrained_parameters <- constrain_malaria_params(unconstrained_parameters)
  expect_equal(constrained_parameters$d_in, 0.16)
  expect_equal(constrained_parameters$phi, 0)
  expect_equal(constrained_parameters$eta0, 0)
  expect_equal(constrained_parameters$sigma, 0)
  
  unconstrained_parameters <- list(d_in = 100, phi = 100, eta0 = 100, sigma = -1)
  constrained_parameters <- constrain_malaria_params(unconstrained_parameters)
  expect_equal(constrained_parameters$eta0, 1)
  expect_equal(constrained_parameters$phi, 1)
})

test_that('Parameters are being unconstrained', {
  
  constrained_parameters <- list(d_in = 26/52, phi = 0.25, eta0 = 0.11, sigma = 10)
  unconstrained_parameters <- unconstrain_malaria_params(constrained_parameters)
  re_constrained_parameters <- constrain_malaria_params(unconstrained_parameters)
  print(unconstrained_parameters)
  print(re_constrained_parameters)
  expect_equal(constrained_parameters, re_constrained_parameters)
})

######################## Looking at Covariances ########################

source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/models/eki_malaria.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/samples/samples_malaria.R')

data_path = "C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/data/Malariah_data.rds"
true_data = log(readRDS(data_path))
simulated_data <- matrix(true_data, nrow = num_particles, ncol = 129, byrow = T)

true_parameters <- list(d_in = 0.5, phi = 0.5, eta0 = 0.001, sigma = 0.5)
true_unconstrained_parameters <- unconstrain_malaria_params(true_parameters)
true_constrained_parameters <- constrain_malaria_params(true_unconstrained_parameters)
initial_particles <- initialise_malaria_particles_d_in_only(num_particles, prior_params)
likelihood_samples <- synthetic_malaria_d_in_only(num_particles, initial_particles, true_unconstrained_parameters)
covariances <- calculate_covariances(initial_particles, likelihood_samples)

covariances$C_xx
covariances$C_yy

(simulated_data - likelihood_samples)
ginv(covariances$C_yy + 9*covariances$C_y_given_x_inv) %*% covariances$C_yx

(simulated_data - likelihood_samples) %*% ginv(covariances$C_yy + 9*covariances$C_y_given_x_inv) %*% covariances$C_yx
