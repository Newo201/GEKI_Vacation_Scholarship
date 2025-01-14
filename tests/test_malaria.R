pacman::p_load(pacman, testthat)

source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/models/eki_malaria.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/samples/samples_malaria.R')

##################### Fixtures #########################

prior_params <- list(din.sd = 2,
                     phi.mean = 0, phi.sd = 1,
                     eta0.mean = 0, eta0.sd = 1,
                     sigma.mean = 10, sigma.sd = 4)

######################## Likelihood Mean #####################################

# Not really sure how we test for the correctness of the differential equation solution
## We can look at more informal intuition - e.g. when eta0 is lower dW/dt is lower so
## we expect the number of cases to be lower also

test_that('Solution to first equation gives correct parameters')

test_that('Length of likelihood mean is correct', {
  
  constrained_parameters <- list(d_in = 2/52, phi = 0.25, eta0 = 0.11, sigma = 10)
  unconstrained_parameters <- unconstrain_malaria_params(constrained_parameters)
  l_mean <- likelihood_malaria_mean(unconstrained_parameters)
  expect_equal(length(l_mean), 129)
})

######################### Initialising Particles ############################

test_that('Dimensions of particles are as expected', {
  
  particles <- initialise_malaria_particles(400, prior_params)
  expect_equal(dim(particles), c(400, 4))
  
})

######################### Other #############################################

test_that('Parameters are being constrained', {
  
  unconstrained_parameters <- list(d_in = -0.1, phi = -0.5, eta0 = 3, sigma = -1)
  constrained_parameters <- constrain_malaria_params(unconstrained_parameters)
  expect_gt(constrained_parameters$d_in, 0)
  expect_gt(constrained_parameters$phi, 0)
  expect_lt(constrained_parameters$phi, 1)
  expect_gt(constrained_parameters$eta0, 0)
  expect_lt(constrained_parameters$eta0, 1)
  expect_gt(constrained_parameters$sigma, 0)
  
})

test_that('Parameters are being unconstrained', {
  
  constrained_parameters <- list(d_in = 2/52, phi = 0.25, eta0 = 0.11, sigma = 10)
  unconstrained_parameters <- unconstrain_malaria_params(constrained_parameters)
  re_constrained_parameters <- constrain_malaria_params(unconstrained_parameters)
  expect_equal(constrained_parameters, re_constrained_parameters)
})