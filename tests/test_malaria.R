pacman::p_load(pacman, testthat, deSolve, rootSolve)

source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/models/eki_malaria.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/samples/samples_malaria.R')

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
  
  constrained_parameters <- list(d_in = 2/52, phi = 0.25, eta0 = 0.11, sigma = 10)
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

particles <- initialise_malaria_particles(num_particles, prior_params)
likelihood_samples <- synthetic_malaria(num_particles, particles, c())
sum(is.na(likelihood_samples))

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

constrained_parameters <- list(d_in = 0.001, phi = 0.25, eta0 = 0.05, sigma = 10)
l_mean <- likelihood_malaria_mean(constrained_parameters)

########################### Julia Output #####################

julia_params <- c(2.9203486e7, 66.67, 0.93, 1.6192924956625236, 0.057692307692307696, 
                  0.87, 0.08, 0.67, 1.23, 0.7026905137530658, 0.5074433108573343)

julia_data <- c(12.245185350950244, 12.109247278026844, 12.166246761399602, 12.457850971978136, 
                12.808273579774317, 13.026765817820838, 13.041064329759188, 12.887283673162376, 
                12.67043539502214, 12.478087223535917, 12.330680071567981, 12.207332612950452, 
                12.09036757079448, 12.015259800374738, 12.104482534704172, 12.41353642996914, 
                12.77406479582915, 12.998495486105837, 13.018864458287556, 12.870769735901215, 
                12.659600027523748, 12.471430972126813, 12.32697985770403, 12.204895934938465, 
                12.088855046360026, 12.01409978696837, 12.104105315061306, 12.413934924448888, 
                12.772292857858853, 12.998453925800433, 13.017930623553525, 12.871497740730065, 
                12.659232354536817, 12.4713244016031, 12.326900999896882, 12.204937874288847, 
                12.088715021270401, 12.014227991829804, 12.103790435911076, 12.413846370285114, 
                12.772694274948911, 12.998536803804969, 13.017708745852286, 12.8713846271384, 
                12.6591818223559, 12.471642669878056, 12.326725123922254, 12.2048337329414, 
                12.088905809650011, 12.014154180263143, 12.103833068857048, 12.412879193464994, 
                12.773634310952165, 12.998045541693442, 13.018608633047924, 12.87049685779888, 
                12.659420711466757, 12.471373173015538, 12.326921671352315, 12.204805598783151, 
                12.088918812455734, 12.014077118031853, 12.10414553654831, 12.413341134914038, 
                12.772870794693405, 12.998310513944407, 13.0180926817107, 12.871206936507802, 
                12.659295879065803, 12.47129925151161, 12.326922777476998, 12.204923467973448, 
                12.088735669367805, 12.014186065386692, 12.103877821691604, 12.413793277544357, 
                12.772667110714808, 12.998551763826116, 13.01768301038459, 12.87141833722966, 
                12.659194318399841, 12.4716379572572, 12.326718502054806, 12.204862969414142, 
                12.088873865826601, 12.014186890098054, 12.103776130552612, 12.412951601520998, 
                12.773613839400461, 12.998038400534718, 13.01865090652373, 12.870450676685493, 
                12.659399714906474, 12.471418418347975, 12.326903766157862, 12.204781183894276, 
                12.088956751146199, 12.014087725674882, 12.104126535235824, 12.413160777119764, 
                12.773068914149066, 12.998253637678346, 13.018176572530818, 12.871083719054935, 
                12.659323518076139, 12.471290284363548, 12.32693185045382, 12.204916094100208, 
                12.088747093873392, 12.014163101605428, 12.1039269283379, 12.41374729902032, 
                12.772667044569888, 12.998551800594182, 13.017672744931682, 12.87144057226024, 
                12.659201640482296, 12.471627399974375, 12.326714459266631, 12.204888837561082, 
                12.088844418965742, 12.014219214822338, 12.103726147893742, 12.413051212176144, 
                12.773562582291305, 12.998041518833945, 13.018659445466646, 12.870447739874509, 
                12.6592945111889)

constrained_parameters <- list(d_in = 1.6192924956625236, phi = 0.7026905137530658, 
                               eta0 = 0.5074433108573343, sigma = 10)
r_output <- log(diff(likelihood_malaria_mean(constrained_parameters)))

r_output - julia_data