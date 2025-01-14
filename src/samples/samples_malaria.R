pacman::p_load(pacman, extraDistr)

source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/utils/equations_malaria.R')

#################### Parameter Transformations #############################

constrain_malaria_params <- function(parameters) {
  
  d_in <- exp(parameters$d_in)
  phi <- plogis(parameters$phi) # Logistic transformation
  eta0 <- plogis(parameters$eta0) # Logistic transformation
  sigma <- exp(parameters$sigma)
  
  return(list(d_in = d_in, phi = phi, eta0 = eta0, sigma = sigma))
}

unconstrain_malaria_params <- function(parameters) {
  
  d_in <- log(parameters$d_in)
  phi <- qlogis(parameters$phi) # Logit transformation
  eta0 <- qlogis(parameters$eta0) # Logit transformation
  sigma <- log(parameters$sigma)
  
  return(list(d_in = d_in, phi = phi, eta0 = eta0, sigma = sigma))
  
}

####################### Prior Distributions #######################

# Generate samples from the prior distributions
# I think this might be stored in log space, but I have to check
log_din_prior_sample <- function(parameters, num_samples) {
  sd = parameters$din.sd # 2
  
  # Sample from a half-normal distribution and then apply a log transformation
  # so that the parameter is unconstrained
  return(log(rhnorm(num_samples, sigma = sd)))
}

logit_phi_prior_sample <- function(parameters, num_samples) {
  mean = parameters$phi.mean # 0
  sd = parameters$phi.sd # 1
  return(rnorm(num_samples, mean = mean, sd = sd))
}

logit_eta0_prior_sample <- function(parameters, num_samples) {
  mean = parameters$eta0.mean # 0
  sd = parameters$eta0.sd # 1
  return(rnorm(num_samples, mean = mean, sd = sd))
}

log_sigma_prior_sample <- function(parameters, num_samples) {
  mean = parameters$sigma.mean # 10
  sd = parameters$sigma.sd # 4
  
  return(rnorm(num_samples, mean = mean, sd = sd))
}

###################### Likelihood ##############################

solve_steady_state <- function(parameters) {
  
  I1 = 5
  I2 =10
  S = 29203486-I1-I2
  R = 0
  ICss <- c(S = S, I1 = I1, I2 =I2, R = R)     #init vals of latent states
  
  # Find the Equlibrum status___________________________________
  Eq<- runsteady(y = ICss, times = c(0,2000), func = mtdrift_theta,  parms = parameters)
  
  #store the steady state values________________________________
  S_1 <- Eq$y['S'][[1]]
  I1_1 <- Eq$y['I1'][[1]]
  I2_1 <- Eq$y['I2'][[1]]
  R_1 <- Eq$y['R'][[1]]
  
  # Solve Baseline Model__________________________________________
  #The initial conditions for solving the ODE   
  ICs<- c(S=S_1, I1= I1_1,  I2=I2_1 ,R= R_1,W=0)
  
  return(ICs)
  
}

likelihood_malaria_mean <- function(variable_parameters) {
  
  # Assumes variable parameters input is in constrained form
  
  parameters <- c(N = 29203486,
                      L = 66.67,
                      dimm = 0.93,
                      d_in = variable_parameters$d_in, 
                      d_treat0 = 3/52,
                      p1 = 0.87,
                      p2 = 0.08, 
                      amp = 0.67,
                      R_m = 1.23,
                      phi = variable_parameters$phi,
                      eta0 = variable_parameters$eta0,
                      k = 0.01, 
                      Tau = 17.33333)
  
  start_res=10                                 #Year Resistance Begins
  dt= 1/12                                     #step size one monthe per a year
  true_time <- seq(0, 10.67, by = dt)                  #time step
  
  #The initial conditions for solving the ODE 
  ICs <- solve_steady_state(parameters)
  
  #Solve the ODE using the Defulet Solver LSoda
  out <- (ode(y = ICs, times = true_time, func = mtdrift, parms = parameters))
  #Save first simulation;
  simulation.data<-c(NA)
  simulation.data[1]<-out[,"W"][1]
  for(k in 2:129){
    simulation.data[k]<-out[,"W"][k]-out[,"W"][k-1]
  }
  return(simulation.data)
}

likelihood_malaria <- function(variable_parameters) {
  
  # Assumes parameters are unconstrained
  variable_parameters <- constrain_malaria_params(variable_parameters)
  sigma <- variable_parameters$sigma
  likelihood_mean <- likelihood_malaria_mean(variable_parameters)
  
  return(rnorm(n = length(likelihood_mean), mean = likelihood_mean, sd = sigma))
}