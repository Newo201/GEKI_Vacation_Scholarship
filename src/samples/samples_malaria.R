#################### Parameter Transformations #############################

constrain_malaria_params <- function(parameters) {
  
  d_in <- exp(parameters$d_in) + 0.16 # Transformation to [0.16, infinity)
  phi <- plogis(parameters$phi) # Logistic transformation to [0, 1]
  eta0 <- plogis(parameters$eta0) # Logistic transformation to [0, 1]
  sigma <- exp(parameters$sigma)
  
  return(list(d_in = d_in, phi = phi, eta0 = eta0, sigma = sigma))
}

unconstrain_malaria_params <- function(parameters) {
  
  d_in <- log(parameters$d_in - 0.16)
  phi <- qlogis(parameters$phi) # Logit transformation
  eta0 <- qlogis(parameters$eta0) # Logit transformation
  sigma <- log(parameters$sigma)
  
  return(list(d_in = d_in, phi = phi, eta0 = eta0, sigma = sigma))
  
}

####################### Prior Distributions #######################

# Generate samples from the prior distributions
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
  print(c(mean, sd))
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
                      d_treat = 3/52,
                      p1 = 0.87,
                      p2 = 0.08, 
                      amp = 0.67,
                      R_m = 1.23,
                      phi = variable_parameters$phi,
                      eta0 = variable_parameters$eta0)
  
  dt= 1/12 #step size one month per a year
  true_time <- seq(0, 10.75, by = dt)                  #time step
  
  #The initial conditions for solving the ODE 
  ICs <- solve_steady_state(parameters)
  
  #Solve the ODE using the Default Solver LSoda
  out <- (ode(y = ICs, times = true_time, func = mtdrift, parms = parameters))
  simulation.data <- out[,"W"]
  return(simulation.data)
}

likelihood_malaria <- function(variable_parameters, log_obs = T) {
  
  # Assumes parameters are unconstrained
  variable_parameters <- constrain_malaria_params(variable_parameters)
  # print(variable_parameters)
  sigma <- variable_parameters$sigma
  # Convert mean to log difference
  likelihood_mean <- log(diff(likelihood_malaria_mean(variable_parameters)))
  # Determine if data is stored in log space
  if (log_obs) {
    return(rnorm(n = length(likelihood_mean), mean = likelihood_mean, sd = sigma))
  } else {
    return(rlnorm(n = length(likelihood_mean), mean = likelihood_mean, sd = sigma))
  }

}


