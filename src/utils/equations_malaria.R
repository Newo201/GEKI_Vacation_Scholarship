##Ground Truth
mtdrift_theta <- function(t, ICss, parameters) {
  
  with(as.list(c(ICss, parameters)),{
    
    # Not sure what the purpose of this P is since we don't use it elsewhere
    # P = eta0*(p1*I1 + p2*I2)/(I1 + I2)
    
    a1 = eta0*p1*(1/d_treat) + (1-eta0*p1)*(1/d_in)
    a2 = eta0*p2*(1/d_treat) + (1-eta0*p2)*(1/d_in)
    
    # a1 = eta0*p1*(1/d_treat)
    # a2 = eta0*p2*(1/d_treat)

    # N is used instead of P to represent population size
    lambda_val = R_m*((1/L)+a1)*(I1 + I2)/N  # for Steady state calculation
    # print("Hello")
    # print(lambda_val)
    
    dS =  N/L - (lambda_val + (1/L))*S + (1/dimm)*R
    # print(S)
    # print(dS)
    
    dI1 = lambda_val*S - (a1 + (1/L))*I1
    # print(dI1)
    
    dI2 = lambda_val*R - (a2 + (1/L))*I2
    # print(dI2)
    
    dR = a1*I1 + a2*I2 - (lambda_val +(1/dimm) + (1/L))*R
    # print(dR)
    
    der <-c(dS, dI1, dI2,dR)
    # return the rate of change
    return(list(der))   
  
  })
}



mtdrift <-function(t, ICs, parameters) {
  
  with(as.list(c(ICs, parameters)),{
    
    ############ Resistance - not sure if needed #########################
    
    # #Switch Resistance on at time, start_res
    # if (t < start_res){
    #   res_on = 0 }
    # else
    # {res_on = 1}
    # 
    # #Start ppres0 at time, start_res
    # if (t < start_res+dt)
    # {ppres0 = 1}
    # else
    # {ppres0 = 100*(Tau-(1/d_treat0))/((1/d_in)-(1/d_treat0))}     
    
    #################################################################
    
    # C = p1*I1 + p2*I2
    # P = eta0*(p1*I1 + p2*I2)/(I1 + I2)
    
    a1 = eta0*p1*(1/d_treat) + (1-eta0*p1)*(1/d_in)
    a2 = eta0*p2*(1/d_treat) + (1-eta0*p2)*(1/d_in)
    
    # a1 = eta0*p1*(1/d_treat)
    # a2 = eta0*p2*(1/d_treat)
    
    lambda_val = (amp*cos(2*pi*(t-phi))+R_m)*(a1 + (1/L))*(I1+I2)/N # for actual simulation
    
    dS =  N/L - (lambda_val + 1/L)*S + 1/dimm*R
    dI1 = lambda_val*S - (a1 + (1/L))*I1
    dI2 = lambda_val*R - (a2 + (1/L))*I2
    dR = a1*I1 + a2*I2 - (lambda_val + 1/dimm + 1/L)*R
    dW=lambda_val*S*eta0*p1+lambda_val*R*eta0*p2

    der <-c(dS, dI1, dI2, dR, dW)
    # return the rate of change
    return(list(der))
    
  }) # end with(as.list ...
}