check_c_yy <- function(covariances, column) {
  
  return(mean(diag(covariances$C_yy)))
}

check_c_y_given_x <- function(covariances, column) {
  
  return(mean(diag(covariances$C_y_given_x)))
}

check_c_yx <- function(covariances, column) {
  
  return(mean(covariances$C_yx[, column]))
}

check_stepsize_var <- function(covariances, column) {
  
  stepsize_var <- ginv(covariances$C_yy)
  return(mean(diag(stepsize_var)))
}

check_stepsize <- function(covariances, column) {
  
  stepsize_var <- ginv(covariances$C_yy)
  stepsize <- stepsize_var %*% covariances$C_yx
  return(mean(stepsize[, column]))
}


