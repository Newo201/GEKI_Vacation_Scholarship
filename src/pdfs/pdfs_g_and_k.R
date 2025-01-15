loglike_g_and_k <- function(y, parameters) {
  # Need to convert parameters back to the scale [0, 10]
  a <- pnorm(parameters$a)*10
  b <- pnorm(parameters$b)*10
  g <- pnorm(parameters$g)*10
  k <- pnorm(parameters$k)*10
  
  return(sum(dgk(y, a, b, g, k), log = T))
}