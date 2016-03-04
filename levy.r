
require('sde')
require('actuar')

stable_process <- function(t_in, t_fin, ngrid, epsilon, sigma) {
# Simulation of jumps
  M <- 2 * epsilon ^ (- sigma) / sigma
  N <- rpois(1, (t_fin-t_in) * M)
  time_jump<- runif(N, t_in, t_fin)
  delta_x <- (2*rbinom(N,1,0.5)-1) * rpareto(N, sigma, epsilon)

  time <- seq(t_in, t_fin, length.out = ngrid)
  time <- c(time_jump, time)
  index <- order(time)
  time <- sort(time)
  # mu_epsilon <- (1 - epsilon ^ (1 - sigma)) / (1-sigma)
  mu_epsilon <- 0
  X <- c(0, rep(- mu_epsilon * (t_fin - t_in) / ngrid, ngrid - 1))
  X <- c(delta_x, X)
  X <- X[index]
  X <- cumsum(X)
  plot(time, X)
  return (N)
}

stable_process(0,5,1000,0.01,0.5)



gamma_process <- function(t_in, t_fin, ngrid, epsilon, tau, sigma) {
# Simulation of jumps
  density <- function(x)
    return (x^(-1-sigma) * exp(-tau * x))  
  M <- integrate(density, epsilon, Inf)[[1]]
  
  N <- rpois(1, (t_fin - t_in) * M)
  time_jump <- runif(N, t_in, t_fin)
  delta_x <- rep(0,N)
  for (i in 1:N) {
    x <- rpareto(1, sigma, epsilon)
    while (runif(1, 0, (epsilon^(-sigma) / sigma) * dpareto(x, sigma, epsilon) > density(x)))
      x <- rpareto(1, sigma, epsilon)
    delta_x[i] <- x / M
  }

  time <- seq(t_in, t_fin, length.out = ngrid)
  time <- c(time_jump, time)
  index <- order(time)
  time <- sort(time)
  X <- c(0, rep(- (t_fin - t_in) / ngrid, ngrid - 1))
  mean_jump_dens <- function(x)
     return (x^(-sigma) * exp(-tau * x))  
  mu_epsilon <- integrate(mean_jump_dens, epsilon, 1)[[1]]

  X <- c(0, rep( mu_epsilon * (t_fin - t_in) / ngrid, ngrid - 1))
#  X <- rep(0, ngrid)
  X <- c(delta_x, X)
  X <- X[index]
  X <- cumsum(X)
  plot(time, X)
  return (N)
}

gamma_process(0,1,1000,0.01,0.5,0.5)


density <- function(x)
    return (x^(-1-0.5) * exp(-0.5 * x))  
(M <- integrate(density, 0.01, Inf))
class(M)
as.double(M)
M[[1]]

mean_jump_dens <- function(x)
  return (x^(-0.5) * exp(-0.5 * x)) 
(mu_epsilon <- integrate(mean_jump_dens, 0.01, 1)[[1]])
