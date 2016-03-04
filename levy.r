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
  M <- tau ^ (-sigma) * gamma(sigma) * (1 - dgamma(tau * epsilon, sigma))
  N <- rpois(1, t_fin * M)
  time <- runif(N, t_in, t_fin)
  delta_x < <- rep(0,N)
  for (i in 1:N) {
    x <- rgamma(1, shape = sigma, rate = tau)
    while (x < epsilon)
      x <- rgamma(1, shape = sigma, rate = tau)
    delta[i] <- x
  }

  time <- seq(t_in, t_fin, length.out = ngrid)
  time <- c(time_jump, time)
  index <- order(time)
  time <- sort(time)
  X <- c(0, rep(- (t_fin - t_in) / ngrid, ngrid - 1))
  mu_epsilon <- (1 - epsilon ^ (1 - sigma)) / (1-sigma)
  X <- c(delta_x, X)
  X <- X[index]
  X <- cumsum(X)
  plot(time, X)
  return (N)
}

gamma_process(0,5,1000,0.01,0.5)
