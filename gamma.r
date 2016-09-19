
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gamma_process <- function(t_in, t_fin, ngrid, epsilon, tau, sigma, n_paths=1) {
# Simulation of jumps
#  density <- function(x)
#    return (x^(-1-sigma) * exp(-tau * x))  
#  M <- integrate(density, epsilon, Inf)[[1]]
#  print(M)
#  N <- rpois(1, (t_fin - t_in) * M)
#  time_jump <- runif(N, t_in, t_fin)
#  delta_x <- simulate_gamma(N, sigma, tau, epsilon)

  Xs <- list()
  times <- list()
  time_jumps <- list()
  delta_xs <- list()

  for (i in 1:n_paths) {

# *********** Teh-Fav alg ***************
  density <- function (s)
    return (s^(-1-sigma) * exp(-tau * s)) 
  w <- function (s,t)
    return (t^(-1-sigma) * exp(-tau * s))
  W <- function (s,t)
    return (integrate(function(x) w(x,t), t, Inf)[[1]])
  W_inv <- function (r, t)
    return (t - (1 / tau) * log(1 - (r * tau) / (t ^ (- 1 - sigma) * exp(- tau * t))))

  delta_x <- c()
  t <- epsilon
  flag <- 0
  while (flag == 0){
    r <- rexp(1)
    if (r > W(Inf, t))
      flag <- 1
    else {
      t_prime <- W_inv(r, t)
      if (runif(1) <  density(t_prime) / w(t_prime, t))
        delta_x <- c(delta_x, t_prime)
      t <- t_prime
    }
  }

  N <- length(delta_x)
  time_jump <- runif(N, t_in, t_fin)
# ****************************************

  time <- seq(t_in, t_fin, length.out = ngrid)
  time <- c(time_jump, time)
  index <- order(time)
  time <- sort(time)
  X <- c(0, rep(- (t_fin - t_in) / ngrid, ngrid - 1))
  mean_jump_dens <- function(x)
     return (x^(-sigma) * exp(-tau * x))  
  mu_epsilon <- - integrate(mean_jump_dens, epsilon, 1)[[1]]
#  mu_epsilon <- 0
                                        # BM variance 
  sigma_2eps <- integrate(function(x) exp(-tau * x) * x ^ (1-sigma),0,epsilon)[[1]]
  print(sigma_2eps)
  W <- rep(0, ngrid + N)
  W[2:(ngrid + N)] <- cumsum(W[1:(ngrid + N - 1)] +
                             rnorm(ngrid + N - 1, 0,
                                   sqrt( sigma_2eps * (time[2:(ngrid+N)] - time[1:(ngrid+N-1)])))) 
  
  X <- c(0, rep( mu_epsilon * (t_fin - t_in) / ngrid, ngrid - 1))
#  X <- rep(0, ngrid)
  X <- c(delta_x, X)
  X <- X[index]
  X <- cumsum(X)
  X <- X + W
  #plot(time, X)

  Xs[[i]] <- X
  times[[i]] <- time 
  time_jumps[[i]] <- time_jump 
  delta_xs[[i]] <- delta_x 
}
  if(n_paths == 1)
    plot_process(time, time_jump, delta_x, X)
  else
    plot_many(times, time_jumps, delta_xs, Xs)
  return (N)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

