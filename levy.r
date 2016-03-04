require('sde')
require('actuar')
require('ggplot2')
require('stabledist')
require('rootSolve')

plot_process <- function(time, jump_time, jump_size, X){
  # plots a sample path of the levy process
  index = order(jump_time)
  jump_time = sort(jump_time)
  jump_size = jump_size[index]
  n = length(X)
  tj = length(jump_time)
  X_new = rep(0, n + tj)
  time_new = rep(0, n + tj)
  X_group_ind = rep(0, n+tj)
  # extend X by repeating observations when jump happened
  offset = 0
  next_jump = jump_time[1 + offset]
  for (i in 1:n){
    X_new[i + offset] = X[i]
    time_new[i + offset] = time[i]
    X_group_ind[i+offset] = offset
    if (time[i] == next_jump){
      X_new[i + offset] = X_new[i + offset] - jump_size[1 + offset]
      offset = offset + 1
      X_new[i + offset] = X[i]
      time_new[i + offset] = time[i]
      X_group_ind[i + offset] = offset
      if (offset != tj-1){
        next_jump = jump_time[1 + offset]
      }
    }
  }
  df = data.frame(x = time_new, y = X_new, group = X_group_ind)
  p <- ggplot(df, aes(x = x, y = y, group = group))
  print(p+geom_line())
  return(df)
}

compare_distr <- function(distr_A, distr_B){
  # Kolmogorov-Smirnoff type of cdf comparison
  distr_A <- sort(distr_A)
  distr_B <- sort(distr_B)
  union <- sort(c(distr_A, distr_B))
  nA <- length(distr_A)
  nB <- length(distr_B)
  if (length(unique(union))!= nA + nB){
    print('resample and repeat, something went wrong')
    stop
  }
  dif_in_P <- rep(0, nA + nB)
  counterA = 0
  counterB = 0
  for (i in 1:(nA+nB)){
    if ((counterA != nA -1) && (distr_A[1 + counterA] == union[i])){
      counterA = counterA + 1
    }
    else{
      counterB = counterB + 1
    }
    dif_in_P[i] = abs(counterA/nA - counterB/nB)
  }
  return(list(dist = dif_in_P, x = union))
}

compare_density_plots<-function(distr_A, distr_B, xlim = c(-50, 50)){
  # plots density plots for two distributions on the same axis
  distr_A <- sort(distr_A)
  distr_B <- sort(distr_B)
  nA <- length(distr_A)
  nB <- length(distr_B)
  union <- matrix(c(distr_A, distr_B,rep(1, nA), rep(2, nB)), ncol = 2)
  union <- union[order(union[,1]),]
  df = data.frame(x = union[,1], group = union[,2])
  p <- ggplot(df, aes(x = x))
  print(p + geom_density(aes(group = group, colour = factor(group))) + xlim(xlim))
}

test_stable_package <- function(){
  # not incredibly useful, just to make sure that the package
  # from CRAN is doing what it is supposed to do
  simul_stable <- function(n, alpha){
    W <- rexp(n)
    Phi <- runif(n, min = -0.5 * pi, max = 0.5 * pi)
    S <- sin(alpha * Phi)/cos(Phi)^(1/alpha) * (cos(Phi * (1-alpha))/W)^((1-alpha)/alpha)
    return(S)
  }
  test_simulator <- function(sigma, num_samples, density = F){
    package <- rstable(num_samples, alpha = sigma, beta = 0)
    self <- simul_stable(num_samples, alpha = sigma)
    if (density)
      compare_density_plots(package, self)
    else{
      temp <- compare_distr(package, self)
      plot(temp$x, temp$dist, typ = 'l', xlim = c(-100,100))
    }
  }
  test_simulator(0.5, 10000, F)
  test_simulator(0.5, 10000, T)
}

stable_process <- function(t_in, t_fin, ngrid, epsilon, sigma, show_plot = F) {
  # Simulation of jumps
  M <- 2 * epsilon ^ (- sigma) / sigma
  N <- rpois(1, (t_fin-t_in) * M)
  time_jump<- runif(N, t_in, t_fin)
  delta_x <- (2*rbinom(N,1,0.5)-1) * rpareto(N, sigma, epsilon) / M
  
  time <- seq(t_in, t_fin, length.out = ngrid)
  time <- c(time_jump, time)
  index <- order(time)
  time <- sort(time)
  # mu_epsilon <- (1 - epsilon ^ (1 - sigma)) / (1-sigma)
  mu_epsilon <- 0
  X <- c(0, rep(mu_epsilon * (t_fin - t_in) / ngrid, ngrid - 1))
  X <- c(delta_x, X)
  X <- X[index]
  X <- cumsum(X)
  if (show_plot){
    plot(time, X)
    plot_process(time, time_jump,delta_x, X)
  }
  return(X[length(X)])
}

empirical_stable_distr <- function(epsilon, sigma, num_samples){
  simul_val <- rep(0, num_samples)
  for (i in 1:num_samples){
    simul_val[i] <- stable_process(0,1,1000,epsilon, sigma)
  }
  return(simul_val)
}

compare_stable <- function(epsilon, sigma, num_samples, density = F){
  empirical <- empirical_stable_distr(epsilon, sigma, num_samples)
  exact <- rstable(num_samples, alpha = sigma, beta = 0)
  if (density)
    compare_density_plots(empirical, exact)
  else{
    temp <- compare_distr(empirical, exact)
    plot(temp$x, temp$dist, typ = 'l', xlim = c(-100,100))
  }
}

# run some tests
compare_stable(0.01, 0.5, 10000, F)
compare_stable(0.01, 0.5, 10000, T)

gamma_process <- function(t_in, t_fin, ngrid, epsilon, tau, sigma) {
# Simulation of jumps
  density <- function(x)
    return (x^(-1-sigma) * exp(-tau * x))  
  M <- integrate(density, epsilon, Inf)[[1]]
  M2 <- epsilon^(-sigma)/sigma * exp(-tau * epsilon) - 
    tau/sigma * gamma(1-sigma)/tau^(1-sigma) * (1-pgamma(epsilon, 1-sigma, tau))
  print(M)
  print(M2)
  
  N <- rpois(1, (t_fin - t_in) * M)
  time_jump <- runif(N, t_in, t_fin)
  delta_x <- rep(0,N)
  for (i in 1:N) {
    x <- rpareto(1, sigma, epsilon)
    
    while (runif(1, 0, (epsilon^(-sigma) / sigma) * dpareto(x, sigma, epsilon)) > density(x))
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

gamma_process(0,1,1000,0.01,0.02,0.01)


density <- function(x)
    return (x^(-1-0.5) * exp(-0.5 * x))  
(M <- integrate(density, 0.01, Inf))
class(M)
as.double(M)
M[[1]]

mean_jump_dens <- function(x)
  return (x^(-0.5) * exp(-0.5 * x)) 
(mu_epsilon <- integrate(mean_jump_dens, 0.01, 1)[[1]])
