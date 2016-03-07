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
    print(length(unique(union)))
    print(nA + nB)
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
  nA <- length(distr_A)
  nB <- length(distr_B)
  df = data.frame(x = c(distr_A, distr_B), group = c(rep(1, nA), rep(2, nB)))
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
    return(2 ^(1/alpha) * S)
  }
  test_simulator <- function(sigma, num_samples, density = F){
    package <- rstable(num_samples, alpha = sigma, beta = 0, gamma = 2^(1/sigma))
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
test_stable_package()


#============================== STABLE PROCESS ==============================
stable_process <- function(t_in, t_fin, ngrid, epsilon, sigma, show_plot = F) {
  # Simulation of jumps
  M <- 2 * epsilon ^ (- sigma) / sigma
  N <- rpois(1, (t_fin-t_in) * M)
  time_jump<- runif(N, t_in, t_fin)
  delta_x <- (2*rbinom(N,1,0.5)-1) * (rpareto(N, sigma, epsilon)+epsilon)
  
  time <- seq(t_in, t_fin, length.out = ngrid)
  time <- c(time_jump, time)
  index <- order(time)
  time <- sort(time)
  # mu_epsilon <- (1 - epsilon ^ (1 - sigma)) / (1-sigma)
  mu_epsilon <- 0
  
# BM variance 
  sigma_2eps <- (2 / (2 - sigma)) * epsilon ^ (2 - sigma)
  W <- rep(0, ngrid + N)
  W[2:(ngrid + N)] <- cumsum(W[1:(ngrid + N - 1)] +
                             rnorm(ngrid + N - 1, 0,
                                   sqrt( sigma_2eps * (time[2:(ngrid+N)] - time[1:(ngrid+N-1)])))) 

  X <- c(0, rep(mu_epsilon * (t_fin - t_in) / ngrid, ngrid - 1))
  X <- c(delta_x, X)
  X <- X[index]
  X <- cumsum(X)
  X <- X + W
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
#==========================================================================

berry_esseen_stable <- function(epsilon, sigma){
  sigma_sq <- 2 * epsilon^(2-sigma)/(2-sigma)
  third_moment <- 2 * epsilon^(3-sigma)/(3-sigma)
  BE_bound <- 0.7975 * sigma_sq ^(-3/2) * third_moment
  return(BE_bound)
}

compare_stable_fixed_eps <- function(epsilon, sigma, num_samples, density = F){
  empirical <- empirical_stable_distr(epsilon, sigma, num_samples)
  exact <- rstable(num_samples, alpha = sigma, beta = 0, gamma = 2^(1/sigma))
  if (density)
    compare_density_plots(empirical, exact)
  else{
    temp <- compare_distr(empirical, exact)
    bound <- berry_esseen_stable(epsilon, sigma)
    plot(temp$x, temp$dist, typ = 'l', xlim = c(-100,100), ylim = c(0, max(c(1.1 * bound,temp$dist))))
    print(bound)
    abline(a = bound, b = 0, col = 'red')
  }
}

compare_stable <- function(epsilon_range, sigma, num_samples){
  n <- length(epsilon_range)
  max_discrep <- rep(0, n)
  BE_bound <- rep(0, n)
  
  for (i in seq_along(epsilon_range)){
    print(paste("Starting loop",i,"out of",n))
    eps <- epsilon_range[i]
    empirical <- empirical_stable_distr(eps, sigma, num_samples)
    exact <- rstable(num_samples, alpha = sigma, beta = 0)
    temp <- compare_distr(empirical, exact)
    max_discrep[i] <- max(temp$dist)
    BE_bound[i] <- berry_esseen_stable(eps, sigma)
  }
  
  plot(epsilon_range, BE_bound, type = 'l', col = 'red', ylim = c(0, BE_bound[1]))
  lines(epsilon_range, max_discrep)
}

# run some tests
compare_stable_fixed_eps(0.01, 0.5, 10000, F)
compare_stable_fixed_eps(0.01, 0.5, 10000, T)
?rpareto
compare_stable(seq(1,0.1, by = -0.1), 0.5, 1000)

simulate_gamma2 <- function(iter, sigma, tau, epsilon){
  # ----------------------------- #
  # disregard this abomination!!! #
  # ----------------------------- #
  x_star = -1/tau * log(sigma * epsilon ^ sigma)
  pareto <- function(x){
    dpareto(x,sigma, epsilon)
  }
  f <- function(x)
    return (x^(-1-sigma) * exp(-tau * x))
  M <- integrate(f, epsilon, Inf)[[1]]
#  f <- function(x)
#    return(density(x)/M)
  
  P = integrate(pareto, x_star, Inf)[[1]]
  MaxVal = epsilon ^(-1-sigma) * exp(-tau * epsilon)
  A = x_star * MaxVal

  total_accepted = 0
  samples <- rep(0, iter)
  while (total_accepted < iter){
    if (runif(1,0,1) < A/(P+A) ){
      y <- runif(1,0,x_star)
      if (runif(1,0, MaxVal) < f(y)){
        total_accepted = total_accepted + 1
        samples[total_accepted] = y
      }
    }
    else{
      y <- rpareto(1, sigma, epsilon)
      while (y < x_star)
        y <- rpareto(1, sigma, epsilon)
      if (runif(1,0, dpareto(y, sigma, epsilon)) < f(y)){
        total_accepted = total_accepted + 1
        samples[total_accepted] = y
      }
    }
  }
  return(samples)
}

simulate_gamma <- function(n, sigma, tau, epsilon){
  samples <- rep(0, n)
  for (i in 1:n){
    prop <- rpareto(1, sigma, epsilon)+epsilon
    while (runif(1) > exp(-tau * prop))
      prop <- rpareto(1, sigma, epsilon)+epsilon
    samples[i] <- prop
  }
  return(samples)
}

#============================== GAMMA PROCESS ==============================
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
  delta_x <- simulate_gamma(N, sigma, tau, epsilon)/ M

  time <- seq(t_in, t_fin, length.out = ngrid)
  time <- c(time_jump, time)
  index <- order(time)
  time <- sort(time)
  X <- c(0, rep(- (t_fin - t_in) / ngrid, ngrid - 1))
  mean_jump_dens <- function(x)
     return (x^(-sigma) * exp(-tau * x))  
  mu_epsilon <- integrate(mean_jump_dens, epsilon, 1)[[1]]

# BM variance 
  sigma_2eps <- pgamma(epsilon, tau, rate = 3 - sigma) * gamma(3 - sigma) / tau ^ (3 - sigma)
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
  plot(time, X)
  return (N)
}

gamma_process(0,1,1000,0.01,0.05,0.04)
#============================================================================

# pgamma(epsilon, 2-sigma, tau) * gamma(2-sigma) /tau^(2-sigma)


density <- function(x)
    return (x^(-1-0.5) * exp(-0.5 * x))  
(M <- integrate(density, 0.01, Inf))
class(M)
as.double(M)
M[[1]]

mean_jump_dens <- function(x)
  return (x^(-0.5) * exp(-0.5 * x)) 
(mu_epsilon <- integrate(mean_jump_dens, 0.01, 1)[[1]])
