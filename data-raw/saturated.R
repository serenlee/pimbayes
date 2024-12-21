data_binary_multi_sat <- function(n,p, seed = 777) {
  set.seed(seed)
  px = c(gtools::rdirichlet(1,rep(1,2^p)))
  py.x = runif(2^p, 0, 1)
  pr.xy0 = runif(2^p, 0, 1)
  log.or.r = rnorm(2^p, sd = 1)
  pr.xy1 = expit(log.or.r + logit(pr.xy0))
  x = t(rmultinom(n,1,px))
  y = rbinom(n, 1, prob = x %*% py.x)
  r = vector("numeric",n)
  r[y == 1] = rbinom(sum(y == 1), 1, prob = x[y == 1,] %*% pr.xy1)
  r[y == 0] = rbinom(sum(y == 0), 1, prob = x[y == 0,] %*% pr.xy0)
  ymiss = y
  ymiss[r == 0] = NA
  return(list(x = x, y = ymiss, r = r))
}

small_sat<-data_binary_multi_sat(100,2)
sat<-data_binary_multi_sat(3000,4)


usethis::use_data(small_sat, overwrite = TRUE)
usethis::use_data(sat, overwrite = TRUE)
