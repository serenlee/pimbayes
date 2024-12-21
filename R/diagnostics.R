#' Diagnosing the identification of the model
#'
#' Estimate parameters using small posterior samples from JAGS,
#' calculate the approximate Fisher information matrix
#' and return the condition number and rank of the approximate Fisher information matrix
#' with asymptotic variances of estimated parameters.
#'
#' @param jags_code string. JAGS model code
#' @param data list. Data list of either vector or number
#' @param variable.names vector of string. Names of the variables
#' @param ll_fun function. Log-likelihood function that takes
#' two arguments - data and theta (parameters)
#' @param data_fun function. Function that generates the dataset
#' using theta (parameters).
#' argument of this function is
#' a list of arguments required for creating the dataset.
#'
#' @return list. A list of the condition number, rank, and asymptotic variances.
#' @examples
#' ll_cross <- function(data, theta){
#' se = theta[1]
#' sp = theta[2]
#' pi = theta[3:5]
#' pi = c(pi, 1-sum(pi))
#' eta =  c(se*pi[1]+(1-sp)*pi[3],
#'   se*pi[2]+(1-sp)*pi[4],
#'   sp*pi[3]+(1-se)*pi[1],
#'   (1-se)*pi[2]+sp*pi[4])
#'   log(dmultinom(data, prob = eta))
#'   }
#'
#'   data_cross <- function(args) {
#'     theta = args$theta
#'     n = args$n
#'     se = theta[1]
#'     sp = theta[2]
#'     pi = theta[3:5]
#'     pi = c(pi, 1-sum(pi))
#'     eta =  c(se*pi[1]+(1-sp)*pi[3],
#'       se*pi[2]+(1-sp)*pi[4],
#'       sp*pi[3]+(1-se)*pi[1],
#'       (1-se)*pi[2]+sp*pi[4])
#'       as.vector(rmultinom(1, size = n, prob = eta))
#'      }
#' jags_code <-
#' "model{
#' y[1:4] ~ dmulti(eta[1:4], n)
#' eta[1:4] <- c(se*pi[1]+(1-sp)*pi[3],
#' se*pi[2]+(1-sp)*pi[4],
#' sp*pi[3]+(1-se)*pi[1],
#' (1-se)*pi[2]+sp*pi[4])
#' pi[1:4] ~ ddirch(c(1,1,1,1))
#' se ~ dbeta(25,3)
#' sp ~ dbeta(30,1.5)
#' }"
#' data(small_cross)
#' data = list(y = small_cross, n = sum(small_cross))
#' diagnose_identification(jags_code, data, c("se","sp","pi[1]","pi[2]","pi[3]"),ll_cross,data_cross)
#'
#' @export
diagnose_identification<-function(jags_code, data, variable.names, ll_fun, data_fun) {
  mod <- jags.model(
    textConnection(jags_code),
    data = data,
    n.chains = 1,
    quiet = T
  )
  update(mod, 1000)
  results <- coda.samples(
    mod,
    n.iter = 10000,
    variable.names = variable.names
  )
  theta = apply(results[[1]], 2, median)[variable.names]
  fisher = (find_fisher(20, 2000, ll_cross, data_cross, c(list(theta = theta), data)))
  cond = pracma::cond(fisher)
  rank = pracma::Rank(fisher)
  list(condition_number = cond, rank = rank, var = 1/diag(fisher))
}

ll_cross <- function(data, theta){
  se = theta[1]
  sp = theta[2]
  pi = theta[3:5]
  pi = c(pi, 1-sum(pi))
  eta =  c(se*pi[1]+(1-sp)*pi[3],
                    se*pi[2]+(1-sp)*pi[4],
                   sp*pi[3]+(1-se)*pi[1],
                 (1-se)*pi[2]+sp*pi[4])
  log(dmultinom(data, prob = eta))
}

data_cross <- function(args) {
  theta = args$theta
  n = args$n
  se = theta[1]
  sp = theta[2]
  pi = theta[3:5]
  pi = c(pi, 1-sum(pi))
  eta =  c(se*pi[1]+(1-sp)*pi[3],
           se*pi[2]+(1-sp)*pi[4],
           sp*pi[3]+(1-se)*pi[1],
           (1-se)*pi[2]+sp*pi[4])
  as.vector(rmultinom(1, size = n, prob = eta))
}


find_fisher<-function(M = 1, N, ll_fun, data_fun, args) {
  p = length(args$theta)
  result <- matrix(0, nrow = p, ncol = p)
  for(i in 1:N) {
    result = result + hessian_i(i, M, ll_fun, p, data_fun, args)
  }
  -result/N
}



hessian_i <- function(i, M, ll_fun, p, data_fun, args) {
  data <- data_fun(args)
  result <- matrix(0, nrow = p, ncol = p)
  for(k in 1:M) {
    result = result + hessian(k, p, data, ll_fun, args)
  }
  result/M
}

hessian <- function(k, p, data, ll_fun, args) {
  deltak <- delta(k, p)
  deltaGk <- deltaG(k, deltak, p, data, ll_fun, args)
  inverse_deltak <- matrix(1/deltak, 1, p)
  temp = deltaGk %*% inverse_deltak
  0.5*(0.5 * temp + t(0.5 *temp))
}

deltaG <- function(k, deltak, p, data, ll_fun, args){
  deltaTildek <- delta(k, p)
  g(deltak, deltaTildek, T, p, data, ll_fun, args) - g(deltak, deltaTildek, F, p, data, ll_fun, args)
}

g <- function(deltak, deltaTildek, plus, p, data, ll_fun, args) {
  sign <- ifelse(plus, 1, -1)
  theta = args$theta
  lldiff <- ll_fun(data = data, theta = theta + sign * deltak + deltaTildek) -
    ll_fun(data = data, theta = theta + sign*deltak - deltaTildek)
  matrix(0.5 * lldiff / deltaTildek, nrow = p, ncol = 1)
}

delta<- function(k, p, c = 0.001) {
  bern<-rbinom(p,1,0.5)
  bern[bern == 0] = -c
  bern[bern == 1] = c
  bern
}



