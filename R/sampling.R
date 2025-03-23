#' Posterior sampling of given JAGS (BUGS) code and data
#'
#' This function returns the results of sampling from posteriors of given model
#' and data using either importance sampling with
#' transparent reparameterization (ISTP) or pigeons algorithm.
#'
#' @param jags_code A string of JAGS code that specifies the model.
#' @param data A list of dataset.
#' @param type one of "sat", "cross", "pigeons"
#' @param extract_variable_names (optional) A vector of strings that
#' contain the variable names that will be extracted
#' @param n_chains n_chains for Pigeons.jl The size of final sample is
#' n_chains * 2^(n_rounds)
#' @param n_rounds n_rounds for Pigeons.jl The size of final sample is
#' n_chains * 2^(n_rounds)
#' @param N type sat and cross will generate the posterior sample of size N
#' @param extended_traces logical. If TRUE, all serialized chains are recorded.
#' @return Either a list or julia object. Results including samples.
#' @examples
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
#' results = sampling(jags_code, data, "cross", N = 15000)
#' results$ess
#' head(results$samples)
#'
#' data(small_sat)
#' sigma = 0.5 #Hyperparameter, prior variance of log odds ratios
#' data = c(small_sat, list(n = length(small_sat$r),
#'   p2 = dim(small_sat$x)[2],
#'   tau = 1/sigma^2))
#' jags_code  <- "model{
#'   for (i in 1:n) {
#'     x[i, 1:p2] ~ dmulti(px[1:p2],1)
#'     y[i] ~ dbern(inprod(x[i,1:p2], py_x[1:p2]))
#'     r[i] ~ dbern(pr[i])
#'     pr[i] <- (y[i]) * inprod(x[i,1:p2], pr_xy1[1:p2]) + abs(y[i]-1) * inprod(x[i,1:p2], pr_xy0[1:p2])
#'   }
#'   for(i in 1:p2){
#'     pr_xy1[i] <- ilogit(log_or_r[i] + logit(pr_xy0[i]))
#'   }
#'   px[1:p2] ~ ddirich(Int64wrapper(ones(p2)))
#'   for(i in 1:p2) {
#'     py_x[i] ~ dunif(0, 1)
#'     pr_xy0[i] ~ dunif(0, 1)
#'     log_or_r[i] ~ dnorm(0, tau)
#'   }
#' }"
#' results = sampling(jags_code, data, "sat", N = 15000)
#'
#' results$ess
#' head(results$samples)
#'
#'
#' jags_code <- "model{
#' for (i in 1:n) {
#' r[i] ~ dbern(pr[i])
#' pr[i] <- ilogit(y[i] * alpha1 + alpha0)
#' y[i] ~ dpois(mu)
#' }
#' mu ~ dgamma(1,1)
#' alpha0 ~ dnorm(0, 0.1)
#' alpha1 ~ dnorm(0, tau)
#' }
#' "
#' data(small_count)
#' sigma = 0.5 #Hyperparameter, prior variance of alpha1
#' data = c(small_count, list(n = length(small_count$r), tau = 1/sigma^2))
#' sampling(jags_code, data, "pigeons", n_chains =  4, n_rounds = 4)
#' @export
sampling<-function(jags_code, data, type = c("sat", "cross", "pigeons"), extract_variable_names = NULL,
                   n_chains  = NULL, n_rounds = NULL, N = NULL, extended_traces = F){
  if(type == "pigeons") {
    stopifnot(!is.null(n_chains) && !is.null(n_rounds))
    results_list = pigeons(jags_code, data, extract_variable_names, n_chains, n_rounds, extended_traces)
  } else if (type == "sat") {
    stopifnot(!is.null(N))
    results_list = istp.sat(jags_code,data, N)
  } else {
    stopifnot(!is.null(N))
    results_list = istp.cross(jags_code, data,  N)
  }
  return(results_list)
}


#' Drawing posterior samples by
#' using Importance sampling with Transparent reparameterization (ISTP)
#'
#' Return ESS and posterior samples found by ISTP algorithm.
#'
#'
#' @param posterior_fun function. Generating samples from
#' convenience posteriors posterior_fun(N)
#' @param h_inv function. An inverse function of h(theta) = (phi, lambda).
#' This function takes theta which is vector of paarameters.
#' @param N integer. The posterior sample size N
#' @param weight_fun function. Calculate weights by taking two arguments weight_fun(theta_matrix, deter)
#' @param var_names vector of string. variable names.
#' @return A list of results including samples and ESS.
#' @examples
#' data(small_cross)
#' data = small_cross
#' posterior_fun<-function(N){
#'   eta<- rdirichlet(N,data+1)
#'   se <- rbeta(N,25,3)
#'   sp <- rbeta(N,30,1.5)
#'   phi_lambda_matrix = cbind(eta[,1:3],se,sp)
#'   colnames(phi_lambda_matrix) = c(paste0("eta[",1:3,"]"),"se","sp")
#'   return(phi_lambda_matrix)
#' }
#' weight_fun<-function(theta_matrix, deter){
#'     apply(theta_matrix, 1,
#'       function(theta) ifelse(valid_A(theta),(theta[4]+theta[5]-1)^(-2),0))
#' }
#' valid_A <-function(theta){
#' pi = theta[1:3]
#' pi = c(pi,1-sum(pi))
#' all(pi>=0 & pi <=1)
#' }
#' h_cross_inv <- function(phiLambda) {
#' eta11 = ifelse(length(names(phiLambda)) == 5, phiLambda["eta[1]"],phiLambda[1])
#' eta12 = ifelse(length(names(phiLambda)) == 5, phiLambda["eta[2]"],phiLambda[2])
#' eta21 = ifelse(length(names(phiLambda)) == 5, phiLambda["eta[3]"],phiLambda[3])
#' se = ifelse(length(names(phiLambda)) == 5, phiLambda["se"],phiLambda[4])
#' sp = ifelse(length(names(phiLambda)) == 5, phiLambda["sp"],phiLambda[5])
#' eta22 = 1 - eta11 - eta12 - eta21
#' denom = se+sp-1
#' pi11 = 1/denom*(sp*eta11 + (sp-1)*eta21)
#' pi12 = 1/denom*(sp*eta12 + (sp-1)*eta22)
#' pi21 = 1/denom*((se-1)*eta11 + se*eta21)
#' return(c(pi11,pi12,pi21,se,sp))
#' }
#' var_name = c(paste0("pi","[",c(11,12,21),"]"),"se","sp")
#' results = istp.manual(posterior_fun, h_cross_inv, 15000, weight_fun, var_name)
#' head(results$samples)
#' results$ess
#' @export
istp.manual<-function(posterior_fun, h_inv, N, weight_fun, var_names){
  if(N <= 0) stop("Posterior sample size N > 0")
  phi_lambda_matrix = posterior_fun(N)
  deter = sapply(1:N,
              function(i) det(pracma::jacobian(h_inv, phi_lambda_matrix[i,])))
  theta_matrix = t(apply(phi_lambda_matrix,1,h_inv))
  # Weight
  wht <- weight_fun(theta_matrix, deter)
  index.na = is.na(wht)
  prop.na = sum(index.na)/N
  wht[index.na] = 0
  wht <- wht / sum(wht)
  resample_index<-sample(
    x = 1:N,
    replace = T,
    prob = wht,
    size = N)
  theta_matrix<-theta_matrix[resample_index,]
  colnames(theta_matrix)<-var_names
  list(
    prop.na = prop.na,
    samples = theta_matrix,
    ess = 1 / sum(wht ^ 2))
}


# Internal function
istp.sat<-function(jags_code, data, N){
  prep = prepare_sat(jags_code,data)
  results = istp.manual(posterior_fun = prep$posterior_fun, prep$h_inv, N, prep$weight_fun, prep$var_name)
  return(results)
}

# Internal function
istp.cross<-function(jags_code, data, N){
  prep = prepare_cross(jags_code,data)
  results = istp.manual(posterior_fun = prep$posterior_fun, h_cross_inv, N, prep$weight_fun, prep$var_name)
  return(results)
}

# Internal function
pigeons<-function(jags_code, data, extract_variable_names, n_chains, n_rounds, extended_traces){
  julia_library("Pigeons")
  parse_jags_code(jags_code,data, FALSE)
  julia_command("target= JuliaBUGSPath(model)")
  julia_assign("n_chains",n_chains)
  julia_assign("n_rounds",n_rounds)
  julia_assign("et", extended_traces)
  multithreaded = FALSE
  if(julia_eval("Threads.nthreads()") > 1){
    multithreaded = TRUE
  }
  julia_command("struct varExtractor end")
  parameters = julia_eval("string.(model.parameters)")
  if(is.null(extract_variable_names)){
   julia_command("Pigeons.extract_sample(state, log_potential, extractor::varExtractor) =
    Pigeons.extract_sample(state, log_potential)[1:end-1]")
   julia_command("Pigeons.sample_names(state, log_potential, extractor::varExtractor) =
    Pigeons.sample_names(state, log_potential)[1:end-1]")
  } else{
    stopifnot(all(extract_variable_names %in% parameters))
    index = match(extract_variable_names, parameters)
    julia_assign("index", index)
    julia_command("index = Int64.(index)")
    julia_assign("names", extract_variable_names)
    julia_command("Pigeons.extract_sample(state, log_potential, extractor::varExtractor) =
    Pigeons.extract_sample(state, log_potential)[index]")
    julia_command("Pigeons.sample_names(state, log_potential, extractor::varExtractor) =
    names")
  }
  julia_assign("multithreaded", multithreaded)
  julia_command("pt = pigeons(
        target = target,
        record = [traces; record_default()],
        extended_traces = et,
        extractor = varExtractor(),
        n_chains = Int64(n_chains),
        n_rounds = Int64(n_rounds),
        multithreaded = multithreaded
    )")
  pt = julia_eval("pt")
  return(pt)
}
