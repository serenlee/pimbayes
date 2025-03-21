# Necessary tools for the package

patch_jags_code = function(jags_code) {
  jags_code = gsub("ddirch", "ddirich", jags_code)
  jags_code = gsub("<-", "=", jags_code)
  jags_code
}

check_integer <- function(x){
  if(any(is.na(x))) return(FALSE)
  if(is.integer(x)){
    return(TRUE)
  } else {
    return(all(as.integer(x) == x))
  }
}

parse_jags_code <- function(jags_code, data, is_data_info_required = T) {
  julia_library("JuliaBUGS")
  julia_library("LinearAlgebra")
  julia_command("JuliaBUGS.@register_primitive function Int64wrapper(vec::Vector)
    return Int64.(vec)
end")
  julia_command("JuliaBUGS.@register_primitive function c(a...)
    return [a...]
end")
  julia_assign("data_dict", data)
  integer_list = names(data)[which(sapply(data,check_integer))]
  julia_command("data = NamedTuple{Tuple(keys(data_dict))}(values(data_dict))")
  julia_assign("int_keys",integer_list)
  julia_command("vals = Any[]")
  julia_command("for key in int_keys
                value = data_dict[Symbol(key)]
                push!(vals, Int64.(value))
                end
                ")
  julia_command("data = merge(data, (; zip(Symbol.(int_keys), vals)...))")
  julia_command(paste0('model_def = @bugs("""',patch_jags_code(jags_code),'""",false,false)'))
  julia_command("model = compile(model_def, data)")
  if(!is_data_info_required) return()
  julia_command("parameters = model.parameters")
  parameters = julia_eval("string.(parameters)")
  julia_command("node_data = JuliaBUGS.FlattenedGraphNodeData(model.g)")
  julia_command("sorted_nodes = node_data.sorted_nodes")
  julia_command("evaluation_env = model.evaluation_env")
  julia_command("info_vec = String[]")
  julia_command("is_stochastic_vec = Bool[]")
  julia_command("is_observed_vec = Bool[]")
  julia_command("for (i,vn) in enumerate(sorted_nodes)
                  node_function = node_data.node_function_vals[i];
                  loop_vars = node_data.loop_vars_vals[i];
                  val = Base.invokelatest(node_function, evaluation_env, loop_vars);
                  push!(info_vec, string(val));
                  push!(is_stochastic_vec, node_data.is_stochastic_vals[i]);
                  push!(is_observed_vec, node_data.is_observed_vals[i]);
                end")
  varnames = julia_eval("string.(sorted_nodes)")
  varinfo = julia_eval("info_vec")
  is_stochastic = julia_eval("is_stochastic_vec")
  is_observed = julia_eval("is_observed_vec")
  var_data = data.frame(varnames, varinfo, is_stochastic,is_observed)
  return(list(parameters, var_data))
}

find_index <- function(frame, keyword){
  unlist(sapply(keyword,
         function(key) which(sapply(frame[,2], stringr::str_detect, pattern = key))))
}

extract_name<-function(string){
  unique(sapply(strsplit(string,"\\["), function(x) x[1]))
}

prepare_sat<-function(jags_code,data){
  parsed = parse_jags_code(jags_code,data)
  var_data = parsed[[2]]
  index_of_log_or_r = find_index(var_data,"Normal")
  log_or_r_name = extract_name(var_data[,1][index_of_log_or_r[1]])
  pr_xy1_name = extract_name(var_data[!var_data[,3] & !var_data[,4],1])[1]
  pr_xy0_candidacy_ind = find_index(var_data,c("Uniform","Beta"))
  pr_xy0_candidacy = extract_name(var_data[,1][pr_xy0_candidacy_ind])
  diff = julia_eval(paste0("evaluation_env.",log_or_r_name, "+ map(JuliaBUGS.logit, evaluation_env.",pr_xy0_candidacy[1],") - map(JuliaBUGS.logit, evaluation_env.",pr_xy1_name,")"))
  pr_xy0_name = ifelse(all(abs(diff)<1e-5), pr_xy0_candidacy[1], pr_xy0_candidacy[2])
  py_x_name = ifelse(all(abs(diff)<1e-5), pr_xy0_candidacy[2], pr_xy0_candidacy[1])
  px_ind = find_index(var_data,"Dirichlet")
  px_name = extract_name(var_data[,1][px_ind])
  posterior_fun<-function(N){
    sat_posterior_fun(N, 1/sqrt(data$tau), log2(data$p2), data$r, data$x, data$y)
  }
  h_inv<-function(philambda){
    h.sat.inv(philambda,log2(data$p2))
  }
  p <-log2(data$p2)
  weight_fun<-function(theta_matrix, deter){
    px <- theta_matrix[,1:(2^p-1)]
    px <- cbind(1-rowSums.n(px),px)
    wht <- ddirichlet(px, rep(1,2^p)) *  #px
      rowProd(dunif(theta_matrix[,2^p:(2^p+2^p-1)], 0, 1)) * #qy
      rowProd(dunif(theta_matrix[,(2^p+2^p):(2^p+2^p+2^p-1)], 0, 1)) * abs(deter) #ry0, ry1 is canceled.
  }
  var_name = c(paste0(px_name,"[",2:(data$p2),"]"),
               paste0(py_x_name,"[",1:(data$p2),"]"),
               paste0(pr_xy0_name,"[",1:(data$p2),"]"),
               paste0(log_or_r_name,"[",1:(data$p2),"]"))
  return(list(posterior_fun = posterior_fun, h_inv= h_inv, weight_fun = weight_fun, var_name = var_name))
}




sat_posterior_fun<-function(N, sigma, p, r, x, y){
  pr <- rbeta(N, 1 + sum(r == 1), 1 + sum(r == 0))
  px.r0 <- rdirichlet(N, 1 + colSums.n(x[r == 0,]))[,-1]
  px.r1 <- rdirichlet(N, 1 + colSums.n(x[r == 1,]))[,-1]
  py.xr1 <- rdirichlet(N, 1 + colSums.n(x[r == 1 & y == 1,]))
  ## Convenience prior (lambda given phi)
  log.or.r <- matrix(rnorm(N*2^p, mean = 0, sd = sigma), nrow = N)
  phi_lambda_matrix = cbind(pr, px.r0, px.r1, py.xr1, log.or.r)
  phi_lambda_matrix
}

h.sat<-function(theta, p) {
  px = theta[1:(2^p-1)]
  px = c(1-sum(px), px)
  py.x = theta[2^p:(2^p+2^p-1)]
  pr.xy0 = theta[(2^p+2^p):(2^p+2^p+2^p-1)]
  log.or.r = theta[(2^p+2^p+2^p):(2^p+2^p+2^p+2^p-1)]
  pr.xy1 = expit(logit(pr.xy0) + log.or.r)
  pr = sum(py.x * px * pr.xy1 + (1-py.x) * px * pr.xy0)
  px.r0 = (py.x * px * (1-pr.xy1) + (1-py.x) * px * (1-pr.xy0))/(1-pr)
  px.r1 = (py.x * px * pr.xy1 + (1-py.x) * px * pr.xy0)/pr
  py.xr1 = (py.x * px * pr.xy1)/(px.r1*pr)
  return(c(pr,px.r0[-1],px.r1[-1],py.xr1,log.or.r))
}

h.sat.inv<-function(philambda, p){
  #if(any(is.na(philambda))) return(NA)
  pr = philambda[1]
  px.r0 = philambda[2:(2^p)]
  px.r1 = philambda[(2^p+1):(2^p+2^p-1)]
  py.xr1 = philambda[(2^p+2^p):(2^p+2^p+2^p-1)]
  log.or.r = philambda[(2^p+2^p+2^p):(2^p+2^p+2^p+2^p-1)]
  px.r0 = c((1-sum(px.r0)),px.r0)
  px.r1 = c((1-sum(px.r1)),px.r1)
  px = pr*px.r1 + (1-pr)*px.r0
  py.xr0 = expit(logit(py.xr1) - log.or.r)
  #if(any(pr.yxr0 > 1 | pr.yxr0 <0)) return(NA)
  py.x = (py.xr0 * (1-pr) * px.r0 + py.xr1 * pr * px.r1)/px
  pr.xy0 = (1 - py.xr1) * px.r1 * pr / ((1-py.x)* px)
  return(c(px[-1], py.x, pr.xy0, log.or.r))
}

rowProd<-function(mat){
  if(is.vector(mat)) return(mat)
  apply(mat,1,prod)
}


# Helper function for calculating ESS when importance weighting is applied
calcESS<-function(package.ess,importance.ess,n){
  round((package.ess)*(importance.ess)/n)
}

valid_A <-function(theta){
  pi = theta[1:3]
  pi = c(pi,1-sum(pi))
  all(pi>=0 & pi <=1)
}

# h function
h_cross <- function(theta){
  pi11 = ifelse(length(names(theta)) == 5, theta["pi[1]"],theta[1])
  pi12 = ifelse(length(names(theta)) == 5, theta["pi[2]"],theta[2])
  pi21 = ifelse(length(names(theta)) == 5, theta["pi[3]"],theta[3])
  se = ifelse(length(names(theta)) == 5, theta["se"],theta[4])
  sp = ifelse(length(names(theta)) == 5, theta["sp"],theta[5])
  pi22 = 1-pi11-pi12-pi21
  eta11 = se*pi11+(1-sp)*pi21
  eta12 = se*pi12+(1-sp)*pi22
  eta21 = sp*pi21+(1-se)*pi11
  eta22 = (1-se)*pi12+sp*pi22
  return(c(eta11,eta12,eta21,se,sp))
}

# h inverse function
h_cross_inv <- function(phiLambda) {
  eta11 = ifelse(length(names(phiLambda)) == 5, phiLambda["eta[1]"],phiLambda[1])
  eta12 = ifelse(length(names(phiLambda)) == 5, phiLambda["eta[2]"],phiLambda[2])
  eta21 = ifelse(length(names(phiLambda)) == 5, phiLambda["eta[3]"],phiLambda[3])
  se = ifelse(length(names(phiLambda)) == 5, phiLambda["se"],phiLambda[4])
  sp = ifelse(length(names(phiLambda)) == 5, phiLambda["sp"],phiLambda[5])
  eta22 = 1 - eta11 - eta12 - eta21
  denom = se+sp-1
  pi11 = 1/denom*(sp*eta11 + (sp-1)*eta21)
  pi12 = 1/denom*(sp*eta12 + (sp-1)*eta22)
  pi21 = 1/denom*((se-1)*eta11 + se*eta21)
  return(c(pi11,pi12,pi21,se,sp))
}

colSums.n <-function(matrix){
  if(is.vector(matrix)) return(matrix)
  return(colSums(matrix))
}

rowSums.n <- function(matrix){
  if(is.vector(matrix)) return(matrix)
  return(rowSums(matrix))
}

prepare_cross = function(jags_code,data){
  parsed = parse_jags_code(jags_code, data)
  var_data = parsed[[2]]
  pi_index = find_index(var_data, "Dirichlet")
  pi_name = extract_name(var_data[pi_index[1],1])
  se_index_candidacy = find_index(var_data, c("Beta","Uniform"))
  se_name_candidacy = var_data[se_index_candidacy,1]
  pi = julia_eval(paste0("evaluation_env.",pi_name))
  se_cand1 = julia_eval(paste0("evaluation_env.",se_name_candidacy[1]))
  sp_cand1 = julia_eval(paste0("evaluation_env.",se_name_candidacy[2]))
  eta_name = extract_name(var_data[!var_data$is_stochastic & !var_data$is_observed,1])
  eta = julia_eval(paste0("evaluation_env.",eta_name))
  se_name = ifelse(abs(se_cand1*pi[1]+(1-sp_cand1)*pi[3] - eta[1])<1e-5,
                   se_name_candidacy[1],
                   se_name_candidacy[2])
  sp_name = ifelse(abs(se_cand1*pi[1]+(1-sp_cand1)*pi[3] - eta[1])<1e-5,
                   se_name_candidacy[2],
                   se_name_candidacy[1])
  posterior_fun<-function(N){
    eta<- gtools::rdirichlet(N,data$y+1)
    se <- rbeta(N,25,3)
    sp <- rbeta(N,30,1.5)
    phi_lambda_matrix = cbind(eta[,1:3],se,sp)
    colnames(phi_lambda_matrix) = c(paste0("eta[",1:3,"]"),"se","sp")
    return(phi_lambda_matrix)
  }
  weight_fun<-function(theta_matrix, deter){
    apply(theta_matrix, 1, function(theta) ifelse(valid_A(theta),(theta[4]+theta[5]-1)^(-2),0))
  }
  var_name = c(paste0(pi_name,"[",c(11,12,21),"]"),se_name,sp_name)
  return(list(posterior_fun = posterior_fun, weight_fun = weight_fun, var_name = var_name))
}


