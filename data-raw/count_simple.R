data_count_simple<-function(n, seed = 7){
  set.seed(seed)
  y = rpois(n, lambda = 2)
  r <- rbinom(n, size=1,
              prob=expit(log(2)*((y-2)/sqrt(2))))
  y.obs = rep(NA,n)
  y.obs[r == 1] = y[r == 1]
  return(data.frame(y = y.obs, r = r))
}

small_count<-data_count_simple(100)
count<-data_count_simple(3000)


usethis::use_data(small_count, overwrite = TRUE)
usethis::use_data(count, overwrite = TRUE)
