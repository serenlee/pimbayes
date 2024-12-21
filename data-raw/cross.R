# Example 2 data
data_cross_sectional<-function(n){
  p = 0.3
  q = 0.4
  e = 0.5
  se = 0.95
  sp = 0.95
  pi = c(p*e,(1-p)*e,q*(1-e),(1-q)*(1-e))
  eta = c(se*pi[1]+(1-sp)*pi[3],se*pi[2]+(1-sp)*pi[4],sp*pi[3]+(1-se)*pi[1],(1-se)*pi[2]+sp*pi[4])
  set.seed(777)
  as.vector(rmultinom(1,n,eta))
}

small_cross<-data_cross_sectional(38)
cross<-data_cross_sectional(38000)


usethis::use_data(small_cross, overwrite = TRUE)
usethis::use_data(cross, overwrite = TRUE)
