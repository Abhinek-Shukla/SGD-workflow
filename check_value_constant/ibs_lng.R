ibs_lng <- function(n,cns,alp){
  am <- numeric(1000)
  # JASA Online BM Estimators
  am <- ceiling( cns*(1:1e3)^(2/(1-alp)))
  am <- c(am[am < n], (n+1))
  am <- unique(am)
  length(am)
}