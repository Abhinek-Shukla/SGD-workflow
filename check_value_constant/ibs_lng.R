ibs_lng <- function(n,cns,alp){
  am <- numeric(1000)
  for(m in 1:10000){
    
    am[m] <- floor(cns*m^(2/(1-alp)))
    
  }
  am <- am[am<n]
  am <- unique(am)
  length(am)
}