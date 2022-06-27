ibs_lng <- function(dm,cns,alp){
  am <- numeric(1000)
  for(m in 1:10000){
    
    am[m] <- floor(cns*m^(2/(1-alp)))
    
  }
  am <- unique(am)
  am[dm]+1
}