gradnt_log <- function(y,a,sg){
  
  tmp <- (c( 1 + exp(  y*t(a )%*% (sg))))

  p_thet <- 1/ tmp
  
  return(-y*p_thet*a)
}