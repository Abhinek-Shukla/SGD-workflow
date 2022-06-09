grad_lin <- function(sg,y,x){
  (x %*% (sg) - y) %*% x
} 
