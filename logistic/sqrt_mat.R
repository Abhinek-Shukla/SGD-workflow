sqrt_mat <- function(sigm){
  p <- nrow(sigm)
decomp_sigm<- svd(sigm)
dig_sig <- decomp_sigm$d
v_mat <- decomp_sigm$u
sqrt_sig <- v_mat%*%diag((dig_sig)^(1/2))%*%t(v_mat)

return(sqrt_sig)

}