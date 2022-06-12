sqrt_mat <- function(sigm){
decomp_sigm<- svd(sigm)
dig_sig <- decomp_sigm$d
v_mat <- decomp_sigm$u
sqrt_sig <- v_mat%*%diag(sqrt(dig_sig))%*%t(v_mat)

return(sqrt_sig)

}