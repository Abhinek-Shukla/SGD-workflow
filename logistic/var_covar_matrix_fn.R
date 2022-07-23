
var_covar_matrix_fn <- function(nparm, A = diag(nparm), niter = 5e6){
  library(MASS)
  set.seed(1)
parm <- seq(1 / nparm, 1, length.out = nparm)

sigm_inv <- matrix(rep(0,nparm^2),nrow = nparm, ncol = nparm)
x1 <- mvrnorm(niter, mu = rep(0, nparm), Sigma = A)

for ( rp in 1:niter){
 
   tmp2 <- (c((1 + exp(t(x1[rp,]) %*% parm))*(1 + exp(-t(x1[rp,]) %*% parm))*niter))
  sigm_inv <- sigm_inv + x1[rp,] %*% t(x1[rp,]) / tmp2

  }

return(sigm_inv)

}

sigm_inv_indep_nparm_5 <- var_covar_matrix_fn(nparm = 5)


r <- 0.5
toep_mat <- equiv_mat <- matrix(nrow = 20, ncol =20)

for( i in 1 : 20){
  
  for(j in 1 : 20){
    
    toep_mat[i,j] <- r^(abs(i-j))
    
    
    if(i == j ){
      
      equiv_mat[i, j] <- 1
      
    }else{ 
      
      equiv_mat[i, j] <- r 
      
    }
  }
  
}

sigm_inv_equiv_nparm_5 <- var_covar_matrix_fn(nparm = 5, A = equiv_mat[1 : 5, 1 : 5])
sigm_inv_toep_nparm_5 <- var_covar_matrix_fn(nparm = 5, A = toep_mat[1 : 5, 1 : 5])




sigm_inv_indep_nparm_20 <- var_covar_matrix_fn(nparm = 20)
sigm_inv_equiv_nparm_20 <- var_covar_matrix_fn(nparm = 20, A = equiv_mat)
sigm_inv_toep_nparm_20 <- var_covar_matrix_fn(nparm = 20, A = toep_mat)
save(sigm_inv_indep_nparm_5, sigm_inv_equiv_nparm_5, sigm_inv_toep_nparm_5, sigm_inv_indep_nparm_20, sigm_inv_equiv_nparm_20,  sigm_inv_toep_nparm_20,  file = "out/sigm_inv_logist.RData" )
