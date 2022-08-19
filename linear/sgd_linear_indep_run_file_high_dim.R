rm(list=ls())
set.seed(1)
library(MASS)
library(doParallel)
library(foreach)
library(mcmcse)
library(mvtnorm)
source("grad_lin_and_batch.R")
source("./../ebs_batch_mean.R")
source("./../ibs_jasa_mean.R")
source("./../sqrt_mat.R")

max_sam <- 1e6
nparm <- 20

linear_batch_fn(max_sam = max_sam, nparm = nparm, Rep = 2, ncores_par = 1, nam_matrix = "indep", cns = c( 0.1), cns1 = 0.1, eta_cns = 1 )#max(detectCores() - 5, 1)



#Toeplitz and Equivariance-covariance matrices definitions

r <- 0.5
toep_mat <- equiv_mat <- matrix(nrow = nparm, ncol = nparm)

for( i in 1 : nparm)
{
  
  for(j in 1 : nparm)
  {
    
    toep_mat[i,j] <- r^(abs(i-j))
    
    
    if(i == j )
    {
      
      equiv_mat[i, j] <- 1
      
    }else
    { 
      
      equiv_mat[i, j] <- r 
      
    }
  }
  
}



#Toeplitz variance-covaiance matrix results
linear_batch_fn(max_sam = max_sam, nparm = nparm, Rep = 2, ncores_par = 1, nam_matrix = "toep", cns = c( 0.1), cns1 = 0.1, eta_cns = 1 )


#Equivariance-covariance matrix Results
linear_batch_fn(max_sam = max_sam, nparm = nparm, Rep = 2, ncores_par = 1, nam_matrix = "equiv", cns = c( 0.1), cns1 = 0.1, eta_cns = 1 )







