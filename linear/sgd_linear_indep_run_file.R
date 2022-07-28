rm(list=ls())
set.seed(1)
library(MASS)
library(doParallel)
library(foreach)
library(mcmcse)
source("./../opt_bet_fn.R")
source("grad_lin_and_batch.R")
source("./../ebs_batch_mean.R")
source("./../ibs_jasa_mean.R")
source("./../sqrt_mat.R")

max_sam <- 1e5
nparm <- 20

linear_batch_fn(max_sam = max_sam, nparm = nparm, Rep = 2, ncores_par = 2, nam_matrix = "indep", cns = c(0.001, 0.01))#max(detectCores() - 5, 1)

#foo <- paste("out/linear_", nam_matrix,"_n_",max_sam,"_dim_",nparm,".RData",sep="")
#load(foo)











r <- 0.5
toep_mat <- equiv_mat <- matrix(nrow = nparm, ncol = nparm)

for( i in 1 : nparm){
  
  for(j in 1 : nparm){
    
    toep_mat[i,j] <- r^(abs(i-j))
    
    
    if(i == j ){
      
      equiv_mat[i, j] <- 1
      
    }else{ 
      
      equiv_mat[i, j] <- r 
      
    }
  }
  
}




linear_batch_fn(max_sam = max_sam, A = toep_mat, nparm = nparm, Rep = 30, ncores_par = 15, nam_matrix = "toep")#max(detectCores() - 5, 1)

#foo <- paste("out/linear_", nam_matrix,"_n_",max_sam,"_dim_",nparm,".RData",sep="")
#load(foo)


linear_batch_fn(max_sam = max_sam, A = equiv_mat, nparm = nparm, Rep = 30, ncores_par = 15, nam_matrix = "equiv")#max(detectCores() - 5, 1)

#foo <- paste("out/linear_", nam_matrix,"_n_",max_sam,"_dim_",nparm,".RData",sep="")
#load(foo)





