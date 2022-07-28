rm(list=ls())
set.seed(1)
library(MASS)
library(doParallel)
library(foreach)
library(mcmcse)

source("gradnt_log.R")

source("./../ebs_batch_mean.R")
source("./../ibs_jasa_mean.R")
source("./../sqrt_mat.R")

max_sam <- 1e5
nparm <- 5

#log_batch_fn(max_sam = max_sam, nparm = nparm, Rep = 5, ncores_par = 5, nam_matrix = "indep", cns = c(0.5, 1, 1.5, 2))# Low dimesnional set-up
nam_matrix = "indep"




foo <- paste("out/logistic_", nam_matrix,"_n_",max_sam,"_dim_",nparm,".RData",sep="")
load(foo)


#For dimension 5
cns = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,  1, 1.5, 2)
colMeans(cover_ibs)
colMeans(cover_orc)
cov_ebs_matr <- array(dim = c(3,length(cns),7), dimnames = list(  c("beta1", "beta2", "beta3"), cns, 1:7))
for (k in 1: 7)
{ for(l in 1:length(cns)){
  tmp1 <- (3*l-2)
  tmp2 <- 3*l
  cov_ebs_matr[,l, k  ] <-  colMeans(cover_ebs[k,, tmp1 : tmp2])
  
}
}
cov_ebs_ls_matr <- array(dim = c(3,length(cns),7), dimnames = list(  c("beta1", "beta2", "beta3"), cns, 1:7))
for (k in 1: 7)
{ for(l in 1:length(cns)){
  tmp1 <- (3*l-2)
  tmp2 <- 3*l
  cov_ebs_ls_matr[,l, k  ] <-  colMeans(cover_ebs_ls[k,, tmp1 : tmp2])
  
}
}

