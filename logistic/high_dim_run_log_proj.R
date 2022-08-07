rm(list=ls())
set.seed(1)
library(MASS)
library(doParallel)
library(foreach)
library(mcmcse)

source("gradnt_log_proj.R")

source("./../ebs_batch_mean.R")
source("./../ibs_jasa_mean.R")
source("./../sqrt_mat.R")

max_sam <- 1e6
nparm <- 20

nam_matrix = "indep"


log_batch_fn(max_sam = max_sam,    nparm = nparm, Rep = 1, ncores_par = 1, nam_matrix = "indep", cns = seq(0.5, 4, by =0.5))# High dimesnional set-up


#For dimension 20
cns = seq(0.5, 4, by =0.5)# cns = 1
colMeans(cover_ibs)
colMeans(cover_orc)
cov_ebs_matr <- array(dim = c(3,length(cns),5), dimnames = list(  c("beta1", "beta2", "beta3"), cns, 1:5))
for (k in 1: 5)
{ for(l in 1:length(cns)){
  tmp1 <- (3*l-2)
  tmp2 <- 3*l
  cov_ebs_matr[,l, k  ] <-  colMeans(cover_ebs[k,, tmp1 : tmp2])
  
}
}
cov_ebs_ls_matr <- array(dim = c(3,length(cns),5), dimnames = list(  c("beta1", "beta2", "beta3"), cns, 1:5))
for (k in 1: 5)
{ for(l in 1:length(cns)){
  tmp1 <- (3*l-2)
  tmp2 <- 3*l
  cov_ebs_ls_matr[,l, k  ] <-  colMeans(cover_ebs_ls[k,, tmp1 : tmp2])
  
}
}
