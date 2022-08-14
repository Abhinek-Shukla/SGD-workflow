rm(list=ls())
set.seed(1)
library(MASS)
library(doParallel)
library(foreach)
library(mcmcse)
library(smoothmest)

source("grad_lad_and_batch.R")
source("./../ebs_batch_mean.R")
source("./../ibs_jasa_mean.R")
source("./../sqrt_mat.R")

max_sam <- 1e6
nparm <- 20

lad_batch_fn(max_sam = max_sam, nparm = nparm, Rep = 20, ncores_par = 20, cns = c(0.01, 0.05, 0.1, 0.2, 0.5))#max(detectCores() - 5, 1)

cns = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)# cns = 1
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





