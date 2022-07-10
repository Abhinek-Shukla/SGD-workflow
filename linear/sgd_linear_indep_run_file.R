rm(list=ls())
set.seed(1)
library(MASS)
library(doParallel)
library(foreach)
library(mcmcse)
source("grad_lin_and_batch.R")
source("./../ebs_batch_mean.R")
source("./../ibs_jasa_mean.R")
source("./../sqrt_mat.R")

max_sam <- 1e5
nparm <- 5
linear_batch_fn(max_sam = max_sam, nparm = nparm, Rep = 20, ncores_par = 30)#max(detectCores() - 5, 1)

foo <- paste("out/linear_indep_n_",max_sam,"_dim_",nparm,".RData",sep="")
load(foo)


