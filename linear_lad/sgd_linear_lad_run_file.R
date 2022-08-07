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

lad_batch_fn(max_sam = max_sam, nparm = nparm, Rep = 20, ncores_par = 20, cns = c(0.01, 0.1, 0.5))#max(detectCores() - 5, 1)

#foo <- paste("out/linear_", nam_matrix,"_n_",max_sam,"_dim_",nparm,".RData",sep="")
#load(foo)








