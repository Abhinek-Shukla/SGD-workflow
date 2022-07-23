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

log_batch_fn(max_sam = max_sam, nparm = nparm, Rep = 50, ncores_par = 1, nam_matrix = "indep", cns = c(0.5, 1, 1.5, 2))# Low dimesnional set-up
nam_matrix = "indep"
foo <- paste("out/logistic_", nam_matrix,"_n_",max_sam,"_dim_",nparm,".RData",sep="")
load(foo)
cns = c(0.5, 1, 1.5, 2)
colMeans(cover_ibs)
colMeans(cover_orc)

for (k in 1: (length(cns)*3))
{
  print( mean(cover_ebs[,,k]))
}
for (k in 1: (length(cns)*3))
{
  print( mean(cover_ebs_ls[,,k]))
}

