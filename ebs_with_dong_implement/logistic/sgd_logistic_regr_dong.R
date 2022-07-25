#rm(list=ls())
set.seed(1)
library(MASS)
library(doParallel)
library(foreach)
library(mcmcse)

source("gradnt_log_dong_implement.R")
source("./../ebs_batch_mean_dong.R")
source("./../sqrt_mat.R")

max_sam <- 5e5
nparm <- 20
if(nparm < 10){ an = 20}
if(nparm > 10){an = nparm + 8}

#an = 40
nam_matrix = "indep"

log_batch_fn(max_sam = max_sam, an = an, eta_cns = 2,   nparm = nparm, Rep = 2, ncores_par = 1, nam_matrix = nam_matrix)

foo <- paste("out/logistic_", nam_matrix,"_n_",max_sam,"_dim_",nparm,"_dong.RData",sep="")
load(foo)

colMeans(cover_orc)
colMeans(cover_ebs)
colMeans(cover_ebs_ls)

