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

max_sam <- 1e6
nparm <- 20

linear_batch_fn(max_sam = max_sam, nparm = nparm, Rep = 2, ncores_par = 1, nam_matrix = "indep", cns = c( 0.1), cns1 = 0.1, eta_cns = 1 )#max(detectCores() - 5, 1)



