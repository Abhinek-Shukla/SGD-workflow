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

max_sam <- 1e6
nparm <- 20

nam_matrix = "indep"


log_batch_fn(max_sam = max_sam, eta_typ = 2, eta_cns = 1.5,   nparm = nparm, Rep = 1000, ncores_par = 25, nam_matrix = "indep", cns = seq(0.2, 0.6, by =0.1))# High dimesnional set-up


log_batch_fn(max_sam = max_sam, eta_typ = 3, eta_cns = 2,   nparm = nparm, Rep = 1000, ncores_par = 25, nam_matrix = "indep", cns = seq(0.2, 0.6, by =0.1))# High dimesnional set-up


log_batch_fn(max_sam = max_sam, eta_typ = 4, eta_cns = 3,   nparm = nparm, Rep = 1000, ncores_par = 25, nam_matrix = "indep", cns = seq(0.2, 0.6, by =0.1))# High dimesnional set-up
