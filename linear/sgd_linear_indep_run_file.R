rm(list=ls())
library(MASS)

source("linear/grad_lin_and_batch.R")
source("ebs_batch_mean.R")
source("ibs_jasa_mean.R")
source("sqrt_mat.R")

linear_batch_fn(max_sam = 1e5, burn_in = 1000, nparm = 5, Rep = 1,  eta_cns = 0.5, qlev = 0.95, alp = .51, cns = c(0.1, 1) )
load("linear/linear_indep_rep_1_dim_5.RData")


