rm(list=ls())
set.seed(1)
library(MASS)
library(doParallel)
library(foreach)
library(mcmcse)
library(mvtnorm)
source("logistic_main.R")
source("./../ebs.R")
source("./../ibs.R")
source("./../misc.R")

# First download train.csv and test.csv from the following link 
# https://www.kaggle.com/competitions/santander-customer-transaction-prediction/data

dta_org <- read.csv(file="training.csv", sep=",")
dta            <- dta_org[, -1]
# dta            <- matrix(as.numeric(dta), nrow = nrow(dta), ncol = ncol(dta))
# Transform into -1,1 type of binary setup, Chen AOS (2020).
dta[, 1]       <- 2 * dta[, 1] - 1 


log_batch_fn(max_sam = 2e5, eta_cns = 0.05, alp = .51, cns = c(0.1), 
             cns1 = 0.01, burn_in = 5000, sam_siz = c(1e5,2e5), dta = dta)

load("out/logistic_real_dim_50.RData")


# Joint region Volume Comparison
c( volm_ibs[2], volm_ebs[2, 2], volm_ebs_ls[2, 2])

# Marginal friendly inferences 

# Max Ratio of lengths of intervals among different dimensions 
ibs_ebs_comprs_max <-  max(ratio_ibs_ebs[2, 2, ])
ibs_ebs_ls_comprs_max <-  max(ratio_ibs_ebs_ls[2, 2, ])
ebs_ls_ebs_comprs_max <-  max(ratio_ebs_ls_ebs[2, 2, ])
c(ibs_ebs_comprs_max, ibs_ebs_ls_comprs_max, ebs_ls_ebs_comprs_max)


# min Ratio of lengths of intervals among different dimensions 
ibs_ebs_comprs_min <-  min(ratio_ibs_ebs[2, 2, ])
ibs_ebs_ls_comprs_min <-  min(ratio_ibs_ebs_ls[2, 2, ])
ebs_ls_ebs_comprs_min <-  min(ratio_ebs_ls_ebs[2, 2, ])
c(ibs_ebs_comprs_min, ibs_ebs_ls_comprs_min, ebs_ls_ebs_comprs_min)


# Mean Ratio of lengths of intervals among different dimensions 
ibs_ebs_comprs_mean <-  mean(ratio_ibs_ebs[2, 2, ])
ibs_ebs_ls_comprs_mean <-  mean(ratio_ibs_ebs_ls[2, 2, ])
ebs_ls_ebs_comprs_mean <-  mean(ratio_ebs_ls_ebs[2, 2, ])
c(ibs_ebs_comprs_mean, ibs_ebs_ls_comprs_mean, ebs_ls_ebs_comprs_mean)

# Marginal friendly region volume of cuboids  inflated with respect to joint 
c(marg_volm_ibs[2], marg_volm_ebs[2, 2], marg_volm_ebs_ls[2, 2])




tmp1 <- margn_up_low[[1]]
len_1 <- (tmp1[, 2] - tmp1[, 1])
indx <- c(1 : 5, 46 : 50)
tmp2 <- tmp1[order(len_1), ]

tmp2[indx, ]#IBS

tmp1 <- margn_up_low[[4]]
len_1 <- (tmp1[, 2] - tmp1[, 1])
indx <- c(1 : 5, 46 : 50)
tmp2 <- tmp1[order(len_1), ]

tmp2[indx, ]#EBS


tmp1 <- margn_up_low[[5]]
len_1 <- (tmp1[, 2] - tmp1[, 1])
indx <- c(1 : 5, 46 : 50)
tmp2 <- tmp1[order(len_1), ]

tmp2[indx, ] #EBS + Lugsail
