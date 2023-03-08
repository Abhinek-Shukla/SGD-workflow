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
