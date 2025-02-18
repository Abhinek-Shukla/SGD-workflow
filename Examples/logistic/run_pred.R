set.seed(1)
library(MASS)
library(doParallel)
library(foreach)
library(mcmcse)
library(mvtnorm)
source("logistic_back_funcs.R")
source("pred_main.R")
source("./../ebs.R")


# First download train.csv and test.csv from the following link 
# https://www.kaggle.com/competitions/santander-customer-transaction-prediction/data

dta_org <- read.csv(file="training.csv")
dta            <- dta_org[, -1]
# dta            <- matrix(as.numeric(dta), nrow = nrow(dta), ncol = ncol(dta))
# Transform into -1,1 type of binary setup, Chen AOS (2020).
dta[, 1]       <- 2 * dta[, 1] - 1 

test <- read.csv("testing.csv", header = TRUE)
test <- test[,-1]

log_batch_fn(dta = dta, test = test)


