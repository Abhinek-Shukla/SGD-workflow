set.seed(1)
library(MASS)
library(doParallel)
library(foreach)
library(mcmcse)
library(mvtnorm)
source("./../logistic_back_funcs.R")
source("./../pred_main.R")
source("./../plot_misclass.R")
source("./.././../ebs.R")


# First download train.csv and test.csv from the following link #
# https://www.kaggle.com/competitions/santander-customer-transaction-prediction/data #

dta_org <- read.csv(file="training.csv")
dta     <- dta_org[, -1]

# Transform training response into -1,1 type of binary setup, Chen AOS (2020) #
dta[, 1] <- 2 * dta[, 1] - 1 

# Load test data #
test <- read.csv("testing.csv", header = TRUE)
test <- test[,-1]

cutoffs <- seq(0.01, 0.3, length = 100)

# Feed the training data and obtain misclassfication rates #
log_batch_fn(dta = dta, test = test, eta_cns = .05,
             burn_in = 5000, init = 10000, cutoffs = cutoffs)
load("./out/logistic_pred.Rdata")
plot_misclass(cutoffs, misclass, misclass.lb, "santander")


