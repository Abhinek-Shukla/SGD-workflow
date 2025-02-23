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

dta <- read.csv(file="diabetes_binary.csv")

train_percent <- 50
  
train_indx <- as.integer(0.01*train_percent*nrow(dta))
train      <- dta[1:train_indx, ];
test       <- dta[((train_indx+1):nrow(dta)), ];

train[, 1] <- 2 * train[, 1] - 1 

cutoffs <- seq(0.2, 0.9, length = 100)  
log_batch_fn(dta = train, test = test, eta_cns = 1, 
             init = 10000, cutoffs = cutoffs)

load("./out/logistic_pred.Rdata")
plot_misclass(cutoffs, misclass, misclass.lb, titl = "diabetes")

