set.seed(1)
library(MASS)
library(doParallel)
library(foreach)
library(mcmcse)
library(mvtnorm)
source("logistic_back_funcs.R")
source("pred_main.R")
source("plot_misclass.R")
source("ebs.R")

dta_org <- read.csv(file="covtype.csv")
cutoffs <- seq(0.01, 0.9, length = 500)

train_percent <- 70
  
  dta      <- dta_org[, c(55, c(1:20, 22:28, 30:53))]
  dta[, 1] <- ifelse(dta[, 1] >= 2, 1, 0)
  
  train_indx <- as.integer(0.01*train_percent*nrow(dta))
  train      <- dta[1:train_indx, ];
  test       <- dta[((train_indx+1):nrow(dta)), ];
  
  train[, 1] <- 2 * train[, 1] - 1 
  
log_batch_fn(dta = train, test = test, eta_cns = 0.0001, cutoffs = cutoffs)

load("./out/logistic_pred.Rdata")
plot_misclass(cutoffs, misclass, misclass.lb)