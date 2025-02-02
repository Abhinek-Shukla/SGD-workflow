set.seed(10)
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

# choosing rich enough covariates to avoid imbalance
# at least .01 responses should be available
index <- which(colSums(dta_org[ , 11:54]) > dim(dta_org)[1]/100) + 10

train_percent <- 50
  
# scale the continuous covariates
temp <- dta_org[ , 1:10]
temp <- scale(temp)

# create dataset
# 55 response
# 1:10 continuous covariates
# index rich enough covariates
dta      <- dta_org[, c(55, 1:10, index)]
dta[, 1] <- ifelse(dta[, 1] >= 2, 1, 0)
# replace continuous covariates with centered and scaled versions
dta[, 2:11] <- temp
  
train_indx <- as.integer(0.01*train_percent*nrow(dta))
train      <- dta[1:train_indx, ];
test       <- dta[((train_indx+1):nrow(dta)), ];

train[, 1] <- 2 * train[, 1] - 1 

cutoffs <- seq(0.2, 0.9, length = 100)  
log_batch_fn(dta = train, test = test, eta_cns = 100, cutoffs = cutoffs)

load("./out/logistic_pred.Rdata")
plot_misclass(cutoffs, misclass, misclass.lb)
