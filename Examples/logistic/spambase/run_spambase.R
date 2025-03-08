set.seed(1)
library(MASS)
library(doParallel)
library(foreach)
library(mcmcse)
library(mvtnorm)
library(tidyr)
source("./../logistic_back_funcs.R")
source("./../pred_main.R")
source("./../plot_misclass.R")
source("./.././../ebs.R")

dta_org    <- as.data.frame(read.table(file="spambase.data"),
                         stringsAsFactors = FALSE)
split_list <- strsplit(dta_org$V1, ",")


max_cols     <- max(sapply(split_list, length))
split_matrix <- do.call(rbind, lapply(split_list, 
                       function(x) c(x, rep(NA, max_cols - length(x)))))

df_new <- as.data.frame(split_matrix, stringsAsFactors = FALSE)

head(df_new)

train_percent <- 50

# create dataset
# 58 response
# 1:57 continuous covariates
dta                <- matrix(0, nrow = nrow(df_new), ncol = ncol(df_new))
dta[, 2:ncol(dta)] <- as.numeric(as.matrix(df_new[, 1:(ncol(df_new)-1)]))
dta[, 2:ncol(dta)] <- scale(dta[, 2:ncol(dta)])

# randomly shuffle all rows of dta so as to avoid consecutive 1's or 0's in response #
dta <- dta[sample(nrow(dta)), ]
train_indx <- as.integer(0.01*train_percent*nrow(dta))
train      <- dta[1:train_indx, ];
test       <- dta[((train_indx+1):nrow(dta)), ];

train[, 1] <- 2 * train[, 1] - 1 

cutoffs <- seq(0.2, 0.9, length = 100)  
log_batch_fn(dta = train, test = test, eta_cns = 4.5, burn_in = 500, 
             init = 500, cutoffs = cutoffs)

load("./out/logistic_pred.Rdata")
plot_misclass(cutoffs, misclass, misclass.lb, titl = "spambase")

