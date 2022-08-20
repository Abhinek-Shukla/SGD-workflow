set.seed(1)
library(MASS)
library(doParallel)
library(foreach)
library(mcmcse)
library(mvtnorm)
source("linear_main.R")
source("./../ebs.R")
source("./../ibs.R")
source("./../misc.R")


####################################################################################################################
##     p = 5 runs   
####################################################################################################################
max_sam <- 1e6
nparm <- 5
ncores <- 2
Reps <- 2

#Toeplitz and Equivariance-covariance matrices definitions
rho <- 0.5
foo <- X_mat(nparm = nparm, rho = rho)
equiv_mat <- foo[[1]]
toep_mat <- foo[[2]]

#Identity variance-covaiance matrix case
linear_reps(max_sam = max_sam, nparm = nparm, Rep = Reps, ncores_par = ncores, nam_matrix = "indep", cns = c(0.1, 1))


#Toeplitz variance-covaiance matrix results
linear_reps(max_sam = max_sam, A = toep_mat, nparm = nparm, Rep = Reps, ncores_par = ncores, nam_matrix = "toep")

#Equivariance-covariance matrix Results
linear_reps(max_sam = max_sam, A = equiv_mat, nparm = nparm, Rep = Reps, ncores_par = ncores, nam_matrix = "equiv")

####################################################################################################################


####################################################################################################################
##     p = 20 runs   
####################################################################################################################
max_sam <- 1e6
nparm <- 20
ncores <- 2
Reps <- 2

#Toeplitz and Equivariance-covariance matrices definitions
rho <- 0.5
foo <- X_mat(nparm = nparm, rho = rho)
equiv_mat <- foo[[1]]
toep_mat <- foo[[2]]

# Identity
linear_reps(max_sam = max_sam, nparm = nparm, Rep = Reps, ncores_par = ncores, nam_matrix = "indep", cns = c(0.1), cns1 = 0.1, eta_cns = 1)

# Toeplitz 
linear_batch_fn(max_sam = max_sam, A = toep_mat, nparm = nparm, Rep = Reps, ncores_par = ncores, nam_matrix = "toep", cns = c(0.1), cns1 = 0.1, eta_cns = 1)

# Equivariance
linear_batch_fn(max_sam = max_sam, A = equiv_mat, nparm = nparm, Rep = Reps, ncores_par = ncores, nam_matrix = "equiv", cns = c(0.1), cns1 = 0.1, eta_cns = 1)


