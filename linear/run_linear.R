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
toep_mat <- foo[[2]]


linear_reps(max_sam = max_sam, nparm = nparm, Rep = Reps, ncores_par = ncores, nam_matrix = "indep", cns = c(0.1, 1))
linear_reps(max_sam = max_sam, A = toep_mat, nparm = nparm, Rep = Reps, ncores_par = ncores, nam_matrix = "toep")

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
toep_mat <- foo[[2]]

linear_reps(max_sam = max_sam, nparm = nparm, Rep = Reps, ncores_par = ncores, nam_matrix = "indep", cns = c(0.1), cns1 = 0.1, eta_cns = 1)
linear_reps(max_sam = max_sam, A = toep_mat, nparm = nparm, Rep = Reps, ncores_par = ncores, nam_matrix = "toep", cns = c(0.1), cns1 = 0.1, eta_cns = 1)


