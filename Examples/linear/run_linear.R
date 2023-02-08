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


# p = 5 runs   
max_sam <- 5e6
nparm <- 5
ncores <- 2
Reps <- 2


linear_reps(max_sam = max_sam, nparm = nparm, Rep = Reps, ncores_par = ncores,
            nam_matrix = "indep", cns = c(0.1))

# p = 20 runs   
max_sam <- 5e6
nparm <- 20
ncores <- 2
Reps <- 2

linear_reps(max_sam = max_sam, nparm = nparm, Rep = Reps, ncores_par = ncores, 
            nam_matrix = "indep", cns = c(0.1), cns1 = 0.1, eta_cns = 1)
