set.seed(1)
library(MASS)
library(doParallel)
library(foreach)
library(mcmcse)
library(mvtnorm)
source("linear_back_funcs.R")
source("linear_main.R")
source("./../ebs.R")
source("./../ibs.R")
source("./../misc.R")


# p = 5 runs   
max_sam <-  5e6
nparm <- 5
ncores <- 50
Reps <- 1000

# covariance matrix of X calculations
rho <- 0.5
foo <- X_mat(nparm = nparm, rho = rho)
equiv_mat <- foo[[1]]
toep_mat <- foo[[2]]
st <- Sys.time()
linear_reps(max_sam = max_sam, nparm = nparm, Rep = Reps, ncores_par = ncores,
            nam_matrix = "indep", cns = 0.1)
print(Sys.time() - st)
linear_reps(max_sam = max_sam, A = toep_mat, nparm = nparm, Rep = Reps, ncores_par = ncores,
            nam_matrix = "toep", cns = 0.1)
linear_reps(max_sam = max_sam, A = equiv_mat, nparm = nparm, Rep = Reps, ncores_par = ncores,
            nam_matrix = "equiv", cns = 0.1)


# p = 20 runs   
max_sam <- 5e6
nparm <- 20
ncores <- 50
Reps <- 1000

linear_reps(max_sam = max_sam, nparm = nparm, Rep = Reps, ncores_par = ncores, 
            nam_matrix = "indep", cns = 0.1, cns1 = 0.1, eta_cns = 1)
linear_reps(max_sam = max_sam, nparm = nparm, Rep = Reps, ncores_par = ncores, 
            nam_matrix = "toep", cns = 0.1, cns1 = 0.1, eta_cns = 1)
linear_reps(max_sam = max_sam, nparm = nparm, Rep = Reps, ncores_par = ncores, 
            nam_matrix = "equiv", cns = 0.1, cns1 = 0.1, eta_cns = 1)