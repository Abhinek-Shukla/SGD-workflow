rm(list=ls())
library(MASS)
library(matrixcalc)
source("grad_lin.R")
source("ebs_batch_mean.R")
source("ibs_jasa_mean.R")
source("sqrt_mat.R")
Rep <- 1
cutf <- 1000 #Dropping initial Iterates of SGD
#Sample Size
n <- 5e4+cutf;
#Confidence level 
qlev <- 0.95
#Iterations
Iter <- n;
alp <- .51

nparm <- 5
parm <- rep(5,nparm)

am <- numeric(1000)



sg <- matrix(nrow = Iter, ncol = nparm);
sg_ct <- matrix(nrow = Iter - cutf, ncol = nparm)
#Iterates stored

# Sigma Matrix Stored with Square root
sigm <- 1*diag(nparm)
sqrt_sig <- sqrt_mat(sigm)

forb_ebs <- forb_ibs <- numeric(Rep)
ebs_mean <- ibs_mean <- vector("list", Rep)
sqrt_ebs_mean <- sqrt_ibs_mean <- vector("list", Rep)
inv_ebs_mean <- inv_ibs_mean <- vector("list", Rep)
upp_ebs <- low_ebs <- upp_ibs <- low_ibs <- numeric(Rep)
stand_err_ebs <- stand_err_ibs <- numeric(Rep)

#1000  Replications to obtain stable mses
for(cn in 1 : Rep){
  
  #Data Generated 
  
  x <- matrix(rnorm(n*nparm),nrow=n,ncol=nparm)
  x <- x%*%sqrt_sig
  #noisy Observed Data
  y <- x %*% parm + rnorm(n, mean = 0,sd = 1)
  #Learning Rate
  eta <- numeric(Iter)
  sg[1,] <- rep(0, nparm)
  
  for(i in 2 : Iter){
    eta[i] <- i^(-alp)
    sg[i,] <- sg[i-1,] - .5 * eta[i] * grad_lin(sg[i-1,],y[i],x[i,])
    
  }
  
  sg_ct <- sg[(cutf+1):Iter,]
  asg <- colMeans(  sg_ct)
  n <- n - cutf
  #Variance-Covariance Estimate
  ebs_mean[[cn]] <- ebs_batch_mean(sg_ct,alp)
  ibs_mean[[cn]] <- ibs_jasa_mean(sg_ct,alp)
  
  #Norm Difference 
  forb_ebs[cn] <- sqrt(sum((ebs_mean[[cn]]-sigm)^2))
  forb_ibs[cn] <- sqrt(sum((ibs_mean[[cn]]-sigm)^2))
  
  # Square root of the matrix
  sqrt_ebs_mean[[cn]] <- sqrt_mat(ebs_mean[[cn]])
  sqrt_ibs_mean[[cn]] <- sqrt_mat(ibs_mean[[cn]])
  
  # Inverse of the matrix
  inv_ebs_mean[[cn]] <- solve(ebs_mean[[cn]])
  inv_ibs_mean[[cn]] <- solve(ibs_mean[[cn]])
  
  #Standard Error 
  ons <- (rep(1,nparm))
  stand_err_ebs[cn] <- sqrt(t(ons)%*%ebs_mean[[cn]]%*%ons/n)
  stand_err_ibs[cn] <- sqrt(t(ons)%*%ibs_mean[[cn]]%*%ons/n)
  
  #Upper and Lower Ends of the intervals
  upp_ebs[cn] <- t(ons)%*%asg+qnorm(1-qlev/2)*stand_err_ebs[cn]
  low_ebs[cn] <- t(ons)%*%asg-qnorm(1-qlev/2)*stand_err_ebs[cn]
  upp_ibs[cn] <- t(ons)%*%asg+qnorm(1-qlev/2)*stand_err_ibs[cn]
  low_ibs[cn] <- t(ons)%*%asg-qnorm(1-qlev/2)*stand_err_ibs[cn]
}

#Smart batching EBS
print(ebs_mean)

#JASA Chen batching
print(ibs_mean)


print(forb_ebs)
print(forb_ibs)
