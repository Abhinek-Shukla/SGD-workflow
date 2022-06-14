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
volm_ebs <- volm_ibs <- numeric(Rep)
cover_ebs <- cover_ibs <- cover_orc <- numeric(Rep)
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
  ebs_mean <- ebs_batch_mean(sg_ct,alp)
  ibs_mean <- ibs_jasa_mean(sg_ct,alp)
  
  #Norm Difference 
  forb_ebs[cn] <- sqrt(sum((ebs_mean-sigm)^2))/sqrt(sum(sigm^2))
  forb_ibs[cn] <- sqrt(sum((ibs_mean-sigm)^2))/sqrt(sum(sigm^2))
  
  # Volume of the matrix
  volm_ebs[cn] <- (det(ebs_mean))^(1/nparm)
  volm_ibs[cn] <- (det(ibs_mean))^(1/nparm)
  
#critical value calculation
crt_val <- qchisq(0.95,df=nparm)

cover_ebs[cn] <- as.numeric(n*t(asg-parm)%*%solve(ebs_mean)%*%(asg-parm)<=crt_val)
cover_ibs[cn] <- as.numeric(n*t(asg-parm)%*%solve(ibs_mean)%*%(asg-parm)<=crt_val)
cover_orc[cn] <- as.numeric(n*t(asg-parm)%*%solve(sigm)%*%(asg-parm)<=crt_val)  

}

