rm(list=ls())
library(MASS)
source("grad_lin.R")
source("ebs_batch_mean.R")
source("ibs_jasa_mean.R")
Rep <- 1
#Sample Size
n <- 5e4;

#Iterations
Iter <- n;
alp <- .51

nparm <- 5
parm <- rep(5,nparm)

am <- numeric(1000)
cutf <- 1000 #Dropping initial Iterates of SGD


sg <- matrix(nrow = Iter, ncol = nparm);
sg_ct <- matrix(nrow = Iter - cutf, ncol = nparm)
#Iterates stored

# Sigma Matrix Stored with Square root
sigm <- 1*diag(nparm)
decomp_sigm<- svd(sigm)
dig_sig <- decomp_sigm$d
v_mat <- decomp_sigm$u
sqrt_sig <- v_mat%*%diag(sqrt(dig_sig))%*%t(v_mat)

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

  ebs_mean <- ebs_batch_mean(sg_ct,alp)
 
  ibs_mean <- ibs_jasa_mean(sg_ct,alp)
  
  
}

#Smart batching EBS
print(ebs_mean)

#JASA Chen batching
print(ibs_mean)

