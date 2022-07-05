rm(list=ls())
library(MASS)
library(matrixcalc)
source("grad_lin.R")
source("ebs_batch_mean.R")
source("ibs_jasa_mean.R")
source("sqrt_mat.R")

Rep <- 1000
cutf <- 1000 #Dropping initial Iterates of SGD
#Sample Size
n <- 1e5+cutf;
#Confidence level 
qlev <- 0.95
#Iterations

alp <- .51

nparm <- 5
parm <- seq(1/nparm,1,length.out=nparm)
crt_val <- qchisq(qlev,df=nparm)
am <- numeric(1000)
Iter <- n;


sg <- matrix(nrow = Iter, ncol = nparm);
sg_ct <- matrix(nrow = Iter - cutf, ncol = nparm)
#Iterates stored

# Sigma Matrix Stored with Square root
sigm <- 1*diag(nparm)
sqrt_sig <- sqrt_mat(sigm)

forb_ibs <- volm_ibs <- cover_ibs <- forb_ibs_norm <- numeric(Rep)
cover_orc<- numeric(Rep)

volm_ebs <- forb_ebs <- forb_ebs_norm <- matrix(rep(0,6*Rep),nrow = Rep,ncol = 6)
cover_ebs  <- matrix(rep(0,6*Rep),nrow = Rep,ncol = 6)
volm_ebs_ls <- forb_ebs_ls <- forb_ebs_norm_ls <- matrix(rep(0,6*Rep),nrow = Rep,ncol = 6)
cover_ebs_ls  <- matrix(rep(0,6*Rep),nrow = Rep,ncol = 6)
cns = c(0.1,1)
#1000  Replications to obtain stable results
for(cn in 1 : Rep){
  if(cn>=2){n <- n+cutf}
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
  
  #IBS and Oracle related coverages and volume
  
  ibs_mean     <- ibs_jasa_mean(sg_ct,alp)
  forb_ibs[cn]  <- sqrt(sum((ibs_mean-sigm)^2))/sqrt(sum(sigm^2))
  forb_ibs_norm[cn] <- sqrt(sum((ibs_mean)^2))
  volm_ibs[cn]  <- (det(ibs_mean))^(1/nparm)
  cover_ibs[cn] <- as.numeric(n*t(asg - parm)%*%qr.solve(ibs_mean)%*%(asg - parm) <= crt_val)
  cover_orc[cn] <- as.numeric(n*t(asg - parm)%*%solve(sigm)%*%(asg - parm) <= crt_val)  
  count=1
  #Nine different settings of EBS, c=0.1,1,2
  for( mk in 1:length(cns)){
    for(bk in 1:3){
      ebs_mean            <- ebs_batch_mean(sg_ct,alp,cns[mk],bk,1)
      forb_ebs_norm[cn,count] <- sqrt(sum((ebs_mean)^2))
      
      forb_ebs[cn,count]  <- sqrt(sum((ebs_mean-sigm)^2))/sqrt(sum(sigm^2))
      volm_ebs[cn,count]  <- (det(ebs_mean))^(1/nparm)
      cover_ebs[cn,count] <- as.numeric(n*t(asg-parm)%*%qr.solve(ebs_mean)%*%(asg-parm)<=crt_val)
      
      ebs_mean            <- ebs_batch_mean(sg_ct,alp,cns[mk],bk,2)
      forb_ebs_norm_ls[cn,count] <- sqrt(sum((ebs_mean)^2))
      
      forb_ebs_ls[cn,count]  <- sqrt(sum((ebs_mean-sigm)^2))/sqrt(sum(sigm^2))
      volm_ebs_ls[cn,count]  <- (det(ebs_mean))^(1/nparm)
      cover_ebs_ls[cn,count] <- as.numeric(n*t(asg-parm)%*%qr.solve(ebs_mean)%*%(asg-parm)<=crt_val)
      
      count=count+1           
    }
  }
  
  
}

save(forb_ibs,forb_ebs,forb_ebs_ls,file="forb_details_1e5_lin.RData")
save(volm_ibs,volm_ebs,volm_ebs_ls,file="volm_details_1e5_lin.RData")
save(cover_orc,cover_ibs,cover_ebs,cover_ebs_ls,file="cover_details_1e5_lin.RData")
save(forb_ibs_norm,forb_ebs_norm, forb_ebs_norm_ls,file="forb_norm_details_1e5_lin.RData")

