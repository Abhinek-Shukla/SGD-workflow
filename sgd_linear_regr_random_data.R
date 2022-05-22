rm(list=ls())
library(MASS)
Rep <- 1
#Sample Size
n <- 5e6;
#Iterations
Iter <- n;
alp <- .51
bet <- 0.5416+0.4671*alp-0.6930/log10(n)#(1+alp)/2
nparm <- 5
parm <- rep(5,nparm)
two_seq <- 2^(seq(10:40))
am <- numeric(1000)
cutf <- 1000 #Dropping initial Iterates of SGD
#Equal Batch Size Configuration
bn <- min(two_seq[two_seq >= (Iter - cutf)^bet])
#No. of batches
an <- floor((Iter - cutf)/bn)

sg <- matrix(nrow = Iter, ncol = nparm);
sg_ct <- matrix(nrow=Iter-cutf,ncol=nparm)
#Iterates stored


#1000  Replications to obtain stable mses
for(cn in 1 : Rep){
  
  
  
  #Learning Rate
  eta <- vector()
  
  #Data Generated 
  x <- mvrnorm(n, mu = rep(0, nparm), Sigma = diag(nparm))
  
  #noisy Observed Data
  y <- x%*%parm + rnorm(n, mean = 0,sd = 1)
  
  
  sg[1,] <- rep(0, nparm)
  
  for(i in 2 : Iter){
    eta[i] <- i^(-alp)
    sg[i,] <- sg[i-1,] - .5 * eta[i] * (x[i,] %*% (sg[i-1,]) - y[i]) %*% x[i,]
  }

  sg_ct <- sg[(cutf+1):Iter,]
  tot_mean <- rep(0,nparm)
ebs <- matrix(rep(0,nparm*an),nrow=an,ncol=nparm)
#Equal batch size smart batching (EBS)
for(k in 1 : an){
  strt_pt <- (k-1)*bn + 1
  end_pt <- bn*k
  ebs[k,] <-  colMeans(sg_ct[strt_pt : end_pt,])
  tot_mean <- tot_mean + colSums(sg_ct[strt_pt:end_pt,])
}
tot_mean <- tot_mean/(an*bn)
ebs_mean <- matrix(rep(0,nparm^2),nrow=nparm,ncol=nparm)
for(k in 1 : an){
  ebs_mean <- ebs_mean + (ebs[k,]-tot_mean)%*%t((ebs[k,]-tot_mean))
  
}
ebs_mean <- ebs_mean*bn/an


# JASA Online BM Estimators
for(m in 1:1000){
  
  am[m] <- floor(m^(2/(1-alp)))
  
}
am <- c(am[am <= Iter-cutf],Iter-cutf)

ibs_jasa <- matrix(rep(0,nparm*(length(am)-1)),nrow=(length(am)-1),ncol=nparm)
tot_mean <- rep(0,nparm)
#Equal batch size smart batching (EBS)
for(k in 1 : (length(am)-1)){
  strt_pt <- am[k]
  end_pt <- am[k+1]-1
  if(k==(length(am)-1)){  end_pt <- am[k+1]}
  ibs_jasa[k,] <-  colSums(sg_ct[strt_pt : end_pt,])
  tot_mean <- tot_mean + colSums(sg_ct[strt_pt:end_pt,])
}
tot_mean <- tot_mean/(Iter-cutf)
ibs_jasa_mean <- matrix(rep(0,nparm^2),nrow=nparm,ncol=nparm)

for(k in 1 : (length(am)-1)){

   tmp_vl <- (ibs_jasa[k,]-(am[k+1]-am[k])*tot_mean)
   if(k==(length(am)-1)){   tmp_vl <- (ibs_jasa[k,]-(am[k+1]-am[k]+1)*tot_mean)}
   ibs_jasa_mean <- ibs_jasa_mean + tmp_vl%*%t(tmp_vl)
  
}
ibs_jasa_mean <- ibs_jasa_mean/(Iter-cutf)


}

#Smart batching EBS
print(ebs_mean)

#JASA Chen batching
print(ibs_jasa_mean)