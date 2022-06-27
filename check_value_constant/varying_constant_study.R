rm(list=ls())
setwd("C:/Users/Hp/Documents/GitHub/Batch_Means_Online/check_value_constant")
source("ibs_lng.R")
source("ebs_lng.R")
alp <- .51

sq_n <- seq(1e3,1e6,by=5e3)
leng_ebs <- matrix(nrow=5,ncol=length(sq_n))
for( k in 1:length(sq_n)){
  n <- sq_n[k]
  leng_ebs[1,k] <- ebs_lng(n,cns=0.01,alp)
  leng_ebs[2,k] <- ebs_lng(n,cns=0.1,alp)
  leng_ebs[3,k] <- ebs_lng(n,cns=1,alp)
  leng_ebs[4,k] <- ebs_lng(n,cns=2,alp)
  leng_ebs[5,k] <- ebs_lng(n,cns=5,alp)

  
}

leng_ibs <- matrix(nrow=5,ncol=length(sq_n))

for( k in 1:length(sq_n)){
  nparm <- leng_ebs[1,k]
  leng_ibs[1,k] <- ibs_lng(nparm,cns=0.01,alp)
  nparm <- leng_ebs[2,k]
  leng_ibs[2,k] <- ibs_lng(nparm,cns=0.1,alp)
  nparm <- leng_ebs[3,k]
  leng_ibs[3,k] <- ibs_lng(nparm,cns=1,alp)
  nparm <- leng_ebs[4,k]
  leng_ibs[4,k] <- ibs_lng(nparm,cns=2,alp)
  nparm <- leng_ebs[5,k]
  leng_ibs[5,k] <- ibs_lng(nparm,cns=5,alp)
  #Number of batches in IBS scheme
  
}


sampl_size_log <- log10(sq_n)
plot(leng_ebs[1,],sampl_size_log,type="l")
lines(eng_ebs[1,],log10(leng_ibs[1,]))
