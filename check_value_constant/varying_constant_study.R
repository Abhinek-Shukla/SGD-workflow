rm(list=ls())
setwd("C:/Users/Hp/Documents/GitHub/Batch_Means_Online/check_value_constant")
source("ibs_lng.R")
source("ebs_lng.R")
alp <- .51
sq <- seq(5,100,by=5)

leng_ibs <- matrix(nrow=5,ncol=length(sq))

for( k in 1:length(sq)){
  nparm <- sq[k]
  leng_ibs[1,k] <- ibs_lng(nparm,cns=0.01,alp)
  leng_ibs[2,k] <- ibs_lng(nparm,cns=0.1,alp)
  leng_ibs[3,k] <- ibs_lng(nparm,cns=1,alp)
  leng_ibs[4,k] <- ibs_lng(nparm,cns=2,alp)
  leng_ibs[5,k] <- ibs_lng(nparm,cns=5,alp)
#Number of batches in IBS scheme
  
}

sq_n <- seq(1e3,1e7,by=1e3)
leng_ebs <- matrix(nrow=5,ncol=length(sq_n))
for( k in 1:length(sq_n)){
  n <- sq_n[k]
  leng_ebs[1,k] <- ebs_lng(n,cns=0.01,alp)
  leng_ebs[2,k] <- ebs_lng(n,cns=0.1,alp)
  leng_ebs[3,k] <- ebs_lng(n,cns=1,alp)
  leng_ebs[4,k] <- ebs_lng(n,cns=2,alp)
  leng_ebs[5,k] <- ebs_lng(n,cns=5,alp)

  
}
