rm(list=ls())
source("ibs_lng.R")

alp <- .51
sq <- seq(5,20,by=5)

leng_ibs <- matrix(nrow=5,ncol=length(sq))
leng_ebs <- matrix(nrow=5,ncol=length(sq))
for( k in 1:length(sq)){
  nparm <- sq[k]
  leng_ibs[1,k] <- ibs_lng(nparm,cns=0.01,alp)
  leng_ibs[2,k] <- ibs_lng(nparm,cns=0.1,alp)
  leng_ibs[3,k] <- ibs_lng(nparm,cns=1,alp)
  leng_ibs[4,k] <- ibs_lng(nparm,cns=2,alp)
  leng_ibs[5,k] <- ibs_lng(nparm,cns=5,alp)
#Number of batches in IBS scheme
  
}



