rm(list=ls())
setwd("C:/Users/Hp/Documents/GitHub/Batch_Means_Online/check_value_constant")
source("ibs_lng.R")
source("ebs_lng.R")
alp <- .51
cns_sq <- c(0.1,1)
sq_n <- c(5e4, 1e5, 2e5, 5e5, 8e5, 1e6, 5e6, 1e7 )
leng_ebs <- matrix(nrow = length(cns_sq), ncol = length(sq_n))
leng_ibs <- matrix(nrow = length(cns_sq), ncol = length(sq_n))


for( k in 1:length(sq_n)){
  n <- sq_n[k]
  for (m in 1 : length(cns_sq)){
    
    leng_ebs[m,k] <- ebs_lng(n,cns=cns_sq[m],alp,bet_typ=1)
    
    leng_ibs[m,k] <- ibs_lng(n,cns=cns_sq[m],alp)
  }
  
}
print(leng_ebs)

print(leng_ibs)