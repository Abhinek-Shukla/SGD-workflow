ebs_batch_mean <- function(sg_ct,alp){
  n <- nrow(sg_ct)#Number of SGD iterates
  nparm <- ncol(sg_ct)
library(mcmcse)
bet <- 0.5414 + 0.4669*alp - 0.30195/log10(n)
two_seq <- 2^(seq(10:40))
#Equal Batch Size Configuration
bn <- min(two_seq[two_seq >= n^bet])
#No. of batches
an <- floor(n/bn)


tot_mean <- colMeans(sg_ct)
ebs <- matrix(rep(0,nparm^2),nrow=nparm,ncol=nparm)

ebs <- mcse.multi(sg_ct,size=bn,r=1)$cov
ebs <- ebs*(an-1)*bn/n

add_trm <- numeric(nparm)
add_trm <- (colSums(sg_ct[(an*bn+1):n,])-(n-an*bn)*tot_mean)

ebs <- ebs + add_trm%*%t(add_trm)/n
return(ebs)
}
