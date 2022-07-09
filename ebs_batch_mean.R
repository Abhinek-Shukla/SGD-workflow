library(mcmcse)

ebs_batch_mean <- function(sgd, alp = 0.51, cns = 0.1, bet_typ = 1, lug = 1)
{
	n <- nrow(sgd)#Number of SGD iterates
	nparm <- ncol(sgd)

	if(bet_typ == 1){ bet <- (alp + 1)/2}
	if(bet_typ == 2){ bet <- (2*alp + 1)/3}
	if(bet_typ == 3){ bet <- 0.5414 + 0.4669*alp - 0.30195/log10(n)}

	two_seq <- 2^(seq(10:40))
	
	#Equal Batch Size Configuration
	bn <- min(two_seq[two_seq >= cns*n^bet])
	#No. of batches
	an <- floor(n/bn)
	
	tot_mean <- colMeans(sgd)

	ebs <- mcse.multi(sgd, size = bn, r = lug)$cov
	ebs <- ebs*(an-1)*bn/n

	add_trm <- (colSums(sgd[(an*bn+1):n, ]) - (n-an*bn)*tot_mean)

	ebs <- ebs + add_trm%*%t(add_trm)/n
	return(ebs)
}
