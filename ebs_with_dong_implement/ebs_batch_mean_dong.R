library(mcmcse)

ebs_batch_mean_dong <- function(sgd, an, bn, alp = 0.51,  lug = 1)
{
	n <- nrow(sgd)#Number of SGD iterates
	nparm <- ncol(sgd)

print(dim(sgd))
	ebs <- mcse.multi(sgd, size = bn, r = lug)$cov

	ebs <- ebs*(an-1)/(an)

	return(ebs)
}
