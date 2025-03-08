# Obtains value of bound on MSE of EBS 
b_our <- function(beta1, n, alpha1)
{
  bn <- n^(beta1)
  an <- n/bn
  return( log(n^(-alpha1/4) + an^(-1/2) + bn^(alpha1-1) 
              + bn^(-1/2) * n^(alpha1/2) + an^(-1) + n^(-2*alpha1) * bn) )
}

## Optimal value of beta for alpha and sample size
opt_beta_fn <- function(alpha, m)
{ 
  opt_beta <- optim(.6, fn = b_our, lower = alpha, upper = 1, 
                    n = m, alpha1 = alpha, method = "L-BFGS-B")$par
  return(opt_beta)
}

## EBS estimator function
ebs_batch_mean <- function(sgd, alp = 0.51, cns = 0.1, 
                           bet_typ = 1, lug = 1, poly = 0)
{
	n    <- nrow(sgd)  
	nparm <- ncol(sgd)
	tot_mean <- colMeans(sgd)


	if(bet_typ == 1){ bet <- (alp + 1)/2}
	if(bet_typ == 2){ bet <- (2*alp + 1)/3}
	if(bet_typ == 3){ bet <- opt_beta_fn(alp, n)}

	if(poly == 1){
	  #polynomial batching
	  bn <- floor(cns * n^bet)
	  an <- floor(n/bn) 	
	}else {
	#Smart batching
	two_seq <- 2^(seq(10:40))
	bn      <- min(two_seq[two_seq >= cns * n^bet]) 
	an      <- floor(n/bn) 	
}
	# EBS is simple call to mcse.multi
	ebs <- mcse.multi(sgd, size = bn, r = lug)$cov

	# Adjustment for remaining samples (an*bn to n) since mcse.multi chops the end
	if(an*bn + 1 < n)	
	{
	  ebs     <- ebs * (an - 1)/an
	  add_trm <- (colSums(sgd[(an*bn + 1) : n, ]) - (n - an*bn) * tot_mean)
		ebs     <- ebs + add_trm %*% t(add_trm)/n
  }
	return(ebs)
}
