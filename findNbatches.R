b_our <- function(beta1, n, alpha1)
{
  bn <- n^(beta1)
  an <- n/bn
  return( log(n^(-alpha1/4) + an^(-1/2) + bn^(alpha1-1) + bn^(-1/2) * n^(alpha1/2) + an^(-1) + n^(-2*alpha1) * bn) )
}
## Optimal value of beta for alpha and sample size
opt_beta_fn <- function(alpha, m)
{ 
  opt_beta <- optim(.6, fn = b_our, lower = alpha, upper = 1, n = m, alpha1 = alpha, method = "L-BFGS-B")$par
  return(opt_beta)
}

nbatches <- function(n, alp = 0.51, cns = 1, cns1 = 1)
{
  num_batc <- numeric(length(4))
  am <- ceiling( cns1*(1 : 1e3)^(2/(1 - alp)))
  am <- c(am[am < n], (n+1))
  num_batc[1] <- length(am)

  bet <- (alp + 1)/2
  two_seq <- 2^(seq(10:40))
  bn <- min(two_seq[two_seq >= cns * n^bet])    # Equal Batch Size Configuration
  an <- floor(n/bn)    # No. of batches
  num_batc[2] <- an

   bet <- (2*alp + 1)/3
   bn <- min(two_seq[two_seq >= cns * n^bet])   # Equal Batch Size Configuration
   an <- floor(n/bn)   # No. of batches
   num_batc[3] <- an

   bet <- opt_beta_fn(alp, n)
   two_seq <- 2^(seq(10:40))
   bn <- min(two_seq[two_seq >= cns * n^bet])   # Equal Batch Size Configuration
   an <- floor(n/bn)   # No. of batches
  num_batc[4] <- an

  names(num_batc) <- c("ibs", "beta1", "beta2", "beta3")
  return(num_batc)
}

nbatches(n = 1e6, cns = 0.1, cns1 = .1)
