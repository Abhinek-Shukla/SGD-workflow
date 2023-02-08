# Miscellaneous functions used in the paper

sqrt_mat <- function(sigm)
{
  
  decomp_sigm <- svd(sigm)
  dig_sig     <- decomp_sigm$d
  v_mat       <- decomp_sigm$u
  sqrt_sig    <- v_mat %*% diag((dig_sig)^(1/2)) %*% t(v_mat)

  return(sqrt_sig)

}


# Covmatrix matrix of dim nparm
# Of kind Toeplitz and and equivariance
X_mat <- function(nparm, rho)
{
  equiv <- matrix(rho, nrow = nparm, ncol = nparm)
  diag(equiv) <- rep(1, nparm)

  toep <- matrix(rho, nrow = nparm, ncol = nparm)
  toep <- toep^(abs(col(toep) - row(toep)))

  return(list("equiv" = equiv, "toep" = toep))
}




# Simultaneous confidence intervals 

# Function is from Assessing simultaneous paper by Robertson et.al
# JCGS. Code from JCGS website

###############################################################

## Function produces simultaenous intervals
## n.sim  = to calculate the degrees of freedom
## Sigma = estimated asymptotic covariance matrix
## conf = confidence level for simultaneous regions
## center = the mean estimate around which to create intervals
## epsilon = tolerance level for SI algorithm
## n.star = maximum number of iterations for the algorithm

###############################################################


new.sim.int <-  function(Sigma, conf = .90, center = rep(0,dim(Sigma)[1]),
                epsilon = .001,  n.star = 100)
{

  p          <- dim(Sigma)[1]
  crit.under <- qnorm(1-(1-conf)/2)
  crit.over  <- qnorm(1-(1-conf)/2/p)

  i <- 1
  LB.i  <- crit.under*sqrt(diag(Sigma))
  UB.i  <- crit.over*sqrt(diag(Sigma))
  X.i   <- (LB.i+UB.i)/2
  F.X.i <- pmvnorm(lower=c(-X.i),upper=c(X.i),sigma=Sigma)

  while(abs(F.X.i[1] - conf) > epsilon & i <= n.star)
  {
    i = i + 1
    ## equality should never happen since then it would terminate?
    if(F.X.i[1] - conf < 0)
    {
      LB.i <- X.i
      UB.i <- UB.i
      X.i <- (X.i + UB.i)/2
    }
    if(F.X.i[1] - conf > 0)
    {
      LB.i <- LB.i
      UB.i <- X.i
      X.i <- (X.i + LB.i)/2
    }
    F.X.i <- pmvnorm(lower = c(-X.i), upper = c(X.i), sigma = Sigma)
  }

  out <- list("prob_level_error" = F.X.i[1] - conf, "iterations" = i,
    "pmvt_error" = attributes(F.X.i)$error,
    "ints" = cbind(center - X.i, center + X.i))
  return(out)
}