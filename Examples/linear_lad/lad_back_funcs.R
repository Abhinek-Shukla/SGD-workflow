#######################################
## Function that generates linear
## regression SGD iterates. Plus a few 
## other functions required in the code
#######################################

# Gradient Function for LAD Model
grad_lad <- function(sg,y,x, tau = 0.5)
{
  (tau - as.numeric( y - t(x) %*% sg < 0 )) %*% x
} 

lad_sgd <- function(n, burn_in, parm, nparm = length(parm), 
						alp, sqrt_sig, eta_cns)
{
    # Generating data of Maximum Sample Size
    x <- matrix(rnorm((n + burn_in) * nparm),
                nrow = (n + burn_in), ncol = nparm)
    x <- x %*% sqrt_sig
    y <- x %*% parm + rdoublex(n + burn_in, mu = 0,lambda = 1)
    
    # Learning Rate
    eta   <- numeric(n + burn_in)
    sg    <- matrix(nrow = n + burn_in, ncol = nparm)
    sg[1, ] <- rep(0, nparm)  

    ## Generating full SGD sequence -- this takes time
    for(i in 2:(n + burn_in))
    {
      eta[i] <- eta_cns*i^(-alp)
      sg[i, ] <- sg[i-1, ] + eta[i]*grad_lad(sg[i-1, ], y[i], x[i, ]) 
    }
    sg_ct_full <- sg[(burn_in + 1) : (n + burn_in), ]
    sg_ct <- sg[(burn_in + 1) : (n + burn_in), ]     
    return(sg_ct)
}



# Calculate IBS batch means only (not the batch-means estimator)
ebs_batch <- function(sgd, alp = 0.51, cns = 0.1)
{
  n <- nrow(sgd)  
  nparm <- ncol(sgd)
  tot_mean <- colMeans(sgd)
  
  # Appropriate options for beta 
  bet <- (alp + 1)/2

  #Smart batching
  two_seq <- 2^(seq(10:40))
  bn <- min(two_seq[two_seq >= cns * n^bet]) 
  an <- floor(n/bn) 
  
  ebs <- t(sapply(1 : an, function(i) colMeans( 
                  matrix(sgd[(bn*(i - 1) + 1) : (bn*i), ], ncol = nparm ))))
  ebs <- sqrt(bn) * scale(ebs, center = colMeans(sgd), scale = FALSE)
  return(ebs)
}


# Calculate IBS batch means only (not the batch-means estimator)
ibs_jasa <- function(sgd, alp = 0.51, cns = 1) 
{
  n <- nrow(sgd)
  nparm <- ncol(sgd)
  
  am <- ceiling( cns*(1 : 1e3)^(2/(1 - alp)))
  am <- c(am[am < n], (n+1))
  
  # Batch Sizes
  bm <- diff(am)

  batch_means <- t(sapply(1:(length(am)-1), function(i) colMeans(
    matrix(sgd[am[i] : (am[i+1]-1), ], ncol = nparm ))))
  
  batch_means <- sqrt(bm)*scale(batch_means, center = 
                                  colMeans(sgd), scale = FALSE)
  return(batch_means)
}

