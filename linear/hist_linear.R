##########################################################
## The file shows histogram for different estimators
##########################################################
## Function to calculate the SQRT of a matrix using SVD
sqrt_mat <- function(sigm)
{
  
  decomp_sigm <- svd(sigm)
  dig_sig <- decomp_sigm$d
  v_mat <- decomp_sigm$u
  sqrt_sig <- v_mat %*% diag((dig_sig)^(1/2)) %*% t(v_mat)
  
  return(sqrt_sig)
  
}

# Gradient Function for Linear Model
grad_lin <- function(sg,y,x)
{
  (x %*% (sg) - y) %*% x
} 


##############################################
## File contains the function to calculate
## the EBS estimator (with or without lugsail)
##############################################

## Obtains value of bound on MSE of EBS 
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

## EBS estimator function
ebs_batch_hist <- function(sgd, alp = 0.51, cns = 0.1, bet_typ = 1, lug = 1)
{
  n <- nrow(sgd)  # Number of SGD iterates
  nparm <- ncol(sgd)
  tot_mean <- colMeans(sgd)
  
  # Appropriate options for beta 
  if(bet_typ == 1){ bet <- (alp + 1)/2}
  if(bet_typ == 2){ bet <- (2*alp + 1)/3}
  if(bet_typ == 3){ bet <- opt_beta_fn(alp, n)}
  
  #Smart batching
  two_seq <- 2^(seq(10:40))
  bn <- min(two_seq[two_seq >= cns * n^bet]) 	# Equal Batch Size Configuration
  an <- floor(n/bn) 	# No. of batches
  
  ebs <- t(sapply(1 : an, function(i) colMeans( matrix(sgd[(bn*(i - 1) + 1) : (bn*i - 1), ], ncol = nparm ))))
  ebs <- sqrt(bn) * scale(ebs, center = colMeans(sgd), scale = FALSE)
  return(ebs)
}


##############################################
## Function to calculate the IBS estimator 
##############################################

ibs_jasa_hist <- function(sgd, alp = 0.51, cns = 1) 
{
  n <- nrow(sgd)
  nparm <- ncol(sgd)
  
  # IBS Estimators
  am <- ceiling( cns*(1 : 1e3)^(2/(1 - alp)))
  am <- c(am[am < n], (n+1))
  
  # Batch Sizes
  bm <- diff(am)
  batch_means <- t(sapply(1:(length(am)-1), function(i) colMeans( matrix(sgd[am[i] : (am[i+1]-1), ], ncol = nparm ))))
  
  batch_means <- sqrt(bm)*scale(batch_means, center = colMeans(sgd), scale = FALSE)
  return(batch_means)
}



linear_hist <- function( nparm = 5, A = diag(nparm), cns = c(0.1, 1),
                        eta_cns = 0.5, sam_siz = 1e6, qlev = 0.95, alp = .51, 
                        burn_in = 10000, nam_matrix = "indep", cns1 = 0.1)
{
  #sigm <- qr.solve(A) 
  n <- sam_siz
  parm <- seq(1 / nparm, 1, length.out = nparm)  #Parameter space
  
  sg <- matrix(nrow = n + burn_in, ncol = nparm) 
  sg_ct <- matrix(nrow = n, ncol = nparm) # sgd iterates after dropping burn-in samples

  # Following matrix will be used in generating regressors X
  sqrt_sig <- diag(nparm)
  
  # Declaring variable names  
  cns_ln <- 3*length(cns)    # 3 types of beta 
  
  
  # Generating data of Maximum Sample Size
  x <- matrix(rnorm((n + burn_in) * nparm), nrow = (n + burn_in), ncol = nparm)
  x <- x %*% sqrt_sig
  y <- x %*% parm + rnorm((n + burn_in), mean = 0, sd = 1)
  
  # Learning Rate
  eta <- numeric(n + burn_in)
  
  sg[1, ] <- rep(0, nparm)  # always starting from the zero vector estimate
  
  ## Generating full SGD sequence -- this takes time
  for(i in 2 : (n + burn_in))
  {
    eta[i] <- i^( - alp)
    sg[i, ] <- sg[i - 1, ] - eta_cns * eta[i] * grad_lin(sg[i - 1, ], y[i], x[i, ]) 
  }
  sg_ct_full <- sg[(burn_in + 1) : (n + burn_in), ]      # Removing burn in
    
  sg_ct <- sg_ct_full # main process inside this loop
  asg <- colMeans(sg_ct) # Averaged SGD
  
    ################################
    ### IBS batch means calculations
    ################################
    ibs_hist     <- as.vector(ibs_jasa_hist(sg_ct, alp, cns = cns1))
  
  pdf(file = paste("out/Hist_qq_dim_", nparm, "_sam_", n, ".pdf"))
  par(mfrow = c(2, 2))
    hist(ibs_hist, main = "IBS")
    qqnorm(ibs_hist, pch = 1, frame = FALSE, main = "IBS")
    qqline(ibs_hist, col = "steelblue", lwd = 2)
    # Different settings of equal batch size estimator (EBS), for values of cns and three types of beta
    for(mk in 1 : 1)
    { 
      for(bt_typ in 2 : 2)
      { 
        ###########
        ### EBS calculation 
        ###########
        ebs_hist <- as.vector(ebs_batch_hist(sg_ct, alp, cns[mk], bt_typ))
        hist(ebs_hist, main = "EBS")
        qqnorm(ebs_hist, pch = 1, frame = FALSE, main = "EBS")
        qqline(ebs_hist, col = "steelblue", lwd = 2)
      }
    }  
  dev.off()
}