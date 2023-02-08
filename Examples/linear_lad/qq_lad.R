# QQ plot for different estimators

sqrt_mat <- function(sigm)
{
  
  decomp_sigm <- svd(sigm)
  dig_sig     <- decomp_sigm$d
  v_mat       <- decomp_sigm$u
  sqrt_sig    <- v_mat %*% diag((dig_sig)^(1/2)) %*% t(v_mat)
  
  return(sqrt_sig)
  
}

# Gradient Function for lad Model
grad_lad <- function(sg,y,x, tau = 0.5)
{
  (tau - as.numeric( y - t(x) %*% sg < 0 )) %*% x
} 


## Obtains value of bound on MSE of EBS 
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
  opt_beta <- optim(.6, fn = b_our, lower = alpha, upper = 1, n = m, 
                    alpha1 = alpha, method = "L-BFGS-B")$par
  return(opt_beta)
}

## EBS estimator function
ebs_batch <- function(sgd, alp = 0.51, cns = 0.1, bet_typ = 1, lug = 1)
{
  n <- nrow(sgd)  
  nparm <- ncol(sgd)
  tot_mean <- colMeans(sgd)
  
  if(bet_typ == 1){ bet <- (alp + 1)/2}
  if(bet_typ == 2){ bet <- (2*alp + 1)/3}
  if(bet_typ == 3){ bet <- opt_beta_fn(alp, n)}
  
  #Smart batching
  two_seq <- 2^(seq(10:40))
  bn      <- min(two_seq[two_seq >= cns * n^bet]) 
  an      <- floor(n/bn) 
  
  ebs <- t(sapply(1 : an, function(i) colMeans( 
         matrix(sgd[(bn*(i - 1) + 1) : (bn*i), ],  ncol = nparm ))))
  ebs <- sqrt(bn) * scale(ebs, center = 
                            colMeans(sgd), scale = FALSE)
  return(ebs)
}


# Function to calculate the IBS estimator 

ibs_jasa <- function(sgd, alp = 0.51, cns = 1) 
{
  n     <- nrow(sgd)
  nparm <- ncol(sgd)
  
  # IBS Estimators
  am <- ceiling( cns*(1 : 1e3)^(2/(1 - alp)))
  am <- c(am[am < n], (n+1))
  
  # Batch Sizes
  bm <- diff(am)

  batch_means <- t(sapply(1:(length(am)-1), function(i) colMeans(
                 matrix(sgd[am[i] : (am[i+1]-1), ], ncol = nparm ))))
  
  batch_means <- sqrt(bm)*scale(batch_means, center = colMeans(sgd), 
                                scale = FALSE)
  return(batch_means)
}



lad_qq <- function( nparm = 5, cns = c(0.1, 1),
                        eta_cns = 0.5, sam_siz = c(1e6), alp = .51, 
                        burn_in = 1000, cns1 = 0.1, Iter = 2, nsampl = 1000)
{

  n    <- sam_siz
  parm <- seq(1 / nparm, 1, length.out = nparm)
  
  sg    <- matrix(nrow = n + burn_in, ncol = nparm) 
  sg_ct <- matrix(nrow = n, ncol = nparm) 

  
  sqrt_sig    <- diag(nparm)
  ibs_batches <- ebs_batches <- vector()
  
  for( itr in 1 : Iter)
  {  
    
    x <- matrix(rnorm((n + burn_in) * nparm), nrow = (n + burn_in)
                , ncol = nparm)
    x <- x %*% sqrt_sig
    y <- x %*% parm + rdoublex(n + burn_in, mu = 0, lambda = 1)
    
    eta     <- numeric(n + burn_in)
    sg[1, ] <- rep(0, nparm)  
    
    ## Generating full SGD sequence -- this takes time
    for(i in 2 : (n + burn_in))
    {
      eta[i]  <- i^( - alp)
      sg[i, ] <- sg[i - 1, ] + eta_cns * eta[i] * grad_lad(
                sg[i - 1, ], y[i], x[i, ])
    }
    sg_ct_full <- sg[(burn_in + 1) : (n + burn_in), ]     
      
    sg_ct <- sg_ct_full 
    asg   <- colMeans(sg_ct) 
    
      # IBS batch means calculations
      
      ibs_batches <- c(ibs_batches, 
                     as.vector(ibs_jasa(sg_ct, alp, cns = cns1)))
    
      # EBS calculation 
     
      ebs_batches <- c(ebs_batches, 
                     as.vector(ebs_batch(sg_ct, alp, cns = cns, bet_typ = 2)))
         
      
  }

  pdf(file = paste("out/qq_dim_lad",nparm,"_sam_",n,".pdf", sep = ""),
  height = 5, width = 5)
   
  p1 <- pnorm(seq(-.5,.5, length = nsampl))
  ibs_new <- quantile(scale(ibs_batches), probs = p1)
  ebs_new <- quantile(scale(ebs_batches), probs = p1)
  theory <- qnorm(p1)
  plot(theory, ibs_new, xlab = "Theoretical Quantiles", 
    ylab = "Sample Quantiles", ylim = range(theory), type = "n")
  lines(theory, ibs_new, pch = 1, col = "steelblue")
  lines(theory, ebs_new, pch = 1, col = "brown", lwd = 2)
  lines(theory, theory, col = "black", lwd = 1, lty = 2)
  legend("bottomright", legend = c("IBS", "EBS"), 
    col = c("steelblue", "brown"), 
    lty = 1, box.lty=0, box.lwd=1)
  
  dev.off()
}