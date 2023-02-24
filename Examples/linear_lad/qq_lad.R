#######################################
## File runs simulation for qq_plots
## for the batch means from ebs and ibs
## estimators and compares to normality
#######################################
source("lad_back_funcs.R")

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