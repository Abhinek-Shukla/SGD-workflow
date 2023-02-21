#######################################
## File runs simulation for qq_plots
## for the batch means from ebs and ibs
## estimators and compares to normality
#######################################
source("linear_back_funcs.R")

linear_qq <- function(nparm = 5, cns = c(0.1, 1),
                        eta_cns = 0.5, sam_siz = c(1e6), 
                        alp = .51, burn_in = 1000, cns1 = 0.1, 
                        Iter = 2, nsampl = 1000)
{

  n    <- sam_siz
  parm <- seq(1/nparm, 1, length.out = nparm)  
  sqrt_sig <- diag(nparm)
  
  # memory allocation is challenging here
  ibs_batches <- vector() 
  ebs_batches <- vector()
  
  for(itr in 1:Iter)
  {  
    
    sg_ct <- linear_sgd(n = n, burn_in = burn_in, parm = parm, alp = alp)
    asg <- colMeans(sg_ct)

    # batch means calculations
    ibs_batches <- c(ibs_batches, as.vector(ibs_jasa(sg_ct, alp, cns = cns1)))
    ebs_batches <- c(ebs_batches, as.vector(ebs_batch(sg_ct, alp, cns = cns)))  
  }

  # Save plot 
  pdf(file = paste("out/qq_dim_", nparm, "_sam_", n, ".pdf"),
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