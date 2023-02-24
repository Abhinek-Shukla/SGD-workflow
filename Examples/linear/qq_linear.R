#######################################
## File runs simulation for qq_plots
## for the batch means from ebs and ibs
## estimators and compares to normality
#######################################
source("linear_back_funcs.R")

linear_qq <- function(nparm = 5, cns = c(0.1, 1),
                        eta_cns = 0.5, sam_siz = c(1e6), 
                        alp = .51, burn_in = 1000, cns1 = 0.1, 
                        Iter = 2)
{

  n    <- sam_siz
  parm <- seq(1/nparm, 1, length.out = nparm)  
  sqrt_sig <- diag(nparm)
  
  # memory allocation is challenging here
  ibs_batches <- vector() 
  ebs_batches <- vector()
  
  for(itr in 1:Iter)
  {  
    
    sg_ct <- linear_sgd(n = n, burn_in = burn_in, parm = parm, 
                alp = alp, sqrt_sig = sqrt_sig, eta_cns = eta_cns)
    asg <- colMeans(sg_ct)

    # batch means calculations
    ibs_batches <- c(ibs_batches, as.vector(ibs_jasa(sg_ct, alp, cns = cns1)))
    ebs_batches <- c(ebs_batches, as.vector(ebs_batch(sg_ct, alp, cns = cns)))  
  }

  file <- paste0("out/qq_dim_", nparm, "_sam_", n, ".Rdata")
  save(ibs_batches, ebs_batches, file = file)
}