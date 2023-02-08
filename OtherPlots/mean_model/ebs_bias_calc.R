###############################################################
## bias calculation for EBS and Lugsail-EBS estimator

###############################################################

## lug = 0 corresponds to EBS
## lug = 1 corresponds to EBS with double batch size
## al  = alpha in the learning rate
## n   = sample size
## bn  = batch size
## an  = number of batches

#########################################################

ebs_bias <- function(n, al, lug = 0)
{
  # al <- 0.51
  bn <- floor(0.1*n^((1+al)/2)) ## batch size
  if(lug == 1) { bn = 2*bn}
  an <- floor(n/bn)
  
  # packages for parallel computation
  require(foreach)
  require(doParallel)
  
  cl <- makeCluster(25) 
  registerDoParallel(cl)
  
  aa <- foreach(j = 1:(an-1), .combine=cbind)%dopar% {
    
    abc = 0
    
    for (k in (j+1):an) 
    {
      for (q in ((k-1)*bn+1):(k*bn)) 
      {
        abc = abc + q^(-al)*sum( (1-q^(-al))^( q -
                    c(((j-1)*bn+1):(j*bn)) ) )
      }
    }
    abc
  }
  stopCluster(cl) # stopping the cluster
  
  return(-sum(aa)/n) ## the bias
}

# values of alpha
alpha = c(0.51, 0.75)

# sample sizes
sam_siz = c(5e4, 1e5, 2e5, 5e5, 8e5, 1e6)

# empty bias matrix
out_ebs <- matrix(NA, nrow = length(sam_siz), ncol = 4)

# computing bias

for (j in 1:length(alpha)) 
{
  al = alpha[j]

for (i in 1:length(sam_siz)) 
{ 
  print(sam_siz[i])     
  
  res1 = ebs_bias(sam_siz[i], al = al, lug = 0) 
  res2 = ebs_bias(sam_siz[i], al = al, lug = 1)
  out_ebs[i,(2*j-1):(2*j)] = c(res1, 2*res2 - res1) 
}
  
}

# saving the output
save(out_ebs, sam_siz, file = "ebs_lugsail_bias.RData")


