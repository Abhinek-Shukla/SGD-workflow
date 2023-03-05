########################################
## Function that returns probability
## of success for test data
###############################################################

## dta     = dataset
## cns     = multiplicative constants for batch size calculations EBS
## eta_cns = multicative constant in learning rate
## qlev    = confidence level
## alp     = alpha in the learning rate

###############################################################

log_batch_fn <- function(dta, test, eta_cns = .05, alp = .51, 
                         burn_in = 5000, qlev = 0.95, cns = 0.1)
{
  nparm = (ncol(dta) - 1)
  max_sam = nrow(dta)  
  init <- 10000
  n <- max_sam - burn_in - init

  sg    <- matrix(nrow = n + burn_in - init, ncol = nparm) 
  sg_ct <- matrix(nrow = n - init, ncol = nparm) 

  sg_ct_full  <- logistic_sgd(n = n, burn_in = burn_in, 
                dta = dta, init = init, alp = alp, eta_cns = eta_cns,
                epochs = 1)

  asg   <- colMeans(sg_ct_full)

  ebs <- ebs_batch_mean(sg_ct_full, alp, cns, 1, 1)

  y.test <- test[ ,1]
  x.test <- as.matrix(test[ ,-1])

  n.test <- length(y.test)
  pi.hat <- numeric(length = n.test)
  lb <- numeric(length = n.test)
  lb.lug <- numeric(length = n.test)

  for(i in 1:n.test)
  {
    xi <- x.test[i, ]
    foo <- exp( sum(xi*asg))
    pi.hat[i] <- foo/(1 + foo)

    st.err <- sqrt((pi.hat[i] * (1 - pi.hat[i]) )^2 * xi %*% ebs %*% xi/n)
    st.err.lug <- sqrt((pi.hat[i] * (1 - pi.hat[i]) )^2 * xi %*% ebs %*% xi/n)
    lb[i] <- pi.hat[i] - qnorm((1 - (1-qlev)/2)) * st.err
    lb.lug[i] <- pi.hat[i] - qnorm((1 - (1-qlev)/2)) * st.err.lug
  }

  cutoffs <- seq(0.01, .3, length = 50)

  misclass <- numeric(length = length(cutoffs))
  misclass.lb <- numeric(length = length(cutoffs))
  for(p in 1:length(cutoffs))
  {
    y.hat <- (pi.hat > cutoffs[p])
    y.hat.lb <- (lb > cutoffs[p])

    misclass[p] <- sum(y.test != y.hat)/length(y.hat)
    misclass.lb[p] <- sum(y.test != y.hat.lb) /length(y.hat)
  }
  save(cutoffs, misclass, misclass.lb, file = "out/logistic_pred.Rdata")
}




