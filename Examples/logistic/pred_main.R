# Functions required to run comparative inference for SGD


###############################################################

## max_sam = the maximum number of data size
## nparm   = dimension of the regression problem 
## cns     = multiplicative constants for batch size calculations EBS
## eta_cns = multicative constant in learning rate
## qlev    = confidence level
## alp     = alpha in the learning rate
## cns1    = multiplicative constants for batch size calculation for IBS

###############################################################

log_batch_fn <- function(max_sam = 2e5, eta_cns = 0.05, alp = .51, 
                         burn_in = 5000, sam_siz = c(1e5,2e5), dta,
                         qlev = 0.95, cns = 0.1, cns1 = 0.01)
{
nparm = (ncol(dta) - 1)
max_sam = nrow(dta)  
init <- 10000
n <- max_sam - burn_in - init

sg    <- matrix(nrow = n + burn_in - init, ncol = nparm) 
sg_ct <- matrix(nrow = n - init, ncol = nparm) 

cns_ln <- 3*length(cns)    



#Maximum Likelihood Estimators

y1          <- (y + 1) / 2 
dt_fr       <- data.frame(y1,x)
# model_train <- glm(formula = y1 ~ x + 0, family = binomial, data = dt_fr)
# mle_s       <- as.vector(model_train$coefficients)

sg_ct_full  <- logistic_sgd(n = n, burn_in = burn_in, 
              dta = dta, init = init, alp = alp, eta_cns = eta_cns)

asg   <- colMeans(sg_ct_full)

ebs <- ebs_batch_mean(sg_ct_full, alp, cns, 1, 1)
lug.ebs <- ebs_batch_mean(sg_ct_full, alp, cns, 1, 2)

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

cutoff <- .2
y.hat <- (pi.hat > cutoff)
y.hat.lb <- (lb > cutoff)
ind <- which(y.test == 0)
(sum( (y.hat[ind] - y.test[ind]) != 0)/length(ind))
(sum( (y.hat.lb[ind] - y.test[ind]) != 0)/length(ind))


fil_nam <- paste("out/logistic_real","_dim_", nparm, ".RData", sep = "")
save(n, margn_up_low, init, burn_in,  marg_volm_ibs, marg_volm_ebs, 
     marg_volm_ebs_ls, ratio_ibs_ebs, ratio_ibs_ebs_ls, ratio_ebs_ls_ebs, 
     volm_ibs, volm_ebs, volm_ebs_ls, file = fil_nam)
}


test <- read.csv("testing.csv", header = TRUE)
test <- test[,-1]
y.test <- test[ ,1]
x.test <- as.matrix(test[ ,-1])


