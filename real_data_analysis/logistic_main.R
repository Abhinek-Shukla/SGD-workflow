##########################################################
## The file contains the functions required to run
## comparative inference for SGD for different estimators
##########################################################

### Gradient Function for Logistic Model ###
gradnt_log <- function(y, a, sg)
{
  
  tmp <- as.numeric(( 1 + exp(  y * sum(a * sg))))
  
  p_thet <- 1/ tmp
  
  return(-y * p_thet)
}


###############################################################
## Function to obtain performance of estimators in Linear Reg
## max_sam = the maximum number of data size
## Rep = number of replications for the simulation
## nparm = dimension of the regression problem (without intercept models only)
## A = covariance matrix for generating X
## cns = multiplicative constants for batch size calculations
## eta_cns = multicative constant in learning rate
## sam_siz = the sample sizes at which performance will be checked
## qlev = confidence level
## alp = alpha in the learning rate
## burn_in = how much of the beginning to discard
## cns1 = multiplicative constants for batch size calculation for unequal batch size
###############################################################

log_batch_fn <- function(max_sam = 2e5, eta_cns = 0.05, alp = .51, burn_in = 5000, sam_siz = c(1e5,2e5), dta, qlev = 0.95, cns = c(0.1), cns1 = 0.01)
{
nparm = (ncol(dta) - 1)
max_sam = nrow(dta)
sam_siz <- sam_siz[sam_siz <= max_sam]    # grid of sample sizes we will check at
init <- 10000
n <- sam_siz[length(sam_siz)] - burn_in         # Largest sample size

sg    <- matrix(nrow = n + burn_in - init, ncol = nparm) 
sg_ct <- matrix(nrow = n - init, ncol = nparm) # SGD iterates after dropping burn-in samples

cns_ln <- 3*length(cns)    # 3 types of beta 


# For IBS estimator
volm_ibs             <- numeric(length(sam_siz))
marg_volm_ibs        <- numeric(length(sam_siz))
marg_sim_cov_ibs     <- numeric(length(sam_siz))
forb_ibs             <- numeric(length(sam_siz))

# For EBS estimator  
volm_ebs             <- matrix(rep(0, cns_ln * length(sam_siz)), nrow = length(sam_siz), ncol = cns_ln, dimnames = list(sam_siz, 1 : cns_ln))
marg_volm_ebs        <- matrix(rep(0, cns_ln * length(sam_siz)), nrow = length(sam_siz), ncol = cns_ln, dimnames = list(sam_siz, 1 : cns_ln))  
marg_sim_cov_ebs     <- matrix(rep(0, cns_ln * length(sam_siz)), nrow = length(sam_siz), ncol = cns_ln, dimnames = list(sam_siz, 1 : cns_ln))
forb_ebs             <- matrix(rep(0, cns_ln * length(sam_siz)), nrow = length(sam_siz), ncol = cns_ln, dimnames = list(sam_siz, 1 : cns_ln))

# For EBS with Lugsail
volm_ebs_ls          <- matrix(rep(0, cns_ln * length(sam_siz)), nrow = length(sam_siz), ncol = cns_ln, dimnames = list(sam_siz, 1 : cns_ln))
marg_volm_ebs_ls     <- matrix(rep(0, cns_ln * length(sam_siz)), nrow = length(sam_siz), ncol = cns_ln, dimnames = list(sam_siz, 1 : cns_ln))
marg_sim_cov_ebs_ls  <- matrix(rep(0, cns_ln * length(sam_siz)), nrow = length(sam_siz), ncol = cns_ln, dimnames = list(sam_siz, 1 : cns_ln))
forb_ebs_ls          <- matrix(rep(0, cns_ln * length(sam_siz)), nrow = length(sam_siz), ncol = cns_ln, dimnames = list(sam_siz, 1 : cns_ln))

# For relative volume comparisons  
ratio_ibs_ebs        <- array(rep(0, nparm * cns_ln * length(sam_siz)), dim = c(length(sam_siz), cns_ln, nparm), dimnames = list(sam_siz, 1 : cns_ln, 1 : nparm))
ratio_ibs_ebs_ls     <- array(rep(0, nparm * cns_ln * length(sam_siz)), dim = c(length(sam_siz), cns_ln, nparm), dimnames = list(sam_siz, 1 : cns_ln, 1 : nparm))
ratio_ebs_ls_ebs     <- array(rep(0, nparm * cns_ln * length(sam_siz)), dim = c(length(sam_siz), cns_ln, nparm), dimnames = list(sam_siz, 1 : cns_ln, 1 : nparm))

##########################################
### Dependent and Independent Variable ###
##########################################
y <- as.vector(dta[1 : (n + burn_in), 1])
x <- as.matrix(dta[1 : (n + burn_in), -1])

#############################
### Getting initial guess ###
#############################
y1 <- (y[1 : init] + 1) / 2 ### Conversion for applying GLM ###
x1 <- x[ 1 : init, ]
dt_fr <- data.frame(y1,x1)
model_train <- glm(formula = y1 ~ x1 + 0, family = binomial, data = dt_fr)
sg[1,] <- as.vector(model_train$coefficients) #MLE


####################
### SGD iterates ###
####################
eta <- numeric(n + burn_in - init)
for(i in 2 : (n + burn_in - init)){
  
  eta[i] <- i^( - alp)
  sg[i,] <- sg[(i - 1),] - eta_cns * eta[i] * gradnt_log( y[(init + i)], x[(init + i),], sg[(i - 1),]) * x[(init + i),]
  
}

####################################
### Maximum Likeliood Estimators ###
####################################
y1 <- (y + 1) / 2 
dt_fr <- data.frame(y1,x)
model_train <- glm(formula = y1 ~ x + 0, family = binomial, data = dt_fr)
mle_s <- as.vector(model_train$coefficients) #MLE

sg_ct_full <- sg[(burn_in + 1) : (n + burn_in - init),]
crt_val <- qchisq(qlev, df = nparm)

for(smpl in 1 : length(sam_siz)) 
{ 
  sg_ct <- sg_ct_full[1 : (sam_siz[smpl] - init - burn_in) , ]  # main process inside this loop
  asg <- colMeans(sg_ct)                    # Averaged SGD
  
  # constant next to the volume of an ellipsoidal confidence region
  tmp_vol <- 2*(pi^(nparm/2)) / (nparm * gamma(nparm/2)) * (crt_val/sam_siz[smpl]) ^ (nparm/2) 
  
  margn_up_low <- list()
  
  ########################
  ### IBS calculations ###
  ########################
  ibs_mean <- ibs_jasa_mean(sg_ct, alp, cns = cns1)
  volm_ibs[smpl]  <- (tmp_vol * (det(ibs_mean)) ^ (1 / 2)) ^ (1 / nparm)
  tmp_ibs <- new.sim.int(ibs_mean/sam_siz[smpl], conf = 0.95, center = asg)$ints   # Marginal Simultaneous for IBS
  leng_ibs <- tmp_ibs[, 2] - tmp_ibs[, 1]
  marg_volm_ibs[smpl] <- (prod(leng_ibs) ^(1 / nparm) / volm_ibs[smpl])
  margn_up_low[[1]] <- tmp_ibs
  # Different settings of equal batch size estimator (EBS), for values of cns and three types of beta
  count <- 1
  cnt <- 2
  for(mk in 1 : length(cns))
  { 
    for(bt_typ in 1 : 3)
    { 
    #######################
    ### EBS calculation ###
    #######################
    ebs_mean <- ebs_batch_mean(sg_ct, alp, cns[mk], bt_typ, 1)
    volm_ebs[smpl, count]  <- (tmp_vol * (det(ebs_mean) ) ^ (1 / 2))^( 1/ nparm)
    tmp_ebs <- new.sim.int(ebs_mean/sam_siz[smpl], conf = 0.95, center = asg)$ints
    leng_ebs <- tmp_ebs[, 2] - tmp_ebs[, 1]
    ratio_ibs_ebs[smpl, count, ] <- leng_ibs/leng_ebs
    marg_volm_ebs[smpl, count] <- (prod(leng_ebs) ^ (1 / nparm)  / volm_ebs[smpl, count])
    margn_up_low[[cnt]] <- tmp_ebs 
    cnt <- cnt + 1
    ###############################
    ### EBS Lugsail calculation ###
    ###############################
    ebs_mean <- ebs_batch_mean(sg_ct, alp, cns[mk], bt_typ, 2)
    volm_ebs_ls[smpl, count]  <- (tmp_vol * (det(ebs_mean) ) ^ (1 / 2)) ^ (1 / nparm)
    tmp_ebs_ls <- new.sim.int(ebs_mean/sam_siz[smpl], conf = 0.95, center = asg)$ints
    leng_ebs_ls <- tmp_ebs_ls[, 2] - tmp_ebs_ls[, 1]
    marg_volm_ebs_ls[smpl, count] <- (prod(leng_ebs_ls) ^ (1 / nparm) / volm_ebs_ls[smpl, count])         
    ratio_ibs_ebs_ls[smpl, count, ] <- leng_ibs/leng_ebs_ls
    margn_up_low[[cnt]] <- tmp_ebs_ls 
    cnt <- cnt + 1
    ######################################
    ### EBS Lugsail v/s without Lugsail ### 
    ######################################
    ratio_ebs_ls_ebs[smpl, count, ] <- leng_ebs_ls/leng_ebs
      
    count = count + 1 
    
    }
  }  
}


fil_nam <- paste("out/logistic_real","_dim_", nparm, ".RData", sep = "")
save(n, margn_up_low, init, burn_in,  marg_volm_ibs, marg_volm_ebs, marg_volm_ebs_ls, ratio_ibs_ebs, ratio_ibs_ebs_ls, ratio_ebs_ls_ebs,  volm_ibs, volm_ebs, volm_ebs_ls, file = fil_nam)
}

