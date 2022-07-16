rm(list=ls())
set.seed(1)
library(MASS)
library(doParallel)
library(foreach)
library(mcmcse)


source("./../ebs_batch_mean.R")
source("./../ibs_jasa_mean.R")
source("./../sqrt_mat.R")


gradnt_log <- function(y,a,sg){
  
  tmp <- (c( 1 + exp(  y*t(a )%*% (sg))))
  
  p_thet <- 1/ tmp
  
  return(-y*p_thet*a)
}
##########################################################################################################
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}


Rep = 1
nparm = 20
A = diag( nparm)# Variance covariance matrix of rows of X
cns = seq(0.2, 0.6, by =0.1) # 0.05 for fifty dimension value
ncores_par = 1
eta_cns = 1
sam_siz = 1e5# Performance of EBS will be achieved  4e5 onwards
qlev = 0.95
alp = .505#0.505
burn_in = 5000
niter = 5e6
n <- sam_siz


parm <- seq(1 / nparm, 1, length.out = nparm)

sigm_inv <- matrix(rep(0,nparm^2),nrow = nparm, ncol = nparm)
x1 <- mvrnorm(niter, mu = rep(0, nparm), Sigma = A)

for ( rp in 1:niter){
  tmp2 <- (c((1 + exp(t(x1[rp,]) %*% parm))*(1 + exp(-t(x1[rp,]) %*% parm))*niter))
  sigm_inv <- sigm_inv + x1[rp,] %*% t(x1[rp,]) / tmp2
}

sigm <- qr.solve(sigm_inv)

crt_val <- qchisq(qlev, df = nparm)


sg <- matrix(nrow = n + burn_in, ncol = nparm);
sg_ct <- matrix(nrow = n , ncol = nparm)
#Iterates stored

# Sigma Matrix Stored with Square root

sqrt_sig <- sqrt_mat(A)

cns_ln <- 3*length(cns)# 3 multiplied due to three types of beta under study
forb_ibs <- volm_ibs <- cover_ibs <- forb_ibs_norm <- matrix(rep(0,length(sam_siz)*Rep), nrow = Rep, ncol = length(sam_siz), dimnames = list( 1 : Rep, sam_siz))
cover_orc<- matrix(rep(0,length(sam_siz)*Rep),nrow = Rep, ncol = length(sam_siz), dimnames = list( 1 : Rep, sam_siz))

volm_ebs         <- forb_ebs <- forb_ebs_norm <- array(rep(0, cns_ln * Rep * length(sam_siz)), dim = c(length(sam_siz), Rep, cns_ln), dimnames = list(sam_siz, 1:Rep, 1:cns_ln))
cover_ebs        <- array(rep(0, cns_ln * Rep * length(sam_siz)), dim=c(length(sam_siz), Rep, cns_ln), dimnames = list(sam_siz, 1:Rep, 1:cns_ln))
volm_ebs_ls      <- array(rep(0, cns_ln * Rep * length(sam_siz)), dim=c(length(sam_siz), Rep, cns_ln), dimnames = list(sam_siz, 1:Rep, 1:cns_ln))
forb_ebs_norm_ls <- cover_ebs_ls  <- forb_ebs_ls <- array(rep(0, cns_ln * Rep * length(sam_siz)), dim=c(length(sam_siz), Rep, cns_ln), dimnames = list(sam_siz, 1:Rep, 1:cns_ln))

#Data Generated of Maximum Sample Size
for( cn in 1 : Rep ){           
  
  x <- matrix(rnorm((n + burn_in) * nparm), nrow = (n + burn_in), ncol = nparm)
  x <- x %*% sqrt_sig
  p1 <- 1/(1 + exp(-x %*% parm))# Prob of success for 1 (1-p1 is for -1)
  y<-2*rbinom(n + burn_in, size = 1, prob = p1) - 1
  #Learning Rate
  eta <- numeric(n + burn_in)
  sg[1,] <- rep(0, nparm)
  
  for(i in 2 : (n + burn_in)){
    
    eta[i] <- i^( - alp)
    sg[i,] <- sg[i - 1,] - eta_cns * eta[i] * gradnt_log( y[i], x[i,], sg[i - 1,])
    if(sum(is.na(gradnt_log( y[i], x[i,], sg[i - 1,]))) >= 1){print(c(sg[i-1],i))}
  }
  
  sg_ct_full <- sg[(burn_in + 1) : (n + burn_in),]
  
  
  sg_ct <- sg_ct_full[1:sam_siz, ]
  asg <- colMeans(sg_ct)
  
  
  #Oracle coverage
  cover_orc[cn,1] <- as.numeric(sam_siz  * t(asg - parm) %*% solve(sigm) %*% (asg - parm) <= crt_val)  
  print(as.numeric(sam_siz  * t(asg - parm) %*% solve(sigm) %*% (asg - parm) <= crt_val))
  
  
  #IBS  related coverages and volume
  
  ibs_mean     <- ibs_jasa_mean(sg_ct, alp, cns = 0.05 )# cns = 0.05 or less is running for IBS high dim logistic sample size 1e5
  
  forb_ibs[cn,1]  <-  norm(ibs_mean - sigm, "F")/ norm(sigm, "F")  #sqrt(sum((ibs_mean - sigm) ^ 2))/sqrt(sum(sigm ^ 2))
  forb_ibs_norm[cn,1] <- norm(ibs_mean, "F")
  volm_ibs[cn,1]  <- (det(ibs_mean)) ^ (1 / nparm)
  cover_ibs[cn,1] <- as.numeric(sam_siz  * t(asg - parm) %*% qr.solve(ibs_mean) %*% (asg - parm) <= crt_val)
  
  
  
  count = 1
  #Different settings of EBS, for values of cns and three types of beta
  for( mk in 1 : length(cns)){ #Different values of constant cns
    for(bt_typ in 1 : 3){
      ebs_mean <- ebs_batch_mean(sg_ct, alp, cns[mk], bt_typ, 1)
      forb_ebs_norm[1, cn, count] <- sqrt(sum((ebs_mean) ^ 2))
      
      forb_ebs[1, cn, count]  <- sqrt(sum((ebs_mean - sigm) ^ 2))/sqrt(sum(sigm ^ 2))
      volm_ebs[1, cn, count]  <- (det(ebs_mean) ) ^ (1 / nparm)
      cover_ebs[1, cn, count] <- as.numeric(sam_siz  * t(asg - parm ) %*% qr.solve(ebs_mean ) %*% (asg - parm) <= crt_val)
      
      ebs_mean <- ebs_batch_mean(sg_ct, alp, cns[mk], bt_typ, 2)
      forb_ebs_norm_ls[1, cn, count] <- sqrt(sum((ebs_mean) ^ 2))
      
      forb_ebs_ls[1, cn, count]  <- sqrt(sum((ebs_mean - sigm) ^ 2))/sqrt(sum(sigm ^ 2))
      volm_ebs_ls[1, cn, count]  <- (det(ebs_mean) ) ^ (1 / nparm)
      cover_ebs_ls[1, cn, count] <- as.numeric(sam_siz  * t(asg - parm ) %*% qr.solve(ebs_mean ) %*% (asg - parm) <= crt_val)
      print(sam_siz  * t(asg - parm ) %*% qr.solve(ebs_mean ) %*% (asg - parm))
      count = count + 1           
    }
  }
  
}
#list(forb_ibs[cn,],forb_ibs_norm[cn,], volm_ibs[cn,], cover_ibs[cn,], cover_orc[cn,],  forb_ebs_norm[, cn, ], forb_ebs[, cn, ], volm_ebs[, cn, ],  cover_ebs[, cn, ], forb_ebs_norm_ls[, cn, ], forb_ebs_ls[, cn, ], volm_ebs_ls[, cn, ],  cover_ebs_ls[, cn, ]  )



m_cover_ibs <- mean(cover_ibs) 
m_cover_ebs <- m_cover_ebs_ls <- vector()
for (k in 1: (length(cns)*3))
{
  m_cover_ebs[k] <- mean(cover_ebs[,,k])
}
for (k in 1: (length(cns)*3))
{
  m_cover_ebs_ls[k] <- mean(cover_ebs_ls[,,k])
}
m_cover_orc <- mean(cover_orc)

save(m_cover_ebs, m_cover_ibs, m_cover_orc, m_cover_ebs_ls, file = paste("out/logistic_",n, "sam", Rep, "Reps.RData", sep =))

#fil_nam <- paste("out/linear_", nam_matrix, "_n_",max_sam,"_dim_",nparm,".RData",sep="")
#save(forb_ibs_norm,forb_ebs_norm, forb_ebs_norm_ls,cover_orc,cover_ibs,cover_ebs,cover_ebs_ls,volm_ibs,volm_ebs,volm_ebs_ls,forb_ibs,forb_ebs,forb_ebs_ls,file=fil_nam)