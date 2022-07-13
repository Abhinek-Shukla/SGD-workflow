rm(list=ls())
set.seed(1)
library(MASS)
library(doParallel)
library(foreach)
library(mcmcse)

source("grad_lin_and_batch.R")
source("./../ebs_batch_mean.R")
source("./../ibs_jasa_mean.R")
source("./../sqrt_mat.R")


grad_lin <- function(sg,y,x){
  (x %*% (sg) - y) %*% x
} 

##########################################################################################################
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}


  Rep = 1
  nparm = 20
  A = diag(nparm)
  cns = c(0.1)
  ncores_par = 1
  eta_cns = 0.5
  sam_siz = 1e5
  qlev = 0.95
  alp = .51
  burn_in = 1000
sigm <- qr.solve(A) 

n <- sam_siz


parm <- seq(1 / nparm, 1, length.out = nparm)
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
                          
                          cn <- 1
                          x <- matrix(rnorm((n + burn_in) * nparm), nrow = (n + burn_in), ncol = nparm)
                          x <- x %*% sqrt_sig
                          #noisy Observed Data
                          y <- x %*% parm + rnorm((n + burn_in), mean = 0, sd = 1)
                          #Learning Rate
                          eta <- numeric(n + burn_in)
                          sg[1,] <- rep(0, nparm)
                          
                          for(i in 2 : (n + burn_in)){
                            
                            eta[i] <- i^( - alp)
                            sg[i,] <- sg[i - 1,] - eta_cns * eta[i] * grad_lin(sg[i - 1,], y[i], x[i,])
                            
                          }
                          
                          sg_ct_full <- sg[(burn_in + 1) : (n + burn_in),]
                          
                          for ( smpl in 1 : 1) 
                          { 
                            sg_ct <- sg_ct_full[1:sam_siz[smpl], ]
                            asg <- colMeans(sg_ct)
                            
                            #IBS and Oracle related coverages and volume
                            
                            ibs_mean     <- ibs_jasa_mean(sg_ct, alp)
                            #print((ibs_mean))
                            forb_ibs[cn,smpl]  <-  norm(ibs_mean - sigm, "F")/ norm(sigm, "F")  #sqrt(sum((ibs_mean - sigm) ^ 2))/sqrt(sum(sigm ^ 2))
                            forb_ibs_norm[cn,smpl] <- norm(ibs_mean, "F")
                            volm_ibs[cn,smpl]  <- (det(ibs_mean)) ^ (1 / nparm)
                            cover_ibs[cn,smpl] <- as.numeric(sam_siz[smpl]  * t(asg - parm) %*% qr.solve(ibs_mean) %*% (asg - parm) <= crt_val)
                            
                            
                            cover_orc[cn,smpl] <- as.numeric(sam_siz[smpl]  * t(asg - parm) %*% solve(sigm) %*% (asg - parm) <= crt_val)  
                            
                            count = 1
                            #Different settings of EBS, for values of cns and three types of beta
                            for( mk in 1 : length(cns)){ #Different values of constant cns
                              for(bt_typ in 1 : 3){
                                ebs_mean <- ebs_batch_mean(sg_ct, alp, cns[mk], bt_typ, 1)
                                forb_ebs_norm[smpl, cn, count] <- sqrt(sum((ebs_mean) ^ 2))
                                
                                forb_ebs[smpl, cn, count]  <- sqrt(sum((ebs_mean - sigm) ^ 2))/sqrt(sum(sigm ^ 2))
                                volm_ebs[smpl, cn, count]  <- (det(ebs_mean) ) ^ (1 / nparm)
                                cover_ebs[smpl, cn, count] <- as.numeric(sam_siz[smpl]  * t(asg - parm ) %*% qr.solve(ebs_mean ) %*% (asg - parm) <= crt_val)
                                
                                ebs_mean <- ebs_batch_mean(sg_ct, alp, cns[mk], bt_typ, 2)
                                forb_ebs_norm_ls[smpl, cn, count] <- sqrt(sum((ebs_mean) ^ 2))
                                
                                forb_ebs_ls[smpl, cn, count]  <- sqrt(sum((ebs_mean - sigm) ^ 2))/sqrt(sum(sigm ^ 2))
                                volm_ebs_ls[smpl, cn, count]  <- (det(ebs_mean) ) ^ (1 / nparm)
                                cover_ebs_ls[smpl, cn, count] <- as.numeric(sam_siz[smpl]  * t(asg - parm ) %*% qr.solve(ebs_mean ) %*% (asg - parm) <= crt_val)
                                
                                count = count + 1           
                              }
                            }
                            
                          }
                          #list(forb_ibs[cn,],forb_ibs_norm[cn,], volm_ibs[cn,], cover_ibs[cn,], cover_orc[cn,],  forb_ebs_norm[, cn, ], forb_ebs[, cn, ], volm_ebs[, cn, ],  cover_ebs[, cn, ], forb_ebs_norm_ls[, cn, ], forb_ebs_ls[, cn, ], volm_ebs_ls[, cn, ],  cover_ebs_ls[, cn, ]  )
                       

print(cover_orc)
print(cover_ibs)
print(cover_ebs)
print(cover_ebs_ls)
#fil_nam <- paste("out/linear_", nam_matrix, "_n_",max_sam,"_dim_",nparm,".RData",sep="")
#save(forb_ibs_norm,forb_ebs_norm, forb_ebs_norm_ls,cover_orc,cover_ibs,cover_ebs,cover_ebs_ls,volm_ibs,volm_ebs,volm_ebs_ls,forb_ibs,forb_ebs,forb_ebs_ls,file=fil_nam)