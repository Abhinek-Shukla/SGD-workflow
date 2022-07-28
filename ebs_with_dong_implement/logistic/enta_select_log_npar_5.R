rm(list=ls())
set.seed(1)

library(MASS)
library(mcmcse)
source("./../ebs_batch_mean_dong.R")
source("./../sqrt_mat.R")

gradnt_log <- function(y,a,sg){
  
  tmp <- (c( 1 + exp(  y*t(a )%*% (sg))))
  
  p_thet <- 1/ tmp
  
  return(-y*p_thet*a)
}


comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}





  an = 25
  Rep = 2
  nparm = 5
  A = diag(nparm)
  eta_sq <- seq(0.5, 10, by = 0.5)
  sam_siz = c( 1e5, 2e5)#, 5e5, 8e5, 1e6)
  qlev = 0.95
  alp = .51
  burn_in = 10000
  nam_matrix = "indep"
  
  n <- sam_siz[length(sam_siz)] 
  
  parm <- seq(1 / nparm, 1, length.out = nparm)
  load("out/sigm_inv_logist.RData")
  
  sigm_inv <- get(paste("sigm_inv_",nam_matrix,"_nparm_", nparm, sep= ""))
  
  sigm <- qr.solve(sigm_inv)
  
  crt_val <- qchisq(qlev, df = nparm)
  
  
  sg <- matrix(nrow = n + burn_in, ncol = nparm);
  sg_ct <- matrix(nrow = n , ncol = nparm)
  #Iterates stored
  
  # Sigma Matrix Stored with Square root
  
  sqrt_sig <- sqrt_mat(A)
  
  cover_orc <- matrix(rep(0,length(sam_siz)*length(eta_sq)),nrow = length(eta_sq), ncol = length(sam_siz), dimnames = list( 1 : length(eta_sq), sam_siz))
  
  volm_ebs <- forb_ebs <- forb_ebs_norm <- matrix(rep(0,length(sam_siz)*length(eta_sq)),nrow = length(eta_sq), ncol = length(sam_siz), dimnames = list( 1 : length(eta_sq), sam_siz))
  cover_ebs <- volm_ebs_ls <-  matrix(rep(0,length(sam_siz)*length(eta_sq)),nrow = length(eta_sq), ncol = length(sam_siz), dimnames = list( 1 : length(eta_sq), sam_siz))
  forb_ebs_norm_ls <- cover_ebs_ls  <- forb_ebs_ls <- matrix(rep(0,length(sam_siz)*length(eta_sq)),nrow = length(eta_sq), ncol = length(sam_siz), dimnames = list( 1 : length(eta_sq), sam_siz))
  
  library(doParallel) 
  library(foreach)
  
  registerDoParallel(cores = 1)  
  
  final_values <- foreach(cn = 1 : Rep, .combine = "comb", .packages = c("MASS","mcmcse"), .export = c( 'sqrt_mat'), .multicombine=TRUE,
                          .init=list(list(), list(),list(), list(),list(), list(),list(), list(),list() )) %dopar% {
                 
      for (enta in 1 : length(eta_sq)){
        eta_cns <- eta_sq[enta]
    
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
        
        
        for ( smpl in 1 : length(sam_siz)) 
        {
          
          
          sg_ct <- sg_ct_full[1 : sam_siz[smpl], ]
          asg <- colMeans(sg_ct)
          
          
          #Oracle coverage
          cover_orc[enta, smpl] <- as.numeric(sam_siz[smpl]  * t(asg - parm) %*% solve(sigm) %*% (asg - parm) <= crt_val)  
          
          
          crt_val_ebs <- qf(qlev, df1 = nparm, df2 = (an - nparm))
          
          
          bn <- floor(sam_siz[smpl]/an)
          updt_crt <- ( sam_siz[smpl] * nparm / ((an - nparm) * bn)) 
          
          crt_val_ebs <- crt_val_ebs * updt_crt
          
          ebs_mean <- ebs_batch_mean_dong(sg_ct, an, bn,  alp, 1)
          forb_ebs_norm[enta, smpl] <- sqrt(sum((ebs_mean) ^ 2))
        
          forb_ebs[enta, smpl]  <- sqrt(sum((ebs_mean - sigm) ^ 2))/sqrt(sum(sigm ^ 2))
          volm_ebs[enta, smpl]  <- (det(ebs_mean) ) ^ (1 / nparm)
          
          cover_ebs[enta, smpl] <- as.numeric(sam_siz[smpl]  * t(asg - parm ) %*% qr.solve(ebs_mean ) %*% (asg - parm) <= crt_val_ebs)
          
          ebs_mean <- ebs_batch_mean_dong(sg_ct, an, bn,  alp, 2)
          forb_ebs_norm_ls[enta, smpl] <- sqrt(sum((ebs_mean) ^ 2))
          
          forb_ebs_ls[enta, smpl]  <- sqrt(sum((ebs_mean - sigm) ^ 2))/sqrt(sum(sigm ^ 2))
          volm_ebs_ls[enta, smpl]  <- (det(ebs_mean) ) ^ (1 / nparm)
          cover_ebs_ls[enta, smpl] <- as.numeric(sam_siz[smpl]  * t(asg - parm ) %*% qr.solve(ebs_mean ) %*% (asg - parm) <= crt_val_ebs)

         
        }
        
        
      }
                            list(cover_orc, forb_ebs_norm, forb_ebs, volm_ebs,  cover_ebs, forb_ebs_norm_ls, forb_ebs_ls, volm_ebs_ls,   cover_ebs_ls  )    
                          }
  
  
  
  cover_orc <- matrix(rep(0,length(sam_siz)*length(eta_sq)),nrow = length(eta_sq), ncol = length(sam_siz), dimnames = list( 1 : length(eta_sq), sam_siz))
  
  volm_ebs <- forb_ebs <- forb_ebs_norm <- matrix(rep(0,length(sam_siz)*length(eta_sq)),nrow = length(eta_sq), ncol = length(sam_siz), dimnames = list( 1 : length(eta_sq), sam_siz))
  cover_ebs <- volm_ebs_ls <-  matrix(rep(0,length(sam_siz)*length(eta_sq)),nrow = length(eta_sq), ncol = length(sam_siz), dimnames = list( 1 : length(eta_sq), sam_siz))
  forb_ebs_norm_ls <- cover_ebs_ls  <- forb_ebs_ls <- matrix(rep(0,length(sam_siz)*length(eta_sq)),nrow = length(eta_sq), ncol = length(sam_siz), dimnames = list( 1 : length(eta_sq), sam_siz))
  
  for( k in 1 : Rep){
    
    cover_orc <- cover_orc + final_values[[1]][[k]]
    forb_ebs_norm <- forb_ebs_norm + final_values[[2]][[k]]
    forb_ebs <- forb_ebs + final_values[[3]][[k]]
    volm_ebs <- volm_ebs +  final_values[[4]][[k]]
    cover_ebs <- cover_ebs + final_values[[5]][[k]]
    forb_ebs_norm_ls <- forb_ebs_norm_ls + final_values[[6]][[k]]
    forb_ebs_ls <-  forb_ebs_ls +  final_values[[7]][[k]]
    volm_ebs_ls <- volm_ebs_ls +  final_values[[8]][[k]]
    cover_ebs_ls <- cover_ebs_ls +final_values[[9]][[k]]
  
    }
  
  

  
  fil_nam <- paste("out/log_enta_select_dim_",nparm,".RData",sep="")
  save(forb_ebs_norm, forb_ebs_norm_ls,cover_orc,cover_ebs,cover_ebs_ls,volm_ebs,volm_ebs_ls,forb_ebs,forb_ebs_ls,file=fil_nam)
  
  
  
 