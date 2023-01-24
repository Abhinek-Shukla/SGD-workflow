gradnt_log <- function(y,a,sg){
  
  tmp <- (c( 1 + exp(  y*t(a )%*% (sg))))
  
  p_thet <- 1/ tmp
  
  return(-y*p_thet*a)
}


comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}


log_batch_fn <- function(max_sam = 1e5, Rep = 1, nparm = 5, A = diag(nparm), cns = c(0.1, 1), ncores_par = 1, eta_cns = 0.5, eta_typ = 1,  sam_siz = c(1e5,2e5,5e5,8e5,1e6,5e6,1e7), qlev = 0.95, alp = .505, burn_in = 10000, nam_matrix = "indep"){
  
  sam_siz <- sam_siz[sam_siz <= max_sam]
  n <- sam_siz[length(sam_siz)] 
  
  parm <- seq(1 / nparm, 1, length.out = nparm)
  load("out/sigm_inv_logist.RData")
  
  sigm_inv <- get(paste("sigm_inv_",nam_matrix,"_nparm_", nparm, sep= ""))
  
  sigm <- qr.solve(sigm_inv)
  
  
  
  
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
  
  
  
  library(doParallel) 
  library(foreach)
  
  registerDoParallel(cores = min(Rep, ncores_par))  
  
  
  
  
  #1000  Replications to obtain stable results under Parallel Setting
  final_values <- foreach(cn = 1 : Rep, .combine = "comb", .packages = c("MASS","mcmcse"), .export = c( 'ebs_batch_mean', 'ibs_jasa_mean', 'sqrt_mat', 'gradnt_log'), .multicombine=TRUE,
                          .init=list(list(), list(),list(), list() )) %dopar% {
                            crt_val <- qnorm((1+qlev)/2)
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
                              #print(asg)
                              ons <- rep(1,nparm)
                              
                              #Oracle coverage
                              print(sam_siz[smpl]) 
                              tmp_sd <- sqrt((t(ons) %*% (sigm) %*% ons)/sam_siz[smpl])
                              low_orc <- t(ons) %*% asg - crt_val * tmp_sd 
                              up_orc <-  t(ons) %*% asg + crt_val * tmp_sd 
                              cover_orc[cn, smpl] <- as.numeric (t(ons) %*% parm <= up_orc & t(ons) %*% parm >= low_orc)
                             print("oracle")
                               print(c(t(ons) %*% asg, t(ons) %*% parm, crt_val * tmp_sd ))
                              
                              #IBS  related coverages and volume
                              
                              ibs_mean     <- ibs_jasa_mean(sg_ct, alp, cns = 0.01)
                              
                           
                              tmp_sd <- sqrt((t(ons) %*% ibs_mean%*%ons)/sam_siz[smpl])
                              low_ibs <- t(ons) %*% asg - crt_val * tmp_sd
                              up_ibs <- t(ons) %*% asg + crt_val * tmp_sd
                              cover_ibs[cn, smpl] <- as.numeric (t(ons)%*%parm<=up_ibs&t(ons)%*%parm>=low_ibs)
                              print("ibs")
                              print(c(t(ons) %*% asg, t(ons) %*% parm, crt_val * tmp_sd ))
                              
                              
                              count = 1
                              #Different settings of EBS, for values of cns and three types of beta
                              for( mk in 1 : length(cns)){ #Different values of constant cns
                                for(bt_typ in 3 : 3){
                                  
                                  ebs_mean <- ebs_batch_mean(sg_ct, alp, cns[mk], bt_typ, 1)
                                    
                                     tmp_sd <- sqrt((t(ons) %*% ebs_mean%*%ons)/sam_siz[smpl])
                                     low_ebs <- t(ons) %*% asg - crt_val * tmp_sd
                                     up_ebs <- t(ons) %*% asg + crt_val * tmp_sd
                                     cover_ebs[smpl,  cn, count] <- as.numeric (t(ons)%*%parm<=up_ebs&t(ons)%*%parm>=low_ebs)
                                     print("EBS")
                                     print(c(t(ons) %*% asg, t(ons) %*% parm, crt_val * tmp_sd ))
                                     
                                     
                                     
                                  ebs_mean <- ebs_batch_mean(sg_ct, alp, cns[mk], bt_typ, 2)
                                  tmp_sd <- sqrt((t(ons) %*% ebs_mean%*%ons)/sam_siz[smpl])
                                  low_ebs <- t(ons) %*% asg - crt_val * tmp_sd
                                  up_ebs <- t(ons) %*% asg + crt_val * tmp_sd
                                  cover_ebs_ls[smpl,  cn, count] <- as.numeric (t(ons)%*%parm<=up_ebs&t(ons)%*%parm>=low_ebs)
                                  print("Lugsail")
                                  print(c(t(ons) %*% asg, t(ons) %*% parm, crt_val * tmp_sd ))
                                  
                                  
                                  
                                  count = count + 1     
                                  
                                  
                                  
                                  
                                }
                              } 
                            }
                            
                            
                            list(  cover_orc[cn,], cover_ibs[cn,],   cover_ebs[, cn, ],  cover_ebs_ls[, cn, ]  )
                            
                          }
  
  for( k in 1 : Rep){
    
   
    cover_orc[k,] <- final_values[[1]][[k]]
    cover_ibs[k,] <- final_values[[2]][[k]]
    cover_ebs[, k, ] <- final_values[[3]][[k]]
    cover_ebs_ls[, k, ] <- final_values[[4]][[k]]
  }
  
  
  
  fil_nam <- paste("out/logistic_proj_new_", nam_matrix, "_n_",max_sam,"_dim_",nparm , ".RData",sep="")
  save(cover_orc,cover_ibs,cover_ebs,cover_ebs_ls,file=fil_nam)
}