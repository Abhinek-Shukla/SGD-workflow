# Functions required to run comparative inference for SGD

# Gradient Function for LAD Model
grad_lad <- function(sg,y,x, tau = 0.5)
{
  (tau - as.numeric( y - t(x) %*% sg < 0 )) %*% x
} 

# Combining function to get outputs from  parallel
comb <- function(x, ...) 
{
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}


###############################################################

## A          = covariance matrix for generating X
## cns        = multiplicative constants for batch size calculations EBS
## eta_cns    = multiplicative constant in learning rate
## sam_siz    = the sample sizes at which performance will be checked
## qlev       = confidence level
## alp        = alpha in the learning rate
## nam_matrix = string descriptor for A
## cns1       = multiplicative constants for batch size calculation IBS

###############################################################

lad_reps <- function(max_sam = 1e5, Rep = 1, nparm = 5, A = diag(nparm), 
                     cns = c(0.1, 1), ncores_par = 1, eta_cns = 0.5,
                     sam_siz = c(5e4, 1e5, 2e5, 5e5, 8e5, 1e6, 5e6, 1e7), 
                     qlev = 0.95,  alp = .51, burn_in = 1000, 
                     nam_matrix = "indep", cns1 = 1, tau = 0.5)
{
 
  # Oracle Sigma Matrix for Lad Model
  sigm    <- qr.solve(A) * (tau*(1-tau))/ (1/2^2)
 
  sam_siz <- sam_siz[sam_siz <= max_sam]    
  n       <- sam_siz[length(sam_siz)]          
  
  parm    <- seq(1 / nparm, 1, length.out = nparm)  
  crt_val <- qchisq(qlev, df = nparm)
  
  
  sg    <- matrix(nrow = n + burn_in, ncol = nparm) 
  sg_ct <- matrix(nrow = n, ncol = nparm) 
  
  
  sqrt_sig <- sqrt_mat(A)

  
  cns_ln <- 3*length(cns)    
  # For oracle
  volm_orc         <- matrix(rep(0, length(sam_siz)*Rep), 
                                 nrow = Rep, ncol = length(sam_siz), 
                                 dimnames = list( 1 : Rep, sam_siz))
  marg_volm_orc    <- matrix(rep(0, length(sam_siz)*Rep), 
                                 nrow = Rep, ncol = length(sam_siz), 
                                 dimnames = list( 1 : Rep, sam_siz))
  cover_orc        <- matrix(rep(0, length(sam_siz)*Rep),
                                 nrow = Rep, ncol = length(sam_siz),
                                 dimnames = list( 1 : Rep, sam_siz))
  marg_sim_cov_orc <- matrix(rep(0, length(sam_siz)*Rep), 
                                 nrow = Rep, ncol = length(sam_siz),
                                 dimnames = list( 1 : Rep, sam_siz))

  # For IBS estimator
  volm_ibs             <- matrix(rep(0, length(sam_siz)*Rep),
                                 nrow = Rep, ncol = length(sam_siz), 
                                 dimnames = list( 1 : Rep, sam_siz))
  marg_volm_ibs        <- matrix(rep(0, length(sam_siz)*Rep),
                                 nrow = Rep, ncol = length(sam_siz),
                                 dimnames = list( 1 : Rep, sam_siz))
  cover_ibs            <- matrix(rep(0, length(sam_siz)*Rep), 
                                 nrow = Rep, ncol = length(sam_siz),
                                 dimnames = list( 1 : Rep, sam_siz))
  marg_sim_cov_ibs     <- matrix(rep(0, length(sam_siz)*Rep), 
                                 nrow = Rep, ncol = length(sam_siz),
                                 dimnames = list( 1 : Rep, sam_siz))
  forb_ibs             <- matrix(rep(0, length(sam_siz)*Rep),
                                 nrow = Rep, ncol = length(sam_siz),
                                 dimnames = list( 1 : Rep, sam_siz))
  
  # For EBS estimator
  volm_ebs         <- array(rep(0, cns_ln * Rep * length(sam_siz)), 
                            dim = c(length(sam_siz), Rep, cns_ln), 
                            dimnames = list(sam_siz, 1 : Rep, 1 : cns_ln))
  marg_volm_ebs    <- array(rep(0, cns_ln * Rep * length(sam_siz)),
                            dim = c(length(sam_siz), Rep, cns_ln), 
                            dimnames = list(sam_siz, 1 : Rep, 1 : cns_ln))  
  cover_ebs        <- array(rep(0, cns_ln * Rep * length(sam_siz)),
                            dim = c(length(sam_siz), Rep, cns_ln), 
                            dimnames = list(sam_siz, 1 : Rep, 1 : cns_ln)) 
  marg_sim_cov_ebs <- array(rep(0, cns_ln * Rep * length(sam_siz)), 
                            dim = c(length(sam_siz), Rep, cns_ln),
                            dimnames = list(sam_siz, 1 : Rep, 1 : cns_ln))
  forb_ebs         <- array(rep(0, cns_ln * Rep * length(sam_siz)), 
                            dim = c(length(sam_siz), Rep, cns_ln),
                            dimnames = list(sam_siz, 1 : Rep, 1 : cns_ln))

  # For EBS with Lugsail
  volm_ebs_ls          <- array(rep(0, cns_ln * Rep * length(sam_siz)), 
                                dim = c(length(sam_siz), Rep, cns_ln), 
                                dimnames = list(sam_siz, 1 : Rep, 1 : cns_ln))
  marg_volm_ebs_ls     <- array(rep(0, cns_ln * Rep * length(sam_siz)), 
                                dim = c(length(sam_siz), Rep, cns_ln), 
                                dimnames = list(sam_siz, 1 : Rep, 1 : cns_ln))
  cover_ebs_ls         <- array(rep(0, cns_ln * Rep * length(sam_siz)),
                                dim = c(length(sam_siz), Rep, cns_ln), 
                                dimnames = list(sam_siz, 1 : Rep, 1 : cns_ln))
  marg_sim_cov_ebs_ls  <- array(rep(0, cns_ln * Rep * length(sam_siz)), 
                                dim = c(length(sam_siz), Rep, cns_ln), 
                                dimnames = list(sam_siz, 1 : Rep, 1 : cns_ln))
  forb_ebs_ls          <- array(rep(0, cns_ln * Rep * length(sam_siz)), 
                                dim = c(length(sam_siz), Rep, cns_ln), 
                                dimnames = list(sam_siz, 1 : Rep, 1 : cns_ln))
 
  # For relative volume comparisons
  ratio_ibs_ebs   <- array(rep(0, nparm * cns_ln * Rep * length(sam_siz)), 
                          dim = c(length(sam_siz), Rep, cns_ln, nparm),
                          dimnames = list(sam_siz, 1 : Rep, 1 : cns_ln, 1 : nparm))
  ratio_ibs_ebs_ls <- array(rep(0, nparm * cns_ln * Rep * length(sam_siz)),
                            dim = c(length(sam_siz), Rep, cns_ln, nparm), 
                            dimnames = list(sam_siz, 1 : Rep, 1 : cns_ln, 1 : nparm))
  ratio_ebs_ls_ebs <- array(rep(0, nparm * cns_ln * Rep * length(sam_siz)),
                            dim = c(length(sam_siz), Rep, cns_ln, nparm), 
                            dimnames = list(sam_siz, 1 : Rep, 1 : cns_ln, 1 : nparm))

  
  registerDoParallel(cores = min(Rep, ncores_par))  
  
  final_values <- foreach(cn = 1:Rep, .combine = "comb",
                          .packages = c("MASS","mcmcse", "mvtnorm", "smoothmest"),
                          .export = c( 'ebs_batch_mean', 'ibs_jasa_mean', 'sqrt_mat',
                                       'grad_lad', 'new.sim.int', 'opt_beta_fn', 
                                       'b_our'), .multicombine=TRUE,
                            .init = list(list(), list(), list(), list(),
                                         list(), list(), list(), list(), list(),
                                         list(), list(), list(), list(), list(), 
                                         list(), list(), list(), list(), list(), 
                                         list(), list(), list() )) %dopar% 
  {
    
    # Generating data of Maximum Sample Size
    x <- matrix(rnorm((n + burn_in) * nparm),
                nrow = (n + burn_in), ncol = nparm)
    x <- x %*% sqrt_sig
    y <- x %*% parm + rdoublex(n + burn_in, mu = 0,lambda = 1)
    
    # Learning Rate
    eta <- numeric(n + burn_in)
    
    sg[1, ] <- rep(0, nparm)  

    ## Generating full SGD sequence -- this takes time
    for(i in 2 : (n + burn_in))
    {
      eta[i] <- i^( - alp)
      sg[i, ] <- sg[i - 1, ] + eta_cns * 
        eta[i] * grad_lad(sg[i - 1, ], y[i], x[i, ]) 
    }
    sg_ct_full <- sg[(burn_in + 1) : (n + burn_in), ]
    


    for(smpl in 1:length(sam_siz)) 
    { 
      sg_ct <- sg_ct_full[1 : sam_siz[smpl], ]  
      asg <- colMeans(sg_ct)                    
      
      # constant next to the volume of an ellipsoidal confidence region
      tmp_vol <- 2*(pi^(nparm/2)) / (nparm * gamma(nparm/2)) 
      * (crt_val/sam_siz[smpl]) ^ (nparm/2) 

      # Oracle calculations
      
      volm_orc[cn, smpl]         <- tmp_vol * (det(sigm)) ^ (1 / 2)
      cover_orc[cn, smpl]        <- as.numeric(sam_siz[smpl] 
                                    * t(asg - parm) %*% qr.solve(sigm) %*% 
                                    (asg - parm) <= crt_val)  
      
      tmp_orc                    <- new.sim.int(sigm/sam_siz[smpl], conf = 0.95, 
                                    center = asg)$ints
      leng_orc                   <- tmp_orc[, 2] - tmp_orc[, 1]      
      marg_volm_orc[cn, smpl]    <- (prod(leng_orc) 
                                  / volm_orc[cn, smpl])^(1 / nparm)      
      marg_sim_cov_orc[cn, smpl] <- as.numeric(sum(tmp_orc[, 1] <= parm &
                                       tmp_orc[, 2] >= parm) == nparm)
      # IBS calculations
     
      ibs_mean     <- ibs_jasa_mean(sg_ct, alp, cns = cns1)

      forb_ibs[cn, smpl]  <-  norm(ibs_mean - sigm, "F")/ norm(sigm, "F")  
      volm_ibs[cn, smpl]  <- tmp_vol * (det(ibs_mean)) ^ (1 / 2)
      cover_ibs[cn, smpl] <- as.numeric(sam_siz[smpl]  * t(asg - parm) 
                              %*% qr.solve(ibs_mean) %*%
                                (asg - parm) <= crt_val)
      
      tmp_ibs  <- new.sim.int(ibs_mean/sam_siz[smpl], conf = 0.95,
                             center = asg)$ints   
      leng_ibs <- tmp_ibs[, 2] - tmp_ibs[, 1]
      
      
      marg_sim_cov_ibs[cn, smpl] <- as.numeric(sum(tmp_ibs[, 1] 
                                    <= parm & tmp_ibs[, 2] >= 
                                                     parm) == nparm)
      marg_volm_ibs[cn, smpl]    <- (prod(leng_ibs) / 
                                    volm_ibs[cn, smpl])^(1 / nparm)
      count <- 1
     
      for(mk in 1 : length(cns))
      { 
        for(bt_typ in 1 : 3)
        { 
         # EBS calculation 
          
          ebs_mean <- ebs_batch_mean(sg_ct, alp, cns[mk], bt_typ, 1)
          
          forb_ebs[smpl, cn, count]  <- norm(ebs_mean - sigm, "F")
                                        / norm(sigm, "F")
          volm_ebs[smpl, cn, count]  <- tmp_vol * (det(ebs_mean) ) 
                                        ^ (1 / 2)
          cover_ebs[smpl, cn, count] <- as.numeric(sam_siz[smpl]  
                                        * t(asg - parm) %*% qr.solve(ebs_mean) 
                                        %*% (asg - parm) <= crt_val)

          
          tmp_ebs  <- new.sim.int(ebs_mean/sam_siz[smpl], conf = 0.95, 
                                  center = asg)$ints
          leng_ebs <- tmp_ebs[, 2] - tmp_ebs[, 1]
          
          ratio_ibs_ebs[smpl, cn, count, ]  <- leng_ibs/leng_ebs
          marg_sim_cov_ebs[smpl, cn, count] <- as.numeric(sum(tmp_ebs[, 1] 
                                               <= parm &  tmp_ebs[, 2] >= parm)
                                               == nparm)
          marg_volm_ebs[smpl, cn, count]    <- (prod(leng_ebs) /
                                               volm_ebs[smpl, cn, count])
                                               ^(1 / nparm)
         
          # EBS Lugsail calculation 
          
          ebs_mean <- ebs_batch_mean(sg_ct, alp, cns[mk], bt_typ, 2)
          
          forb_ebs_ls[smpl, cn, count]  <- norm(ebs_mean - sigm, 
                                                "F")/ norm(sigm, "F")
          volm_ebs_ls[smpl, cn, count]  <- tmp_vol * (det(ebs_mean) )
                                           ^ (1 / 2)
          cover_ebs_ls[smpl, cn, count] <- as.numeric(sam_siz[smpl]  
                                           * t(asg - parm) %*% qr.solve(ebs_mean) 
                                           %*% (asg - parm) <= crt_val)

          tmp_ebs_ls  <- new.sim.int(ebs_mean/sam_siz[smpl],
                         conf = 0.95, center = asg)$ints
          leng_ebs_ls <- tmp_ebs_ls[, 2] - tmp_ebs_ls[, 1]
          
          marg_sim_cov_ebs_ls[smpl, cn, count] <- as.numeric(sum(tmp_ebs_ls[, 1]
                                                  <= parm &  tmp_ebs_ls[, 2] 
                                                  >= parm) == nparm)
          marg_volm_ebs_ls[smpl, cn, count]    <- (prod(leng_ebs_ls) / 
                                                   volm_ebs_ls[smpl, cn,
                                                   count])^(1 / nparm)          
          
          ratio_ibs_ebs_ls[smpl, cn, count, ] <- leng_ibs/leng_ebs_ls
          ratio_ebs_ls_ebs[smpl, cn, count, ] <- leng_ebs_ls/leng_ebs
          
          count = count + 1           
        }
      }  
    }
    list(forb_ibs[cn, ], 
      volm_ibs[cn, ], 
      cover_ibs[cn, ], 
      cover_orc[cn, ], 
      forb_ebs[, cn, ], 
      volm_ebs[, cn, ],  
      cover_ebs[, cn, ], 
      forb_ebs_ls[, cn, ], 
      volm_ebs_ls[, cn, ],  
      cover_ebs_ls[, cn, ],  
      ratio_ibs_ebs[, cn, , ], 
      ratio_ibs_ebs_ls[, cn, , ], 
      ratio_ebs_ls_ebs[, cn, , ], 
      marg_sim_cov_orc[cn, ], 
      marg_sim_cov_ibs[cn, ], 
      marg_sim_cov_ebs[, cn, ], 
      marg_sim_cov_ebs_ls[, cn, ], 
      volm_orc[cn, ], 
      marg_volm_orc[cn, ], 
      marg_volm_ibs[cn, ], 
      marg_volm_ebs[, cn, ], 
      marg_volm_ebs_ls[, cn, ])
    
  }
  
  for( k in 1 : Rep)
  {
    forb_ibs[k, ]              <- final_values[[1]][[k]]
    volm_ibs[k, ]              <- final_values[[2]][[k]]
    cover_ibs[k, ]             <- final_values[[3]][[k]]
    cover_orc[k, ]             <- final_values[[4]][[k]]
    forb_ebs[, k, ]            <- final_values[[5]][[k]]
    volm_ebs[, k, ]            <- final_values[[6]][[k]]
    cover_ebs[, k, ]           <- final_values[[7]][[k]]
    forb_ebs_ls[, k, ]         <- final_values[[8]][[k]]
    volm_ebs_ls[, k, ]         <- final_values[[9]][[k]]
    cover_ebs_ls[, k, ]        <- final_values[[10]][[k]]
    ratio_ibs_ebs[, k, , ]     <- final_values[[11]][[k]]
    ratio_ibs_ebs_ls[, k, , ]  <- final_values[[12]][[k]]
    ratio_ebs_ls_ebs[, k, , ]  <- final_values[[13]][[k]]
    marg_sim_cov_orc[k, ]      <- final_values[[14]][[k]]
    marg_sim_cov_ibs[k, ]      <- final_values[[15]][[k]]
    marg_sim_cov_ebs[, k, ]    <- final_values[[16]][[k]]
    marg_sim_cov_ebs_ls[, k, ] <- final_values[[17]][[k]]
    volm_orc[k, ]              <- final_values[[18]][[k]]
    marg_volm_orc[k, ]         <- final_values[[19]][[k]]
    marg_volm_ibs[k, ]         <- final_values[[20]][[k]]
    marg_volm_ebs[, k, ]       <- final_values[[21]][[k]]
    marg_volm_ebs_ls[, k, ]    <- final_values[[22]][[k]]
  }
  
  
  
  fil_nam <- paste("out/lad_", nam_matrix, "_n_", n, "_dim_", nparm, 
                   ".RData", sep = "")
  save(marg_sim_cov_orc, marg_sim_cov_ibs, marg_sim_cov_ebs, marg_sim_cov_ebs_ls,
       volm_orc, marg_volm_orc, marg_volm_ibs, marg_volm_ebs, marg_volm_ebs_ls, 
       ratio_ibs_ebs, ratio_ibs_ebs_ls, ratio_ebs_ls_ebs, cover_orc, cover_ibs,
       cover_ebs, cover_ebs_ls, volm_ibs, volm_ebs, volm_ebs_ls, forb_ibs,
       forb_ebs, forb_ebs_ls, file = fil_nam)
}





