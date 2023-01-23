# rm(list = ls())
###########################################################################
## Results for Dimension 5
###########################################################################
nparm <- 5
sq_n <- c("5e4","1e5","2e5","5e5","8e5","1e6","5e6")
names_est <- c("Oracle", "IBS", "c.1b1", "c.1b2", "c.1b3", "Lc.1b1", "Lc.1b2", "Lc.1b3")
dnames <-  list(sq_n, names_est)

cover_all         <- matrix(nrow = length(sq_n), ncol = 8, dimnames = dnames)
marg_covg         <- matrix(nrow = length(sq_n), ncol = 8, dimnames = dnames)
marg_covg_volm    <- matrix(nrow = length(sq_n), ncol = 8, dimnames = dnames)
forb_norm_all     <- matrix(nrow = length(sq_n), ncol = 8, dimnames = dnames)
volm_all          <- matrix(nrow = length(sq_n), ncol = 8, dimnames = dnames)
sd_volm           <- matrix(nrow = length(sq_n), ncol = 8, dimnames = dnames)
frob_all          <- matrix(nrow = length(sq_n), ncol = 8, dimnames = dnames)
sd_forb           <- matrix(nrow = length(sq_n), ncol = 8, dimnames = dnames)
sd_cover          <- matrix(nrow = length(sq_n), ncol = 8, dimnames = dnames)
sd_marg_covg      <- matrix(nrow = length(sq_n), ncol = 8, dimnames = dnames)
sd_marg_covg_volm <- matrix(nrow = length(sq_n), ncol = 8, dimnames = dnames)
sd_forb_norm      <- matrix(nrow = length(sq_n), ncol = 8, dimnames = dnames)

interv_ratio      <- array(dim = c(length(sq_n), 6, nparm), dimnames = list(1 : length(sq_n), 1 : 6, 1 : nparm))
sd_interv_ratio   <-  array(dim = c(length(sq_n), 6, nparm), dimnames = list(1 : length(sq_n), 1 : 6, 1 : nparm))

load("out/lad_indep_n_5e+06_dim_5.RData")



for(n in 1 : length(sq_n))
{
  ##### Calculate COVERAGE RATES  #####
  cover_all[n, 1] <- mean(cover_orc[,n])
  cover_all[n, 2] <- mean(cover_ibs[,n])
  cover_all[n, 3 : 8] <- c(colMeans(cover_ebs[n, , ]), colMeans(cover_ebs_ls[n, , ]))
  Reps <- length(cover_orc)
  sd_cover[n, 1] <- sd( cover_orc[, n]) / sqrt(Reps)
  sd_cover[n, 2] <- sd( cover_ibs[, n]) / sqrt(Reps)
  sd_cover[n, 3 : 8] <- c(apply(cover_ebs[n, , ], 2, sd), apply(cover_ebs_ls[n, , ], 2, sd)) / sqrt(Reps)

  ##### Calculate Frobenius Norm  #####    
  frob_all[n, 1] <- 0
  frob_all[n, 2] <- mean( forb_ibs[, n])
  frob_all[n, 3 : 8] <- c(colMeans(forb_ebs[n, , ]), colMeans(forb_ebs_ls[n, , ]))
  
  sd_forb[n, 1] <- 0
  sd_forb[n, 2] <- sd( forb_ibs[, n]) / sqrt(Reps)
  sd_forb[n, 3 : 8] <- c(apply(forb_ebs[n, , ], 2, sd), apply(forb_ebs_ls[n, , ], 2, sd)) / sqrt(Reps)
 
 ##### MARGINAL  SIMULTANEOUS COVERAGE #####

  marg_covg[n, 1] <- mean( marg_sim_cov_orc[,n])
  marg_covg[n, 2] <- mean( marg_sim_cov_ibs[,n])
  marg_covg[n, 3 : 8] <- c(colMeans(marg_sim_cov_ebs[n, , ]), colMeans(marg_sim_cov_ebs_ls[n, , ]))
  Reps <- length(cover_orc)
  sd_marg_covg[n, 1] <- sd( marg_sim_cov_orc[, n]) / sqrt(Reps)
  sd_marg_covg[n, 2] <- sd( marg_sim_cov_ibs[, n]) / sqrt(Reps)
  sd_marg_covg[n, 3 : 8] <- c(apply(marg_sim_cov_ebs[n, , ], 2, sd), apply(marg_sim_cov_ebs_ls[n, , ], 2, sd)) / sqrt(Reps)
  
  
  ##### MEASURE OF VOLUME COMPROMISED DUE TO MARGINALISING SIMULTANEOUS ASSESSMENT #####
  marg_covg_volm[n, 1] <- mean( marg_volm_orc[,n])
  marg_covg_volm[n, 2] <- mean( marg_volm_ibs[,n])
  marg_covg_volm[n, 3 : 8] <- c(colMeans(marg_volm_ebs[n, , ]), colMeans(marg_volm_ebs_ls[n, , ]))
  Reps <- length(cover_orc)
  sd_marg_covg_volm[n, 1] <- sd( marg_volm_orc[, n]) / sqrt(Reps)
  sd_marg_covg_volm[n, 2] <- sd( marg_volm_ibs[, n]) / sqrt(Reps)
  sd_marg_covg_volm[n, 3 : 8] <- c(apply(marg_volm_ebs[n, , ], 2, sd), apply(marg_volm_ebs_ls[n, , ], 2, sd)) / sqrt(Reps)
  
  
  for( l in 1 : 3)
  {
    
    for(m in 1 : 5)
    {
      interv_ratio[n, l, m] <- mean( ratio_ibs_ebs[n, ,l ,m])
      interv_ratio[n, 3 + l, m] <- mean( ratio_ibs_ebs_ls[n, ,l ,m])
      #interv_ratio[n, 3 + l, m] <- mean( ratio_ebs_ls_ebs[n, ,l ,m])
      Reps <- length(cover_orc)
      sd_interv_ratio[n, l, m] <- sd( ratio_ibs_ebs[n, ,l ,m]) / sqrt(Reps)
      sd_interv_ratio[n, 3 + l, m] <- sd( ratio_ibs_ebs_ls[n, ,l ,m]) / sqrt(Reps)
     # sd_interv_ratio[n, 3 + l, m] <- sd( ratio_ebs_ls_ebs[n, ,l ,m]) / sqrt(Reps)
      
    }
    
  }

}


##############################
## Making Plots
sq_n <- log10(c(5e4, 1e5, 2e5, 5e5, 8e5,  1e6, 5e6))


col_choic <- c("red", "steelblue", "brown", "tomato", "green", "brown", "tomato", "green")
lin_typ <- c(6, 3, rep(1,3), rep(2,3))
names_var <- c("Oracle", "IBS", expression(paste("c = 0.1, ", beta[1])), expression(paste("c = 0.1, ", beta[2])), expression(paste("c = 0.1, ", beta[3])), expression(paste("LS c = 0.1, ", beta[1])), expression(paste("LS c = 0.1, ", beta[2])), expression(paste("LS c = 0.1, ", beta[3])))

# constant c = 1 is not needed as batch size is too large, beta2 and beta3 are not needed
index <- c(1, 2, 3, 6)

names_var2 <- c("Oracle", "IBS", "EBS", "Lugsail-EBS")

pdf("out/Ladmcover_dim5.pdf", height = 6, width = 7)
plot(sq_n, cover_all[, 1], type = "n", lwd = 2, lty = 6, xlab = "Log Sample Size", ylab = "Multivariate Coverage Rate", col = "red", ylim = range(c(cover_all[, index], 1)))
for(k in index)
{
  lines(sq_n, cover_all[, k],  lwd = 2, col = col_choic[k], lty = lin_typ[k], type = "b")
  arrows(sq_n, cover_all[, k], sq_n, cover_all[, k] + sd_cover[ ,k], length = 0.05, angle=90, col = col_choic[k], lwd = 2, lty = 1 )
  arrows(sq_n, cover_all[, k], sq_n, cover_all[, k] - sd_cover[ ,k], length = 0.05, angle=90, col = col_choic[k], lwd = 2, lty = 1 )   
}
legend("bottom", legend = names_var2,
       col = col_choic[index], lwd = 2, lty = lin_typ[index], cex=1,
       box.lty=0, box.lwd=1, ncol = 2)
dev.off()


pdf("out/Ladfrob_dim5.pdf", height = 6, width = 7)
plot(sq_n, frob_all[, 1], type = "n", lwd = 2, lty = 6, xlab = "Log Sample Size", ylab = "Relative Frobenius Norm", col = "red", ylim = range(c(frob_all[, index], -.2)))
for(k in index)
{
  lines(sq_n, frob_all[, k],  lwd = 2, col = col_choic[k], lty = lin_typ[k], type = "b")
  arrows(sq_n, frob_all[, k], sq_n, frob_all[, k] + sd_forb[ ,k], length = 0.05, angle=90, col = col_choic[k], lwd = 2, lty = 1 )
  arrows(sq_n, frob_all[, k], sq_n, frob_all[, k] - sd_forb[ ,k], length = 0.05, angle=90, col = col_choic[k], lwd = 2, lty = 1 )   
}
legend("bottom", legend = names_var2,
       col = col_choic[index], lwd = 2, lty = lin_typ[index], cex=1,
       box.lty=0, box.lwd=1, ncol = 2)
dev.off()


pdf("out/Ladsimul_dim5.pdf", height = 6, width = 7)
plot(sq_n, marg_covg[, 1], type = "n", lwd = 2, lty = 6, xlab = "Log Sample Size", ylab = "Simultaneous Marginal Coverage Rate", col = "red", ylim = range(c(cover_all[, index], 1)))
for(k in index)
{
  lines(sq_n, marg_covg[, k],  lwd = 2, col = col_choic[k], lty = lin_typ[k], type = "b")
  arrows(sq_n, marg_covg[, k], sq_n, marg_covg[, k] + sd_marg_covg[ ,k], length = 0.05, angle=90, col = col_choic[k], lwd = 2, lty = 1 )
  arrows(sq_n, marg_covg[, k], sq_n, marg_covg[, k] - sd_marg_covg[ ,k], length = 0.05, angle=90, col = col_choic[k], lwd = 2, lty = 1 )   
}

legend("bottom", legend = names_var2,
       col = col_choic[index], lwd = 2, lty = lin_typ[index], cex=1,
       box.lty=0, box.lwd=1, ncol = 2)
dev.off()



pdf("out/Ladrelvol_dim5.pdf", height = 6, width = 7)
plot(sq_n, marg_covg_volm[, 1], type = "n", lwd = 2, lty = 6, xlab = "Log Sample Size", ylab = "Volume Ratios to the pth root", col = "red", ylim = range(c(marg_covg_volm[, index], 1.02)))
for(k in index)
{
  lines(sq_n, marg_covg_volm[, k],  lwd = 2, col = col_choic[k], lty = lin_typ[k], type = "b")
  arrows(sq_n, marg_covg_volm[, k], sq_n, marg_covg_volm[, k] + sd_marg_covg_volm[ ,k], length = 0.05, angle=90, col = col_choic[k], lwd = 2, lty = 1 )
  arrows(sq_n, marg_covg_volm[, k], sq_n, marg_covg_volm[, k] - sd_marg_covg_volm[ ,k], length = 0.05, angle=90, col = col_choic[k], lwd = 2, lty = 1 )   
}

legend("bottom", legend = names_var2,
       col = col_choic[index], lwd = 2, lty = lin_typ[index], cex=1,
       box.lty=0, box.lwd=1, ncol = 2)
dev.off()





#######################################################################
## Results for Dimension 20 
#######################################################################

load("out/lad_indep_n_5e+06_dim_20.RData")
nparm <- 20
sq_n <- c("5e4","1e5","2e5","5e5","8e5","1e6","5e6")

names_est <- c("Oracle", "IBS", "c.1b1", "c.1b2", "c.1b3", "Lc.1b1", "Lc.1b2", "Lc.1b3")
dnames <-  list(sq_n, names_est)

cover_all         <- matrix(nrow = length(sq_n), ncol = 8, dimnames = dnames)
marg_covg         <- matrix(nrow = length(sq_n), ncol = 8, dimnames = dnames)
marg_covg_volm    <- matrix(nrow = length(sq_n), ncol = 8, dimnames = dnames)
forb_norm_all     <- matrix(nrow = length(sq_n), ncol = 8, dimnames = dnames)
volm_all          <- matrix(nrow = length(sq_n), ncol = 8, dimnames = dnames)
sd_volm           <- matrix(nrow = length(sq_n), ncol = 8, dimnames = dnames)
frob_all          <- matrix(nrow = length(sq_n), ncol = 8, dimnames = dnames)
sd_forb           <- matrix(nrow = length(sq_n), ncol = 8, dimnames = dnames)
sd_cover          <- matrix(nrow = length(sq_n), ncol = 8, dimnames = dnames)
sd_marg_covg      <- matrix(nrow = length(sq_n), ncol = 8, dimnames = dnames)
sd_marg_covg_volm <- matrix(nrow = length(sq_n), ncol = 8, dimnames = dnames)
sd_forb_norm      <- matrix(nrow = length(sq_n), ncol = 8, dimnames = dnames)

interv_ratio      <- array(dim = c(length(sq_n), 6, nparm), dimnames = list(1 : length(sq_n), 1 : 6, 1 : nparm))
sd_interv_ratio   <- array(dim = c(length(sq_n), 6, nparm), dimnames = list(1 : length(sq_n), 1 : 6, 1 : nparm))


for(n in 1 : length(sq_n))
{
  ##### Calculate COVERAGE RATES  #####
  cover_all[n, 1] <- mean(cover_orc[,n])
  cover_all[n, 2] <- mean(cover_ibs[,n])
  cover_all[n, 3 : 8] <- c(colMeans(cover_ebs[n, , ]), colMeans(cover_ebs_ls[n, , ]))
  Reps <- length(cover_orc)
  sd_cover[n, 1] <- sd( cover_orc[, n]) / sqrt(Reps)
  sd_cover[n, 2] <- sd( cover_ibs[, n]) / sqrt(Reps)
  sd_cover[n, 3 : 8] <- c(apply(cover_ebs[n, , ], 2, sd), apply(cover_ebs_ls[n, , ], 2, sd)) / sqrt(Reps)

  ##### Calculate Frobenius Norm  #####    
  frob_all[n, 1] <- 0
  frob_all[n, 2] <- mean( forb_ibs[, n])
  frob_all[n, 3 : 8] <- c(colMeans(forb_ebs[n, , ]), colMeans(forb_ebs_ls[n, , ]))
  
  sd_forb[n, 1] <- 0
  sd_forb[n, 2] <- sd( forb_ibs[, n]) / sqrt(Reps)
  sd_forb[n, 3 : 8] <- c(apply(forb_ebs[n, , ], 2, sd), apply(forb_ebs_ls[n, , ], 2, sd)) / sqrt(Reps)
 
 ##### MARGINAL  SIMULTANEOUS COVERAGE #####

  marg_covg[n, 1] <- mean( marg_sim_cov_orc[,n])
  marg_covg[n, 2] <- mean( marg_sim_cov_ibs[,n])
  marg_covg[n, 3 : 8] <- c(colMeans(marg_sim_cov_ebs[n, , ]), colMeans(marg_sim_cov_ebs_ls[n, , ]))
  Reps <- length(cover_orc)
  sd_marg_covg[n, 1] <- sd( marg_sim_cov_orc[, n]) / sqrt(Reps)
  sd_marg_covg[n, 2] <- sd( marg_sim_cov_ibs[, n]) / sqrt(Reps)
  sd_marg_covg[n, 3 : 8] <- c(apply(marg_sim_cov_ebs[n, , ], 2, sd), apply(marg_sim_cov_ebs_ls[n, , ], 2, sd)) / sqrt(Reps)
  
  
  ##### MEASURE OF VOLUME COMPROMISED DUE TO MARGINALISING SIMULTANEOUS ASSESSMENT #####
  marg_covg_volm[n, 1] <- mean( marg_volm_orc[,n])
  marg_covg_volm[n, 2] <- mean( marg_volm_ibs[,n])
  marg_covg_volm[n, 3 : 8] <- c(colMeans(marg_volm_ebs[n, , ]), colMeans(marg_volm_ebs_ls[n, , ]))
  Reps <- length(cover_orc)
  sd_marg_covg_volm[n, 1] <- sd( marg_volm_orc[, n]) / sqrt(Reps)
  sd_marg_covg_volm[n, 2] <- sd( marg_volm_ibs[, n]) / sqrt(Reps)
  sd_marg_covg_volm[n, 3 : 8] <- c(apply(marg_volm_ebs[n, , ], 2, sd), apply(marg_volm_ebs_ls[n, , ], 2, sd)) / sqrt(Reps)
  
  
  for( l in 1 : 3)
  {
    
    for(m in 1 : 5)
    {
      interv_ratio[n, l, m] <- mean( ratio_ibs_ebs[n, ,l ,m])
      interv_ratio[n, 3 + l, m] <- mean( ratio_ibs_ebs_ls[n, ,l ,m])
      #interv_ratio[n, 12 + l, m] <- mean( ratio_ebs_ls_ebs[n, ,l ,m])
      Reps <- length(cover_orc)
      sd_interv_ratio[n, l, m] <- sd( ratio_ibs_ebs[n, ,l ,m]) / sqrt(Reps)
      sd_interv_ratio[n, 3 + l, m] <- sd( ratio_ibs_ebs_ls[n, ,l ,m]) / sqrt(Reps)
     # sd_interv_ratio[n, 12 + l, m] <- sd( ratio_ebs_ls_ebs[n, ,l ,m]) / sqrt(Reps)
      
    }
    
  }

}




sq_n <- log10(c(5e4, 1e5, 2e5, 5e5, 8e5,  1e6, 5e6))


col_choic <- c("red", "steelblue", "brown", "tomato", "green", "brown", "tomato", "green")
lin_typ <- c(6, 3, rep(1,3), rep(2,3))
names_var <- c("Oracle", "IBS", expression(paste("c = 0.1, ", beta[1])), expression(paste("c = 0.1, ", beta[2])), expression(paste("c = 0.1, ", beta[3])), expression(paste("LS c = 0.1, ", beta[1])), expression(paste("LS c = 0.1, ", beta[2])), expression(paste("LS c = 0.1, ", beta[3])))

index <- c(1, 2, 3, 6)

pdf("out/Ladmcover_dim20.pdf", height = 6, width = 7)
plot(sq_n, cover_all[, 1], type = "n", lwd = 2, lty = 6, xlab = "Log Sample Size", ylab = "Multivariate Coverage Rate", col = "red", ylim = range(c(cover_all[, index], 1, -.2)))
for(k in index)
{
  lines(sq_n, cover_all[, k],  lwd = 2, col = col_choic[k], lty = lin_typ[k], type = "b")
  arrows(sq_n, cover_all[, k], sq_n, cover_all[, k] + sd_cover[ ,k], length = 0.05, angle=90, col = col_choic[k], lwd = 2, lty = 1 )
  arrows(sq_n, cover_all[, k], sq_n, cover_all[, k] - sd_cover[ ,k], length = 0.05, angle=90, col = col_choic[k], lwd = 2, lty = 1 )   
}

legend("bottom", legend = names_var2,
       col = col_choic[index], lwd = 2, lty = lin_typ[index], cex=1,
       box.lty=0, box.lwd=1, ncol = 2)
dev.off()



pdf("out/Ladsimul_dim20.pdf", height = 6, width = 7)
plot(sq_n, marg_covg[, 1], type = "n", lwd = 2, lty = 6, xlab = "Log Sample Size", ylab = "Simultaneous Marginal Coverage Rate", col = "red", ylim = range(c(cover_all[, index], 1)))
for(k in index)
{
  lines(sq_n, marg_covg[, k],  lwd = 2, col = col_choic[k], lty = lin_typ[k], type = "b")
  arrows(sq_n, marg_covg[, k], sq_n, marg_covg[, k] + sd_marg_covg[ ,k], length = 0.05, angle=90, col = col_choic[k], lwd = 2, lty = 1 )
  arrows(sq_n, marg_covg[, k], sq_n, marg_covg[, k] - sd_marg_covg[ ,k], length = 0.05, angle=90, col = col_choic[k], lwd = 2, lty = 1 )   
}

legend("bottom", legend = names_var2,
       col = col_choic[index], lwd = 2, lty = lin_typ[index], cex=1,
       box.lty=0, box.lwd=1, ncol = 2)
dev.off()


pdf("out/Ladfrob_dim20.pdf", height = 6, width = 7)
plot(sq_n, frob_all[, 1], type = "n", lwd = 2, lty = 6, xlab = "Log Sample Size", ylab = "Relative Frobenius Norm", col = "red", ylim = range(c(frob_all[, index], -.2)))
for(k in index)
{
  lines(sq_n, frob_all[, k],  lwd = 2, col = col_choic[k], lty = lin_typ[k], type = "b")
  arrows(sq_n, frob_all[, k], sq_n, frob_all[, k] + sd_forb[ ,k], length = 0.05, angle=90, col = col_choic[k], lwd = 2, lty = 1 )
  arrows(sq_n, frob_all[, k], sq_n, frob_all[, k] - sd_forb[ ,k], length = 0.05, angle=90, col = col_choic[k], lwd = 2, lty = 1 )   
}
legend("bottom", legend = names_var2,
       col = col_choic[index], lwd = 2, lty = lin_typ[index], cex=1,
       box.lty=0, box.lwd=1, ncol = 2)
dev.off()

pdf("out/Ladrelvol_dim20.pdf", height = 6, width = 7)
plot(sq_n, marg_covg_volm[, 1], type = "n", lwd = 2, lty = 6, xlab = "Log Sample Size", ylab = "Volume Ratios to the pth root", col = "red", ylim = range(c(marg_covg_volm[, index], 1)))
for(k in index)
{
  lines(sq_n, marg_covg_volm[, k],  lwd = 2, col = col_choic[k], lty = lin_typ[k], type = "b")
  arrows(sq_n, marg_covg_volm[, k], sq_n, marg_covg_volm[, k] + sd_marg_covg_volm[ ,k], length = 0.05, angle=90, col = col_choic[k], lwd = 2, lty = 1 )
  arrows(sq_n, marg_covg_volm[, k], sq_n, marg_covg_volm[, k] - sd_marg_covg_volm[ ,k], length = 0.05, angle=90, col = col_choic[k], lwd = 2, lty = 1 )   
}

legend("bottom", legend = names_var2,
       col = col_choic[index], lwd = 2, lty = lin_typ[index], cex=1,
       box.lty=0, box.lwd=1, ncol = 2)
dev.off()





