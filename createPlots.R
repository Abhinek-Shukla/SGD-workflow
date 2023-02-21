###########################################################
## This file contains all the code that produces
## all the plots and tables in the paper.
## The plots are made in plots/ and the tables are in the
## the .Rout file
###########################################################




###########################################################
## Figure 1: Pictorial demonstration of the batching strategy
###########################################################

n <- seq(1e2, 1e6, by = 200)
alpha <- .51
beta <- (1 + alpha)/2
c <- .1

bnfun <- function(this.n, beta)
{
  upp <- c*this.n^beta
  bn.vector <- 2^seq(1, 100, by = 1)
  gamma <- which(bn.vector >= upp)[1]
  return(c(2^gamma, this.n/(2^gamma)))
}

out <- sapply(n, bnfun, beta = beta)
bn <- out[1, ]
an <- out[2, ]

pdf("plots/batchsize.pdf", height = 5, width = 5)
plot(n, 2*c*n^(beta), type = 'l', 
     ylab = "Batch Size", xlab = "Iteration Length",
     lty = 3, col = "purple", lwd = 1.1)
lines(n, c*n^(beta), lty = 3, col = "purple", lwd = 1.1)
lines(n, bn, col = "black")
dev.off()

pdf("plots/nbatch.pdf", height = 5, width = 5)
plot(n, n^(1-beta)/c, type = 'l',
     ylab = "Number of Batches", xlab = "Iteration Length",
     lty = 3, col = "purple", lwd = 1.1)
lines(n, n^(1 - beta)/c/2, lty = 3, col = "purple", lwd = 1.1)
lines(n, an, col = "black")
dev.off()




###########################################################
## Figure 2: Bias of EBS and Lugsail-EBS for mean model
###########################################################

load("bias/out/ebs_lugsail_bias.RData")

########## plot for alpha = 0.51
pdf(file = "plots/ebs_bias_51.pdf")

# EBS bias plot against the sample size
plot(sam_siz, out_ebs[,1], type = 'b', pch = 19, col = 'black', 
     ylim = c(min(out_ebs), 0), 
     ylab="Bias upto a constant multiple", xlab="Sample size", 
     cex.main=1.25, cex.lab=1.5, cex.axis=1.1)

# Add a line for Lugsail-EBS bias
lines(sam_siz, out_ebs[,2], type = 'b', pch = 18, col = "blue", 
      lty=2)

# Add legends
legend("bottomright", legend=c("EBS", "Lugsail-EBS"),
       col=c("black", "blue"), lty=1:2, cex= 1)

dev.off()


########## plot for alpha = 0.75
pdf(file = "plots/ebs_bias_75.pdf")
# EBS bias plot against the sample size
plot(sam_siz, out_ebs[,3], type = 'b', pch = 19, col = 'black',
     ylim = c(min(out_ebs), 0), ylab="Bias upto a constant multiple", 
     xlab="Sample size", cex.main=1.25, cex.lab=1.5, cex.axis=1.1)

# Add a line for Lugsail-EBS bias
lines(sam_siz, out_ebs[,4], type = 'b', pch = 18, col = "blue", 
      lty=2)

# Add legends
legend("topright", legend=c("EBS", "Lugsail-EBS"),
       col=c("black", "blue"), lty=1:2, cex= 1)

dev.off()


###########################################################
## Figure 3: Simulatenous Confidence Regions
##            visualizations
###########################################################



###########################################################
## Simulation and Examples
###########################################################

############### Section 6.1  ##############################


#### Dimension 5 - independent ###########
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


load("Examples/linear/out/linear_indep_n_5e+06_dim_5.RData")

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


###########################################################
## Figure 4: Coverage and volume for 5 dimensions (independent)
###########################################################
sq_n <- log10(c(5e4, 1e5, 2e5, 5e5, 8e5,  1e6, 5e6))

col_choic <- c("red", "steelblue", "brown", "tomato", "green", "brown", "tomato", "green")
lin_typ <- c(6, 3, rep(1,3), rep(2,3))
names_var <- c("Oracle", "IBS", expression(paste("c = 0.1, ", beta[1])), expression(paste("c = 0.1, ", beta[2])), expression(paste("c = 0.1, ", beta[3])), expression(paste("LS c = 0.1, ", beta[1])), expression(paste("LS c = 0.1, ", beta[2])), expression(paste("LS c = 0.1, ", beta[3])))

# constant c = 1 is not needed and, beta2 and beta3 are  not needed
index <- c(1, 2, 3, 6)

names_var2 <- c("Oracle", "IBS", "EBS", "Lugsail-EBS")

pdf("plots/Lmcover_dim5.pdf", height = 6, width = 7)
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


pdf("plots/Lfrob_dim5.pdf", height = 6, width = 7)
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


pdf("plots/Lsimul_dim5.pdf", height = 6, width = 7)
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



pdf("plots/Lrelvol_dim5.pdf", height = 6, width = 7)
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





#### Dimension 5 - Toeplitz ###########
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


load("Examples/linear/out/linear_toep_n_5e+06_dim_5.RData")

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


###########################################################
## Figure 5: Coverage and volume for 5 dimensions (Toeplitz)
###########################################################
sq_n <- log10(c(5e4, 1e5, 2e5, 5e5, 8e5,  1e6, 5e6))

col_choic <- c("red", "steelblue", "brown", "tomato", "green", "brown", "tomato", "green")
lin_typ <- c(6, 3, rep(1,3), rep(2,3))
names_var <- c("Oracle", "IBS", expression(paste("c = 0.1, ", beta[1])), expression(paste("c = 0.1, ", beta[2])), expression(paste("c = 0.1, ", beta[3])), expression(paste("LS c = 0.1, ", beta[1])), expression(paste("LS c = 0.1, ", beta[2])), expression(paste("LS c = 0.1, ", beta[3])))

# constant c = 1 is not needed and, beta2 and beta3 are  not needed
index <- c(1, 2, 3, 6)

names_var2 <- c("Oracle", "IBS", "EBS", "Lugsail-EBS")

pdf("plots/Lmcover_dim5_toep.pdf", height = 6, width = 7)
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


pdf("plots/Lfrob_dim5_toep.pdf", height = 6, width = 7)
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


pdf("plots/Lsimul_dim5_toep.pdf", height = 6, width = 7)
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



pdf("plots/Lrelvol_dim5_toep.pdf", height = 6, width = 7)
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







#### Dimension 5 - Equivariance ###########
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


load("Examples/linear/out/linear_equiv_n_5e+06_dim_5.RData")

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


###########################################################
## Figure 6: Coverage and volume for 5 dimensions (Equivariance)
###########################################################
sq_n <- log10(c(5e4, 1e5, 2e5, 5e5, 8e5,  1e6, 5e6))

col_choic <- c("red", "steelblue", "brown", "tomato", "green", "brown", "tomato", "green")
lin_typ <- c(6, 3, rep(1,3), rep(2,3))
names_var <- c("Oracle", "IBS", expression(paste("c = 0.1, ", beta[1])), expression(paste("c = 0.1, ", beta[2])), expression(paste("c = 0.1, ", beta[3])), expression(paste("LS c = 0.1, ", beta[1])), expression(paste("LS c = 0.1, ", beta[2])), expression(paste("LS c = 0.1, ", beta[3])))

# constant c = 1 is not needed and, beta2 and beta3 are  not needed
index <- c(1, 2, 3, 6)

names_var2 <- c("Oracle", "IBS", "EBS", "Lugsail-EBS")

pdf("plots/Lmcover_dim5_equiv.pdf", height = 6, width = 7)
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


pdf("plots/Lfrob_dim5_equiv.pdf", height = 6, width = 7)
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


pdf("plots/Lsimul_dim5_equiv.pdf", height = 6, width = 7)
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



pdf("plots/Lrelvol_dim5_equiv.pdf", height = 6, width = 7)
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








###### Dimension 20 - independent ###########
load("Examples/linear/out/linear_indep_n_5e+06_dim_20.RData")
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

pdf("plots/Lmcover_dim20.pdf", height = 6, width = 7)
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



pdf("plots/Lsimul_dim20.pdf", height = 6, width = 7)
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


pdf("plots/Lfrob_dim20.pdf", height = 6, width = 7)
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

pdf("plots/Lrelvol_dim20.pdf", height = 6, width = 7)
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





