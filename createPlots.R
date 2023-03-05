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
     lty = 3, col = "purple", lwd = 1.8)
lines(n, c*n^(beta), lty = 3, col = "purple", lwd = 1.1)
lines(n, bn, col = "black")
dev.off()

pdf("plots/nbatch.pdf", height = 5, width = 5)
plot(n, n^(1-beta)/c, type = 'l',
     ylab = "Number of Batches", xlab = "Iteration Length",
     lty = 3, col = "purple", lwd = 1.1)
lines(n, n^(1 - beta)/c/2, lty = 3, col = "purple", lwd = 1.8)
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

###############################################################################
##generating plots to showcase upper and lower bounds
set.seed(100)
library(mvtnorm)
library(ellipse)


alpha = .1
conf = 1 - alpha
p = 2

##p: perfect correlation
##h: rho = .75
##i: independent
Sigma.pcor = matrix(c(9,5.5,5.5,4),nrow=2)
Sigma.half = matrix(c(9,4,4,4),nrow=2)
Sigma.ind = matrix(c(9,0,0,4),nrow=2)
s.x = 3
s.y = 2


n = 100000
X.p = rmvnorm(n,c(0,0),Sigma.pcor)
X.h = rmvnorm(n,c(0,0),Sigma.half)
X.i = rmvnorm(n,c(0,0),Sigma.ind)

##m: marginal
##b: bonferroni
##f: intermediate
z.m = qnorm(1-alpha)
z.b = qnorm(1 - (alpha/2)/2)

z.f = qnorm(1-alpha/2)

hwu.1 = z.m*s.x
hwo.1 = z.b*s.x
hwu.2 = z.m*s.y
hwo.2 = z.b*s.y

hwf.1 = z.f*s.x 
hwf.2 = z.f*s.y 

checks = 1000

ratio.12 = hwu.1/hwu.2
#xs = seq(hwu.1,hwo.1,length=checks)
#ys = xs/ratio.12

bl.u = c(-hwu.1,-hwu.2)
tl.u = c(-hwu.1,hwu.2)
br.u = c(hwu.1,-hwu.2)
tr.u = c(hwu.1,hwu.2)

bl.o = c(-hwo.1,-hwo.2)
tl.o = c(-hwo.1,hwo.2)
br.o = c(hwo.1,-hwo.2)
tr.o = c(hwo.1,hwo.2)

bl.f = c(-hwf.1,-hwf.2)
tl.f = c(-hwf.1,hwf.2)
br.f = c(hwf.1,-hwf.2)
tr.f = c(hwf.1,hwf.2)

ends.u = rbind(bl.u,br.u,tl.u,tr.u)
ends.o = rbind(bl.o,br.o,tl.o,tr.o)
ends.f = rbind(bl.f,br.f,tl.f,tr.f)


pdf("plots/marginal_ints.pdf", width = 5, height = 5)

#par(mar=c(5,5,4,6),xpd = TRUE)
par(mar = c(5,5,1,1))
plot(ends.o,col="blue4",pch=12, xlab="First component", ylab="Second component",
     main=" ",font.main=1,
     xlim=c(-7,7), ylim=c(-6,5), frame.plot = TRUE)
points(ends.u, col="black",pch=12)
points(ends.f, col="red2",pch=12)
#lines(xs,ys, col="black")

lines(rbind(bl.u,br.u),col="black", lty=c(5))
lines(rbind(bl.u,tl.u),col="black", lty=c(5))
lines(rbind(br.u,tr.u),col="black", lty=c(5))
lines(rbind(tl.u,tr.u),col="black", lty=c(5))


lines(rbind(bl.o,br.o),col="blue", lty=c(2))
lines(rbind(bl.o,tl.o),col="blue", lty=c(2))
lines(rbind(br.o,tr.o),col="blue", lty=c(2))
lines(rbind(tl.o,tr.o),col="blue", lty=c(2))

lines(rbind(bl.f,br.f),col="red", lty=c(9))
lines(rbind(bl.f,tl.f),col="red", lty=c(9))
lines(rbind(br.f,tr.f),col="red", lty=c(9))
lines(rbind(tl.f,tr.f),col="red", lty=c(9))


legend("bottom", legend=c("Upper Bound", "Simultaneous", "Lower Bound",
                          "Ellipse"),
       col=c("blue", "red", 1, "magenta"),
       lty=c(2,9,5,1), pch = c(12,12,12,NA), box.lwd = 3, cex=0.8, ncol = 2, 
       bty = 'n')


##Sigma.ind, Sigma.pcor, Sigma.half
lines(ellipse(Sigma.half, level = .9), lty=1, col="magenta")
dev.off()





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





###########################################################
## Figure 7: QQ-plot for Dim 5 (independent)
###########################################################

load("Examples/linear/out/qq_dim_5_sam_50000.Rdata")
nparm <- 5
n <- 5e4

p1 <- pnorm(seq(-.5,.5, length = 1e2))
ibs_new <- quantile(scale(ibs_batches), probs = p1)
ebs_new <- quantile(scale(ebs_batches), probs = p1)
theory <- qnorm(p1)

pdf(file = paste("plots/qq_dim_", nparm, "_sam_", n, ".pdf"),
height = 5, width = 5)
plot(theory, ibs_new, xlab = "Theoretical Quantiles", 
  ylab = "Sample Quantiles", ylim = range(theory), type = "n")
lines(theory, ibs_new, pch = 1, col = "steelblue")
lines(theory, ebs_new, pch = 1, col = "brown", lwd = 2)
lines(theory, theory, col = "black", lwd = 1, lty = 2)
legend("bottomright", legend = c("IBS", "EBS"), 
  col = c("steelblue", "brown"), 
  lty = 1, box.lty=0, box.lwd=1)
dev.off()



load("Examples/linear/out/qq_dim_5_sam_1e+06.Rdata")
nparm <- 5
n <- 1e6

p1 <- pnorm(seq(-.5,.5, length = 1e2))
ibs_new <- quantile(scale(ibs_batches), probs = p1)
ebs_new <- quantile(scale(ebs_batches), probs = p1)
theory <- qnorm(p1)

pdf(file = paste("plots/qq_dim_", nparm, "_sam_", n, ".pdf"),
height = 5, width = 5)
plot(theory, ibs_new, xlab = "Theoretical Quantiles", 
  ylab = "Sample Quantiles", ylim = range(theory), type = "n")
lines(theory, ibs_new, pch = 1, col = "steelblue")
lines(theory, ebs_new, pch = 1, col = "brown", lwd = 2)
lines(theory, theory, col = "black", lwd = 1, lty = 2)
legend("bottomright", legend = c("IBS", "EBS"), 
  col = c("steelblue", "brown"), 
  lty = 1, box.lty=0, box.lwd=1)
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
  
  
  for(l in 1:3)
  {
    
    for(m in 1:5)
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


###########################################################
## Figure 7: Coverage and volume for 20 dimensions (indep)
###########################################################

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





##########################################################
## LAD Example
##########################################################



##########################################################
## Figure 11: Coverage and Volume for Dimension 20 (indep)
##########################################################
rm(list = ls())
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

interv_ratio      <- array(dim = c(length(sq_n), 6, nparm), dimnames = list(1:length(sq_n), 1:6, 1:nparm))
sd_interv_ratio   <-  array(dim = c(length(sq_n), 6, nparm), dimnames = list(1:length(sq_n), 1:6, 1:nparm))

load("Examples/linear_lad/out/lad_indep_n_5e+06_dim_20.RData")

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

pdf("plots/Ladmcover_dim20.pdf", height = 6, width = 7)
plot(sq_n, cover_all[, 1], type = "n", lwd = 2, lty = 6, xlab = "Log Sample Size", ylab = "Multivariate Coverage Rate", col = "red", ylim = range(c(-.2, cover_all[, index], 1)))
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


pdf("plots/Ladfrob_dim20.pdf", height = 6, width = 7)
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


pdf("plots/Ladsimul_dim20.pdf", height = 6, width = 7)
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



pdf("plots/Ladrelvol_dim20.pdf", height = 6, width = 7)
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






##########################################################
## Figure 12: Coverage and Volume for Dimension 20 (indep)
##########################################################
rm(list = ls())
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

interv_ratio      <- array(dim = c(length(sq_n), 6, nparm), dimnames = list(1:length(sq_n), 1:6, 1:nparm))
sd_interv_ratio   <-  array(dim = c(length(sq_n), 6, nparm), dimnames = list(1:length(sq_n), 1:6, 1:nparm))

load("Examples/linear_lad/out/lad_toep_n_5e+06_dim_20.RData")

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

pdf("plots/Ladmcover_dim20_toep.pdf", height = 6, width = 7)
plot(sq_n, cover_all[, 1], type = "n", lwd = 2, lty = 6, xlab = "Log Sample Size", ylab = "Multivariate Coverage Rate", col = "red", ylim = range(c(-.2, cover_all[, index], 1)))
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


pdf("plots/Ladfrob_dim20_toep.pdf", height = 6, width = 7)
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


pdf("plots/Ladsimul_dim20_toep.pdf", height = 6, width = 7)
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



pdf("plots/Ladrelvol_dim20_toep.pdf", height = 6, width = 7)
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







##########################################################
## Figure 13: Coverage and Volume for Dimension 20 (equiv)
##########################################################
rm(list = ls())
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

interv_ratio      <- array(dim = c(length(sq_n), 6, nparm), dimnames = list(1:length(sq_n), 1:6, 1:nparm))
sd_interv_ratio   <-  array(dim = c(length(sq_n), 6, nparm), dimnames = list(1:length(sq_n), 1:6, 1:nparm))

load("Examples/linear_lad/out/lad_equiv_n_5e+06_dim_20.RData")

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

pdf("plots/Ladmcover_dim20_equiv.pdf", height = 6, width = 7)
plot(sq_n, cover_all[, 1], type = "n", lwd = 2, lty = 6, xlab = "Log Sample Size", ylab = "Multivariate Coverage Rate", col = "red", ylim = range(c(-.2, cover_all[, index], 1)))
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


pdf("plots/Ladfrob_dim20_equiv.pdf", height = 6, width = 7)
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


pdf("plots/Ladsimul_dim20_equiv.pdf", height = 6, width = 7)
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



pdf("plots/Ladrelvol_dim20_equiv.pdf", height = 6, width = 7)
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



##########################################################
## Figure xx: Logistic Example
##########################################################

load("Examples/logistic/out/logistic_pred.Rdata")
pdf("plots/logistic_misclass.pdf", height = 6, width = 6)
plot(cutoffs, misclass, type = 'l', xlab = "Softmax Threshold", ylab = "Misclassification Rate")
lines(cutoffs, misclass.lb, col = "blue")
legend("topright", legend = c("No confidence intervals", "With confidence intervals"),
  col = c("black", "blue"), lty = 1)
dev.off()



##########################################################
## Table xx: Logistic Example
##########################################################


load("Examples/logistic/out/logistic_real_dim_50.RData")


# Joint region Volume Comparison
c( volm_ibs[2], volm_ebs[2, 2], volm_ebs_ls[2, 2])

# Marginal friendly inferences 

# Max Ratio of lengths of intervals among different dimensions 
ibs_ebs_comprs_max <-  max(ratio_ibs_ebs[2, 2, ])
ibs_ebs_ls_comprs_max <-  max(ratio_ibs_ebs_ls[2, 2, ])
ebs_ls_ebs_comprs_max <-  max(ratio_ebs_ls_ebs[2, 2, ])
c(ibs_ebs_comprs_max, ibs_ebs_ls_comprs_max, ebs_ls_ebs_comprs_max)


# min Ratio of lengths of intervals among different dimensions 
ibs_ebs_comprs_min <-  min(ratio_ibs_ebs[2, 2, ])
ibs_ebs_ls_comprs_min <-  min(ratio_ibs_ebs_ls[2, 2, ])
ebs_ls_ebs_comprs_min <-  min(ratio_ebs_ls_ebs[2, 2, ])
c(ibs_ebs_comprs_min, ibs_ebs_ls_comprs_min, ebs_ls_ebs_comprs_min)


# Mean Ratio of lengths of intervals among different dimensions 
ibs_ebs_comprs_mean <-  mean(ratio_ibs_ebs[2, 2, ])
ibs_ebs_ls_comprs_mean <-  mean(ratio_ibs_ebs_ls[2, 2, ])
ebs_ls_ebs_comprs_mean <-  mean(ratio_ebs_ls_ebs[2, 2, ])
c(ibs_ebs_comprs_mean, ibs_ebs_ls_comprs_mean, ebs_ls_ebs_comprs_mean)

# Marginal friendly region volume of cuboids  inflated with respect to joint 
c(marg_volm_ibs[2], marg_volm_ebs[2, 2], marg_volm_ebs_ls[2, 2])




tmp1 <- margn_up_low[[1]]
len_1 <- (tmp1[, 2] - tmp1[, 1])
indx <- c(1 : 5, 46 : 50)
tmp2 <- tmp1[order(len_1), ]

tmp2[indx, ]#IBS

tmp1 <- margn_up_low[[4]]
len_1 <- (tmp1[, 2] - tmp1[, 1])
indx <- c(1 : 5, 46 : 50)
tmp2 <- tmp1[order(len_1), ]

tmp2[indx, ]#EBS


tmp1 <- margn_up_low[[5]]
len_1 <- (tmp1[, 2] - tmp1[, 1])
indx <- c(1 : 5, 46 : 50)
tmp2 <- tmp1[order(len_1), ]

tmp2[indx, ] #EBS + Lugsail
