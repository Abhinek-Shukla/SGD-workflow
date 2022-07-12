rm(list = ls())
sq_n <- c("5e4","1e5","2e5","5e5","8e5","1e6","5e6","1e7")
cover_all <- forb_norm_all  <- matrix(nrow = length(sq_n), ncol = 14)
volm_all <- sd_volm <- matrix(nrow = length(sq_n), ncol = 14)
forb_all <- sd_forb <- matrix(nrow = length(sq_n), ncol = 14)
sd_cover <- sd_forb_norm <- matrix(nrow = length(sq_n), ncol = 14)

load("C:/Users/Hp/Documents/GitHub/Batch_Means_Online/linear/out/linear_indep_n_1e+07_dim_5.RData")

##### COVERAGE RATES
for(n in 1:length(sq_n)){

  cover_all[n, 1] <- mean( cover_orc[,n])
  cover_all[n, 2] <- mean( cover_ibs[,n])
  cover_all[n, 3 : 14] <- c(colMeans(cover_ebs[n, , ]), colMeans(cover_ebs_ls[n, , ]))
  Reps <- length(cover_orc)
  sd_cover[n, 1] <- sd( cover_orc[, n]) / sqrt(Reps)
  sd_cover[n, 2] <- sd( cover_ibs[, n]) / sqrt(Reps)
  sd_cover[n, 3 : 14] <- c(apply(cover_ebs[n, , ], 2, sd), apply(cover_ebs_ls[n, , ], 2, sd)) / sqrt(Reps)
  
  
  
  
}

sq_n <- log10(c(5e4, 1e5, 2e5, 5e5, 8e5,  1e6, 5e6, 1e7))


col_choic <- c("red", "steelblue", "brown", "tomato", "yellow", "slateblue",  "magenta3", "palegreen" , "brown", "tomato", "yellow", "slateblue",  "magenta3", "palegreen" )
lin_typ <- c(6, 3, rep(1,6), rep(2,6))
names_var <- c("Oracle", "IBS", expression(paste("c = 0.1, ", beta[1])), expression(paste("c = 0.1, ", beta[2])), expression(paste("c = 0.1, ", beta[3])), expression(paste("c = 1, ", beta[1])),  expression(paste("c = 1, ", beta[2])) , expression(paste("c = 1, ", beta[3])), expression(paste("LS c = 0.1, ", beta[1])), expression(paste("LS c = 0.1, ", beta[2])), expression(paste("LS c = 0.1, ", beta[3])), expression(paste("LS c = 1, ", beta[1])),  expression(paste("LS c = 1, ", beta[2])) , expression(paste("LS c = 1, ", beta[3])))


pdf("C:/Users/Hp/Documents/GitHub/Batch_Means_Online/linear/out/cover_rates_linear.pdf",         # File name
    width = 11, height = 12, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk")          # Paper size

plot(sq_n, cover_all[, 1], type = "b", lwd = 2, lty = 6, xlab = "Log Sample Size", ylab = "Coverage Rate", col = "red", ylim = c(0.2,1))
lines(sq_n, cover_all[, 2],  lwd = 2, col = "steelblue", lty = 3, type = "b")
lines(sq_n, cover_all[, 3],  lwd = 2, col = "brown", lty = 1, type = "b")
lines(sq_n, cover_all[, 4],  lwd = 2, col = "tomato", lty = 1, type = "b")
lines(sq_n, cover_all[, 5],  lwd = 2, col = "yellow", lty = 1, type = "b")
lines(sq_n, cover_all[, 6],  lwd = 2, col = "slateblue", lty = 1, type = "b")
lines(sq_n, cover_all[, 7],  lwd = 2, col = "magenta3", lty = 1, type = "b")
lines(sq_n, cover_all[, 8],  lwd = 2, col = "palegreen", lty = 1, type = "b")

lines(sq_n, cover_all[, 9],  lwd = 2, col = "brown", lty = 2, type = "b")
lines(sq_n, cover_all[, 10],  lwd = 2, col = "tomato", lty = 2, type = "b")
lines(sq_n, cover_all[, 11],  lwd = 2, col = "yellow", lty = 2, type = "b")
lines(sq_n, cover_all[, 12],  lwd = 2, col = "slateblue", lty = 2, type = "b")
lines(sq_n, cover_all[, 13],  lwd = 2, col = "magenta3", lty = 2, type = "b")
lines(sq_n, cover_all[, 14],  lwd = 2, col = "palegreen", lty = 2, type = "b")

legend("bottom", legend = names_var,
       col = col_choic, lwd = 2, lty = lin_typ, cex=1.2,
       box.lty=1, box.lwd=2, ncol = 5)

dev.off()



##### VOLUME


for(n in 1:length(sq_n)){
  
  volm_all[n, 1] <- 1
  volm_all[n, 2] <- mean( volm_ibs[, n])
  volm_all[n, 3 : 14] <- c(colMeans(volm_ebs[n, , ]), colMeans(volm_ebs_ls[n, , ]))
  
  sd_volm[n, 1] <- 0
  sd_volm[n, 2] <- sd( volm_ibs[, n]) / sqrt(Reps)
  sd_volm[n, 3 : 14] <- c(apply(volm_ebs[n, , ], 2, sd), apply(volm_ebs_ls[n, , ], 2, sd)) / sqrt(Reps)
  
  
  
  
}




col_choic <- c("red", "steelblue", "brown", "tomato", "yellow", "slateblue",  "magenta3", "palegreen" , "brown", "tomato", "yellow", "slateblue",  "magenta3", "palegreen" )
lin_typ <- c(6, 3, rep(1,6), rep(2,6))
names_var <- c("Oracle", "IBS", expression(paste("c = 0.1, ", beta[1])), expression(paste("c = 0.1, ", beta[2])), expression(paste("c = 0.1, ", beta[3])), expression(paste("c = 1, ", beta[1])),  expression(paste("c = 1, ", beta[2])) , expression(paste("c = 1, ", beta[3])), expression(paste("LS c = 0.1, ", beta[1])), expression(paste("LS c = 0.1, ", beta[2])), expression(paste("LS c = 0.1, ", beta[3])), expression(paste("LS c = 1, ", beta[1])),  expression(paste("LS c = 1, ", beta[2])) , expression(paste("LS c = 1, ", beta[3])))


pdf("C:/Users/Hp/Documents/GitHub/Batch_Means_Online/linear/out/volume_linear.pdf",         # File name
    width = 11, height = 12, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk")          # Paper size

plot(sq_n, volm_all[, 1], type = "b", lwd = 2, lty = 6, xlab = "Log Sample Size", ylab = "Coverage Rate", col = "red", ylim = c(0.1,1))
lines(sq_n, volm_all[, 2],  lwd = 2, col = "steelblue", lty = 3, type = "b")
lines(sq_n, volm_all[, 3],  lwd = 2, col = "brown", lty = 1, type = "b")
lines(sq_n, volm_all[, 4],  lwd = 2, col = "tomato", lty = 1, type = "b")
lines(sq_n, volm_all[, 5],  lwd = 2, col = "yellow", lty = 1, type = "b")
lines(sq_n, volm_all[, 6],  lwd = 2, col = "slateblue", lty = 1, type = "b")
lines(sq_n, volm_all[, 7],  lwd = 2, col = "magenta3", lty = 1, type = "b")
lines(sq_n, volm_all[, 8],  lwd = 2, col = "palegreen", lty = 1, type = "b")

lines(sq_n, volm_all[, 9],  lwd = 2, col = "brown", lty = 2, type = "b")
lines(sq_n, volm_all[, 10],  lwd = 2, col = "tomato", lty = 2, type = "b")
lines(sq_n, volm_all[, 11],  lwd = 2, col = "yellow", lty = 2, type = "b")
lines(sq_n, volm_all[, 12],  lwd = 2, col = "slateblue", lty = 2, type = "b")
lines(sq_n, volm_all[, 13],  lwd = 2, col = "magenta3", lty = 2, type = "b")
lines(sq_n, volm_all[, 14],  lwd = 2, col = "palegreen", lty = 2, type = "b")

legend("bottom", legend = names_var,
       col = col_choic, lwd = 2, lty = lin_typ, cex=1.2,
       box.lty=1, box.lwd=2, ncol = 5)

dev.off()


#### RELATIVE FROBENIUS NORM

for(n in 1:length(sq_n)){
  
  forb_all[n, 1] <- 0
  forb_all[n, 2] <- mean( forb_ibs[, n])
  forb_all[n, 3 : 14] <- c(colMeans(forb_ebs[n, , ]), colMeans(forb_ebs_ls[n, , ]))
  
  sd_forb[n, 1] <- 0
  sd_forb[n, 2] <- sd( forb_ibs[, n]) / sqrt(Reps)
  sd_forb[n, 3 : 14] <- c(apply(forb_ebs[n, , ], 2, sd), apply(forb_ebs_ls[n, , ], 2, sd)) / sqrt(Reps)
  
  
  
  
}




col_choic <- c("red", "steelblue", "brown", "tomato", "yellow", "slateblue",  "magenta3", "palegreen" , "brown", "tomato", "yellow", "slateblue",  "magenta3", "palegreen" )
lin_typ <- c(6, 3, rep(1,6), rep(2,6))
names_var <- c("Oracle", "IBS", expression(paste("c = 0.1, ", beta[1])), expression(paste("c = 0.1, ", beta[2])), expression(paste("c = 0.1, ", beta[3])), expression(paste("c = 1, ", beta[1])),  expression(paste("c = 1, ", beta[2])) , expression(paste("c = 1, ", beta[3])), expression(paste("LS c = 0.1, ", beta[1])), expression(paste("LS c = 0.1, ", beta[2])), expression(paste("LS c = 0.1, ", beta[3])), expression(paste("LS c = 1, ", beta[1])),  expression(paste("LS c = 1, ", beta[2])) , expression(paste("LS c = 1, ", beta[3])))


pdf("C:/Users/Hp/Documents/GitHub/Batch_Means_Online/linear/out/relative_frobenius.pdf",         # File name
    width = 11, height = 12, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk")          # Paper size

plot(sq_n, forb_all[, 1], type = "b", lwd = 2, lty = 6, xlab = "Log Sample Size", ylab = "Coverage Rate", col = "red", ylim = c(-0.1,0.8))
lines(sq_n, forb_all[, 2],  lwd = 2, col = "steelblue", lty = 3, type = "b")
lines(sq_n, forb_all[, 3],  lwd = 2, col = "brown", lty = 1, type = "b")
lines(sq_n, forb_all[, 4],  lwd = 2, col = "tomato", lty = 1, type = "b")
lines(sq_n, forb_all[, 5],  lwd = 2, col = "yellow", lty = 1, type = "b")
lines(sq_n, forb_all[, 6],  lwd = 2, col = "slateblue", lty = 1, type = "b")
lines(sq_n, forb_all[, 7],  lwd = 2, col = "magenta3", lty = 1, type = "b")
lines(sq_n, forb_all[, 8],  lwd = 2, col = "palegreen", lty = 1, type = "b")

lines(sq_n, forb_all[, 9],  lwd = 2, col = "brown", lty = 2, type = "b")
lines(sq_n, forb_all[, 10],  lwd = 2, col = "tomato", lty = 2, type = "b")
lines(sq_n, forb_all[, 11],  lwd = 2, col = "yellow", lty = 2, type = "b")
lines(sq_n, forb_all[, 12],  lwd = 2, col = "slateblue", lty = 2, type = "b")
lines(sq_n, forb_all[, 13],  lwd = 2, col = "magenta3", lty = 2, type = "b")
lines(sq_n, forb_all[, 14],  lwd = 2, col = "palegreen", lty = 2, type = "b")

legend("bottom", legend = names_var,
       col = col_choic, lwd = 2, lty = lin_typ, cex=1.2,
       box.lty=1, box.lwd=2, ncol = 5)

dev.off()





#### SELF FROBENIUS NORM


for(n in 1:length(sq_n)){
  
  forb_norm_all[n, 1] <- sqrt(5)
  forb_norm_all[n, 2] <- mean( forb_ibs_norm[, n])
  forb_norm_all[n, 3 : 14] <- c(colMeans(forb_ebs_norm[n, , ]), colMeans(forb_ebs_norm_ls[n, , ]))
  
  sd_forb_norm[n, 1] <- 0
  sd_forb_norm[n, 2] <- sd( forb_ibs_norm[, n]) / sqrt(Reps)
  sd_forb_norm[n, 3 : 14] <- c(apply(forb_ebs_norm[n, , ], 2, sd), apply(forb_ebs_norm_ls[n, , ], 2, sd)) / sqrt(Reps)
  
  
  
  
}




col_choic <- c("red", "steelblue", "brown", "tomato", "yellow", "slateblue",  "magenta3", "palegreen" , "brown", "tomato", "yellow", "slateblue",  "magenta3", "palegreen" )
lin_typ <- c(6, 3, rep(1,6), rep(2,6))
names_var <- c("Oracle", "IBS", expression(paste("c = 0.1, ", beta[1])), expression(paste("c = 0.1, ", beta[2])), expression(paste("c = 0.1, ", beta[3])), expression(paste("c = 1, ", beta[1])),  expression(paste("c = 1, ", beta[2])) , expression(paste("c = 1, ", beta[3])), expression(paste("LS c = 0.1, ", beta[1])), expression(paste("LS c = 0.1, ", beta[2])), expression(paste("LS c = 0.1, ", beta[3])), expression(paste("LS c = 1, ", beta[1])),  expression(paste("LS c = 1, ", beta[2])) , expression(paste("LS c = 1, ", beta[3])))


pdf("C:/Users/Hp/Documents/GitHub/Batch_Means_Online/linear/out/self_frobenius.pdf",         # File name
    width = 11, height = 12, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk")          # Paper size

plot(sq_n, forb_norm_all[, 1], type = "b", lwd = 2, lty = 6, xlab = "Log Sample Size", ylab = "Coverage Rate", col = "red", ylim = c(0.3,2.7))
lines(sq_n, forb_norm_all[, 2],  lwd = 2, col = "steelblue", lty = 3, type = "b")
lines(sq_n, forb_norm_all[, 3],  lwd = 2, col = "brown", lty = 1, type = "b")
lines(sq_n, forb_norm_all[, 4],  lwd = 2, col = "tomato", lty = 1, type = "b")
lines(sq_n, forb_norm_all[, 5],  lwd = 2, col = "yellow", lty = 1, type = "b")
lines(sq_n, forb_norm_all[, 6],  lwd = 2, col = "slateblue", lty = 1, type = "b")
lines(sq_n, forb_norm_all[, 7],  lwd = 2, col = "magenta3", lty = 1, type = "b")
lines(sq_n, forb_norm_all[, 8],  lwd = 2, col = "palegreen", lty = 1, type = "b")

lines(sq_n, forb_norm_all[, 9],  lwd = 2, col = "brown", lty = 2, type = "b")
lines(sq_n, forb_norm_all[, 10],  lwd = 2, col = "tomato", lty = 2, type = "b")
lines(sq_n, forb_norm_all[, 11],  lwd = 2, col = "yellow", lty = 2, type = "b")
lines(sq_n, forb_norm_all[, 12],  lwd = 2, col = "slateblue", lty = 2, type = "b")
lines(sq_n, forb_norm_all[, 13],  lwd = 2, col = "magenta3", lty = 2, type = "b")
lines(sq_n, forb_norm_all[, 14],  lwd = 2, col = "palegreen", lty = 2, type = "b")

legend("bottom", legend = names_var,
       col = col_choic, lwd = 2, lty = lin_typ, cex=1.2,
       box.lty=1, box.lwd=2, ncol = 5)

dev.off()

