rm(list = ls())
nparm <- 5
sq_n <- c("5e4","1e5","2e5","5e5","8e5","1e6","5e6")
cover_all <- marg_covg <- forb_norm_all  <-   matrix(nrow = length(sq_n), ncol = 14)
volm_all <- sd_volm <- forb_all <- sd_forb <- matrix(nrow = length(sq_n), ncol = 14)
sd_cover <- sd_marg_covg <- sd_forb_norm <- matrix(nrow = length(sq_n), ncol = 14)
interv_ratio <- sd_interv_ratio <-  array(dim = c(length(sq_n), 12, nparm), dimnames = list(1 : length(sq_n), 1 : 12, 1 : nparm))

load("C:/Users/Hp/Documents/GitHub/Batch_Means_Online/linear/out/linear_indep_n_5e+06_dim_5.RData")

##### COVERAGE RATES
for(n in 1 : length(sq_n)){
  
  cover_all[n, 1] <- mean( cover_orc[,n])
  cover_all[n, 2] <- mean( cover_ibs[,n])
  cover_all[n, 3 : 14] <- c(colMeans(cover_ebs[n, , ]), colMeans(cover_ebs_ls[n, , ]))
  Reps <- length(cover_orc)
  sd_cover[n, 1] <- sd( cover_orc[, n]) / sqrt(Reps)
  sd_cover[n, 2] <- sd( cover_ibs[, n]) / sqrt(Reps)
  sd_cover[n, 3 : 14] <- c(apply(cover_ebs[n, , ], 2, sd), apply(cover_ebs_ls[n, , ], 2, sd)) / sqrt(Reps)
  
  
  
  
}

sq_n <- log10(c(5e4, 1e5, 2e5, 5e5, 8e5,  1e6, 5e6))


col_choic <- c("red", "steelblue", "brown", "tomato", "yellow", "slateblue",  "magenta3", "palegreen" , "brown", "tomato", "yellow", "slateblue",  "magenta3", "palegreen" )
lin_typ <- c(6, 3, rep(1,6), rep(2,6))
names_var <- c("Oracle", "IBS", expression(paste("c = 0.1, ", beta[1])), expression(paste("c = 0.1, ", beta[2])), expression(paste("c = 0.1, ", beta[3])), expression(paste("c = 1, ", beta[1])),  expression(paste("c = 1, ", beta[2])) , expression(paste("c = 1, ", beta[3])), expression(paste("LS c = 0.1, ", beta[1])), expression(paste("LS c = 0.1, ", beta[2])), expression(paste("LS c = 0.1, ", beta[3])), expression(paste("LS c = 1, ", beta[1])),  expression(paste("LS c = 1, ", beta[2])) , expression(paste("LS c = 1, ", beta[3])))


pdf("C:/Users/Hp/Documents/GitHub/Batch_Means_Online/linear/out/cover_rates_linear_indep.pdf",         # File name
    width = 11, height = 12,   # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk")          # Paper size

plot(sq_n, cover_all[, 1], type = "b", lwd = 2, lty = 6, xlab = "Log Sample Size", ylab = "Coverage Rate", col = "red", ylim = c(0.06,1))
arrows(sq_n, cover_all[, 1], sq_n, cover_all[, 1] + sd_cover[,1], length=0.05, angle=90, col = "red", lwd = 2, lty = 1 )
arrows(sq_n, cover_all[, 1], sq_n, cover_all[, 1] - sd_cover[,1], length=0.05, angle=90, col = "red", lwd = 2, lty = 1 )


lines(sq_n, cover_all[, 2],  lwd = 2, col = "steelblue", lty = 3, type = "b")
arrows(sq_n, cover_all[, 2], sq_n, cover_all[, 2] + sd_cover[,2], length=0.05, angle=90, col = "steelblue", lwd = 2, lty = 1 )
arrows(sq_n, cover_all[, 2], sq_n, cover_all[, 2] - sd_cover[,2], length=0.05, angle=90, col = "steelblue", lwd = 2, lty = 1 )



lines(sq_n, cover_all[, 3],  lwd = 2, col = "brown", lty = 1, type = "b")
arrows(sq_n, cover_all[, 3], sq_n, cover_all[, 3] + sd_cover[, 3], length=0.05, angle=90, col = "brown", lwd = 2, lty = 1 )
arrows(sq_n, cover_all[, 3], sq_n, cover_all[, 3] - sd_cover[, 3], length=0.05, angle=90, col = "brown", lwd = 2, lty = 1 )


lines(sq_n, cover_all[, 4],  lwd = 2, col = "tomato", lty = 1, type = "b")
arrows(sq_n, cover_all[, 4], sq_n, cover_all[, 4] + sd_cover[, 4], length=0.05, angle=90, col = "tomato", lwd = 2, lty = 1 )
arrows(sq_n, cover_all[, 4], sq_n, cover_all[, 4] - sd_cover[, 4], length=0.05, angle=90, col = "tomato", lwd = 2, lty = 1 )



lines(sq_n, cover_all[, 5],  lwd = 2, col = "yellow", lty = 1, type = "b")
arrows(sq_n, cover_all[, 5], sq_n, cover_all[, 5] + sd_cover[, 5], length=0.05, angle=90, col = "yellow", lwd = 2, lty = 1 )
arrows(sq_n, cover_all[, 5], sq_n, cover_all[, 5] - sd_cover[, 5], length=0.05, angle=90, col = "yellow", lwd = 2, lty = 1 )


lines(sq_n, cover_all[, 6],  lwd = 2, col = "slateblue", lty = 1, type = "b")
arrows(sq_n, cover_all[, 6], sq_n, cover_all[, 6] + sd_cover[, 6], length=0.05, angle=90, col = "slateblue", lwd = 2, lty = 1 )
arrows(sq_n, cover_all[, 6], sq_n, cover_all[, 6] - sd_cover[, 6], length=0.05, angle=90, col = "slateblue", lwd = 2, lty = 1 )



lines(sq_n, cover_all[, 7],  lwd = 2, col = "magenta3", lty = 1, type = "b")
arrows(sq_n, cover_all[, 7], sq_n, cover_all[, 7] + sd_cover[, 7], length=0.05, angle=90, col = "magenta3", lwd = 2, lty = 1 )
arrows(sq_n, cover_all[, 7], sq_n, cover_all[, 7] - sd_cover[, 7], length=0.05, angle=90, col = "magenta3", lwd = 2, lty = 1 )


lines(sq_n, cover_all[, 8],  lwd = 2, col = "palegreen", lty = 1, type = "b")
arrows(sq_n, cover_all[, 8], sq_n, cover_all[, 8] + sd_cover[, 8], length=0.05, angle=90, col = "palegreen", lwd = 2, lty = 1 )
arrows(sq_n, cover_all[, 8], sq_n, cover_all[, 8] - sd_cover[, 8], length=0.05, angle=90, col = "palegreen", lwd = 2, lty = 1 )


lines(sq_n, cover_all[, 9],  lwd = 2, col = "brown", lty = 2, type = "b")
arrows(sq_n, cover_all[, 9], sq_n, cover_all[, 9] + sd_cover[, 9], length=0.05, angle=90, col = "brown", lwd = 2, lty = 1 )
arrows(sq_n, cover_all[, 9], sq_n, cover_all[, 9] - sd_cover[, 9], length=0.05, angle=90, col = "brown", lwd = 2, lty = 1 )


lines(sq_n, cover_all[, 10],  lwd = 2, col = "tomato", lty = 2, type = "b")
arrows(sq_n, cover_all[, 10], sq_n, cover_all[, 10] + sd_cover[, 10], length=0.05, angle=90, col = "tomato", lwd = 2, lty = 1 )
arrows(sq_n, cover_all[, 10], sq_n, cover_all[, 10] - sd_cover[, 10], length=0.05, angle=90, col = "tomato", lwd = 2, lty = 1 )


lines(sq_n, cover_all[, 11],  lwd = 2, col = "yellow", lty = 2, type = "b")
arrows(sq_n, cover_all[, 11], sq_n, cover_all[, 11] + sd_cover[, 11], length=0.05, angle=90, col = "yellow", lwd = 2, lty = 1 )
arrows(sq_n, cover_all[, 11], sq_n, cover_all[, 11] - sd_cover[, 11], length=0.05, angle=90, col = "yellow", lwd = 2, lty = 1 )



lines(sq_n, cover_all[, 12],  lwd = 2, col = "slateblue", lty = 2, type = "b")
arrows(sq_n, cover_all[, 12], sq_n, cover_all[, 12] + sd_cover[, 12], length=0.05, angle=90, col = "slateblue", lwd = 2, lty = 1 )
arrows(sq_n, cover_all[, 12], sq_n, cover_all[, 12] - sd_cover[, 12], length=0.05, angle=90, col = "slateblue", lwd = 2, lty = 1 )

lines(sq_n, cover_all[, 13],  lwd = 2, col = "magenta3", lty = 2, type = "b")
arrows(sq_n, cover_all[, 13], sq_n, cover_all[, 13] + sd_cover[, 13], length=0.05, angle=90, col = "magenta3", lwd = 2, lty = 1 )
arrows(sq_n, cover_all[, 13], sq_n, cover_all[, 13] - sd_cover[, 13], length=0.05, angle=90, col = "magenta3", lwd = 2, lty = 1 )

lines(sq_n, cover_all[, 14],  lwd = 2, col = "palegreen", lty = 2, type = "b")
arrows(sq_n, cover_all[, 14], sq_n, cover_all[, 14] + sd_cover[, 14], length=0.05, angle=90, col = "palegreen", lwd = 2, lty = 1 )
arrows(sq_n, cover_all[, 14], sq_n, cover_all[, 14] - sd_cover[, 14], length=0.05, angle=90, col = "palegreen", lwd = 2, lty = 1 )


legend("bottom", legend = names_var,
       col = col_choic, lwd = 2, lty = lin_typ, cex=1.2,
       box.lty=1, box.lwd=2, ncol = 5)

dev.off()



##### VOLUME


for(n in 1 : length(sq_n)){
  
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


pdf("C:/Users/Hp/Documents/GitHub/Batch_Means_Online/linear/out/volume_linear_indep.pdf",         # File name
    width = 11, height = 12, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk")          # Paper size

plot(sq_n, volm_all[, 1], type = "b", lwd = 2, lty = 6, xlab = "Log Sample Size", ylab = "Volume", col = "red", ylim = c(0.1,1.35))
arrows(sq_n, volm_all[, 1], sq_n, volm_all[, 1] + sd_volm[,1], length=0.05, angle=90, col = "red", lwd = 2, lty = 1 )
arrows(sq_n, volm_all[, 1], sq_n, volm_all[, 1] - sd_volm[,1], length=0.05, angle=90, col = "red", lwd = 2, lty = 1 )


lines(sq_n, volm_all[, 2],  lwd = 2, col = "steelblue", lty = 3, type = "b")
arrows(sq_n, volm_all[, 2], sq_n, volm_all[, 2] + sd_volm[,2], length=0.05, angle=90, col = "steelblue", lwd = 2, lty = 1 )
arrows(sq_n, volm_all[, 2], sq_n, volm_all[, 2] - sd_volm[,2], length=0.05, angle=90, col = "steelblue", lwd = 2, lty = 1 )


lines(sq_n, volm_all[, 3],  lwd = 2, col = "brown", lty = 1, type = "b")
arrows(sq_n, volm_all[, 3], sq_n, volm_all[, 3] + sd_volm[, 3], length=0.05, angle=90, col = "brown", lwd = 2, lty = 1 )
arrows(sq_n, volm_all[, 3], sq_n, volm_all[, 3] - sd_volm[, 3], length=0.05, angle=90, col = "brown", lwd = 2, lty = 1 )


lines(sq_n, volm_all[, 4],  lwd = 2, col = "tomato", lty = 1, type = "b")
arrows(sq_n, volm_all[, 4], sq_n, volm_all[, 4] + sd_volm[, 4], length=0.05, angle=90, col = "tomato", lwd = 2, lty = 1 )
arrows(sq_n, volm_all[, 4], sq_n, volm_all[, 4] - sd_volm[, 4], length=0.05, angle=90, col = "tomato", lwd = 2, lty = 1 )


lines(sq_n, volm_all[, 5],  lwd = 2, col = "yellow", lty = 1, type = "b")
arrows(sq_n, volm_all[, 5], sq_n, volm_all[, 5] + sd_volm[, 5], length=0.05, angle=90, col = "yellow", lwd = 2, lty = 1 )
arrows(sq_n, volm_all[, 5], sq_n, volm_all[, 5] - sd_volm[, 5], length=0.05, angle=90, col = "yellow", lwd = 2, lty = 1 )


lines(sq_n, volm_all[, 6],  lwd = 2, col = "slateblue", lty = 1, type = "b")
arrows(sq_n, volm_all[, 6], sq_n, volm_all[, 6] + sd_volm[, 6], length=0.05, angle=90, col = "slateblue", lwd = 2, lty = 1 )
arrows(sq_n, volm_all[, 6], sq_n, volm_all[, 6] - sd_volm[, 6], length=0.05, angle=90, col = "slateblue", lwd = 2, lty = 1 )


lines(sq_n, volm_all[, 7],  lwd = 2, col = "magenta3", lty = 1, type = "b")
arrows(sq_n, volm_all[, 7], sq_n, volm_all[, 7] + sd_volm[, 7], length=0.05, angle=90, col = "magenta3", lwd = 2, lty = 1 )
arrows(sq_n, volm_all[, 7], sq_n, volm_all[, 7] - sd_volm[, 7], length=0.05, angle=90, col = "magenta3", lwd = 2, lty = 1 )



lines(sq_n, volm_all[, 8],  lwd = 2, col = "palegreen", lty = 1, type = "b")
arrows(sq_n, volm_all[, 8], sq_n, volm_all[, 8] + sd_volm[, 8], length=0.05, angle=90, col = "palegreen", lwd = 2, lty = 1 )
arrows(sq_n, volm_all[, 8], sq_n, volm_all[, 8] - sd_volm[, 8], length=0.05, angle=90, col = "palegreen", lwd = 2, lty = 1 )



lines(sq_n, volm_all[, 9],  lwd = 2, col = "brown", lty = 2, type = "b")
arrows(sq_n, volm_all[, 9], sq_n, volm_all[, 9] + sd_volm[, 9], length=0.05, angle=90, col = "brown", lwd = 2, lty = 1 )
arrows(sq_n, volm_all[, 9], sq_n, volm_all[, 9] - sd_volm[, 9], length=0.05, angle=90, col = "brown", lwd = 2, lty = 1 )


lines(sq_n, volm_all[, 10],  lwd = 2, col = "tomato", lty = 2, type = "b")
arrows(sq_n, volm_all[, 10], sq_n, volm_all[, 10] + sd_volm[, 10], length=0.05, angle=90, col = "tomato", lwd = 2, lty = 1 )
arrows(sq_n, volm_all[, 10], sq_n, volm_all[, 10] - sd_volm[, 10], length=0.05, angle=90, col = "tomato", lwd = 2, lty = 1 )


lines(sq_n, volm_all[, 11],  lwd = 2, col = "yellow", lty = 2, type = "b")
arrows(sq_n, volm_all[, 11], sq_n, volm_all[, 11] + sd_volm[, 11], length=0.05, angle=90, col = "yellow", lwd = 2, lty = 1 )
arrows(sq_n, volm_all[, 11], sq_n, volm_all[, 11] - sd_volm[, 11], length=0.05, angle=90, col = "yellow", lwd = 2, lty = 1 )


lines(sq_n, volm_all[, 12],  lwd = 2, col = "slateblue", lty = 2, type = "b")
arrows(sq_n, volm_all[, 12], sq_n, volm_all[, 12] + sd_volm[, 12], length=0.05, angle=90, col = "slateblue", lwd = 2, lty = 1 )
arrows(sq_n, volm_all[, 12], sq_n, volm_all[, 12] - sd_volm[, 12], length=0.05, angle=90, col = "slateblue", lwd = 2, lty = 1 )


lines(sq_n, volm_all[, 13],  lwd = 2, col = "magenta3", lty = 2, type = "b")
arrows(sq_n, volm_all[, 13], sq_n, volm_all[, 13] + sd_volm[, 13], length=0.05, angle=90, col = "magenta3", lwd = 2, lty = 1 )
arrows(sq_n, volm_all[, 13], sq_n, volm_all[, 13] - sd_volm[, 13], length=0.05, angle=90, col = "magenta3", lwd = 2, lty = 1 )


lines(sq_n, volm_all[, 14],  lwd = 2, col = "palegreen", lty = 2, type = "b")
arrows(sq_n, volm_all[, 14], sq_n, volm_all[, 14] + sd_volm[, 14], length=0.05, angle=90, col = "palegreen", lwd = 2, lty = 1 )
arrows(sq_n, volm_all[, 14], sq_n, volm_all[, 14] - sd_volm[, 14], length=0.05, angle=90, col = "palegreen", lwd = 2, lty = 1 )




legend("bottom", legend = names_var,
       col = col_choic, lwd = 2, lty = lin_typ, cex=1.2,
       box.lty=1, box.lwd=2, ncol = 5)

dev.off()


#### RELATIVE FROBENIUS NORM

for(n in 1 : length(sq_n)){
  
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


pdf("C:/Users/Hp/Documents/GitHub/Batch_Means_Online/linear/out/relative_frobenius_indep.pdf",         # File name
    width = 11, height = 12, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk")          # Paper size

plot(sq_n, forb_all[, 1], type = "b", lwd = 2, lty = 6, xlab = "Log Sample Size", ylab = "Relative Frobenius", col = "red", ylim = c(-0.1,0.9))
arrows(sq_n, forb_all[, 1], sq_n, forb_all[, 1] + sd_forb[,1], length=0.05, angle=90, col = "red", lwd = 2, lty = 1 )
arrows(sq_n, forb_all[, 1], sq_n, forb_all[, 1] - sd_forb[,1], length=0.05, angle=90, col = "red", lwd = 2, lty = 1 )


lines(sq_n, forb_all[, 2],  lwd = 2, col = "steelblue", lty = 3, type = "b")
arrows(sq_n, forb_all[, 2], sq_n, forb_all[, 2] + sd_forb[,2], length=0.05, angle=90, col = "steelblue", lwd = 2, lty = 1 )
arrows(sq_n, forb_all[, 2], sq_n, forb_all[, 2] - sd_forb[,2], length=0.05, angle=90, col = "steelblue", lwd = 2, lty = 1 )


lines(sq_n, forb_all[, 3],  lwd = 2, col = "brown", lty = 1, type = "b")
arrows(sq_n, forb_all[, 3], sq_n, forb_all[, 3] + sd_forb[, 3], length=0.05, angle=90, col = "brown", lwd = 2, lty = 1 )
arrows(sq_n, forb_all[, 3], sq_n, forb_all[, 3] - sd_forb[, 3], length=0.05, angle=90, col = "brown", lwd = 2, lty = 1 )


lines(sq_n, forb_all[, 4],  lwd = 2, col = "tomato", lty = 1, type = "b")
arrows(sq_n, forb_all[, 4], sq_n, forb_all[, 4] + sd_forb[, 4], length=0.05, angle=90, col = "tomato", lwd = 2, lty = 1 )
arrows(sq_n, forb_all[, 4], sq_n, forb_all[, 4] - sd_forb[, 4], length=0.05, angle=90, col = "tomato", lwd = 2, lty = 1 )


lines(sq_n, forb_all[, 5],  lwd = 2, col = "yellow", lty = 1, type = "b")
arrows(sq_n, forb_all[, 5], sq_n, forb_all[, 5] + sd_forb[, 5], length=0.05, angle=90, col = "yellow", lwd = 2, lty = 1 )
arrows(sq_n, forb_all[, 5], sq_n, forb_all[, 5] - sd_forb[, 5], length=0.05, angle=90, col = "yellow", lwd = 2, lty = 1 )


lines(sq_n, forb_all[, 6],  lwd = 2, col = "slateblue", lty = 1, type = "b")
arrows(sq_n, forb_all[, 6], sq_n, forb_all[, 6] + sd_forb[, 6], length=0.05, angle=90, col = "slateblue", lwd = 2, lty = 1 )
arrows(sq_n, forb_all[, 6], sq_n, forb_all[, 6] - sd_forb[, 6], length=0.05, angle=90, col = "slateblue", lwd = 2, lty = 1 )


lines(sq_n, forb_all[, 7],  lwd = 2, col = "magenta3", lty = 1, type = "b")
arrows(sq_n, forb_all[, 7], sq_n, forb_all[, 7] + sd_forb[, 7], length=0.05, angle=90, col = "magenta3", lwd = 2, lty = 1 )
arrows(sq_n, forb_all[, 7], sq_n, forb_all[, 7] - sd_forb[, 7], length=0.05, angle=90, col = "magenta3", lwd = 2, lty = 1 )


lines(sq_n, forb_all[, 8],  lwd = 2, col = "palegreen", lty = 1, type = "b")
arrows(sq_n, forb_all[, 8], sq_n, forb_all[, 8] + sd_forb[, 8], length=0.05, angle=90, col = "palegreen", lwd = 2, lty = 1 )
arrows(sq_n, forb_all[, 8], sq_n, forb_all[, 8] - sd_forb[, 8], length=0.05, angle=90, col = "palegreen", lwd = 2, lty = 1 )



lines(sq_n, forb_all[, 9],  lwd = 2, col = "brown", lty = 2, type = "b")
arrows(sq_n, forb_all[, 9], sq_n, forb_all[, 9] + sd_forb[, 9], length=0.05, angle=90, col = "brown", lwd = 2, lty = 1 )
arrows(sq_n, forb_all[, 9], sq_n, forb_all[, 9] - sd_forb[, 9], length=0.05, angle=90, col = "brown", lwd = 2, lty = 1 )



lines(sq_n, forb_all[, 10],  lwd = 2, col = "tomato", lty = 2, type = "b")
arrows(sq_n, forb_all[, 10], sq_n, forb_all[, 10] + sd_forb[, 10], length=0.05, angle=90, col = "tomato", lwd = 2, lty = 1 )
arrows(sq_n, forb_all[, 10], sq_n, forb_all[, 10] - sd_forb[, 10], length=0.05, angle=90, col = "tomato", lwd = 2, lty = 1 )


lines(sq_n, forb_all[, 11],  lwd = 2, col = "yellow", lty = 2, type = "b")
arrows(sq_n, forb_all[, 11], sq_n, forb_all[, 11] + sd_forb[, 11], length=0.05, angle=90, col = "yellow", lwd = 2, lty = 1 )
arrows(sq_n, forb_all[, 11], sq_n, forb_all[, 11] - sd_forb[, 11], length=0.05, angle=90, col = "yellow", lwd = 2, lty = 1 )


lines(sq_n, forb_all[, 12],  lwd = 2, col = "slateblue", lty = 2, type = "b")
arrows(sq_n, forb_all[, 12], sq_n, forb_all[, 12] + sd_forb[, 12], length=0.05, angle=90, col = "slateblue", lwd = 2, lty = 1 )
arrows(sq_n, forb_all[, 12], sq_n, forb_all[, 12] - sd_forb[, 12], length=0.05, angle=90, col = "slateblue", lwd = 2, lty = 1 )



lines(sq_n, forb_all[, 13],  lwd = 2, col = "magenta3", lty = 2, type = "b")
arrows(sq_n, forb_all[, 13], sq_n, forb_all[, 13] + sd_forb[, 13], length=0.05, angle=90, col = "magenta3", lwd = 2, lty = 1 )
arrows(sq_n, forb_all[, 13], sq_n, forb_all[, 13] - sd_forb[, 13], length=0.05, angle=90, col = "magenta3", lwd = 2, lty = 1 )


lines(sq_n, forb_all[, 14],  lwd = 2, col = "palegreen", lty = 2, type = "b")
arrows(sq_n, forb_all[, 14], sq_n, forb_all[, 14] + sd_forb[, 14], length=0.05, angle=90, col = "palegreen", lwd = 2, lty = 1 )
arrows(sq_n, forb_all[, 14], sq_n, forb_all[, 14] - sd_forb[, 14], length=0.05, angle=90, col = "palegreen", lwd = 2, lty = 1 )


legend("bottom", legend = names_var,
       col = col_choic, lwd = 2, lty = lin_typ, cex=1.2,
       box.lty=1, box.lwd=2, ncol = 5)

dev.off()





#### SELF FROBENIUS NORM


for(n in 1 : length(sq_n)){
  
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


pdf("C:/Users/Hp/Documents/GitHub/Batch_Means_Online/linear/out/self_frobenius_indep.pdf",         # File name
    width = 11, height = 12, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk")          # Paper size

plot(sq_n, forb_norm_all[, 1], type = "b", lwd = 2, lty = 6, xlab = "Log Sample Size", ylab = "Self Frobenius", col = "red", ylim = c(0.2,4.6))
arrows(sq_n, forb_norm_all[, 1], sq_n, forb_norm_all[, 1] + sd_forb_norm[,1], length=0.05, angle=90, col = "red", lwd = 2, lty = 1 )
arrows(sq_n, forb_norm_all[, 1], sq_n, forb_norm_all[, 1] - sd_forb_norm[,1], length=0.05, angle=90, col = "red", lwd = 2, lty = 1 )


lines(sq_n, forb_norm_all[, 2],  lwd = 2, col = "steelblue", lty = 3, type = "b")
arrows(sq_n, forb_norm_all[, 2], sq_n, forb_norm_all[, 2] + sd_forb_norm[,2], length=0.05, angle=90, col = "steelblue", lwd = 2, lty = 1 )
arrows(sq_n, forb_norm_all[, 2], sq_n, forb_norm_all[, 2] - sd_forb_norm[,2], length=0.05, angle=90, col = "steelblue", lwd = 2, lty = 1 )


lines(sq_n, forb_norm_all[, 3],  lwd = 2, col = "brown", lty = 1, type = "b")
arrows(sq_n, forb_norm_all[, 3], sq_n, forb_norm_all[, 3] + sd_forb_norm[, 3], length=0.05, angle=90, col = "brown", lwd = 2, lty = 1 )
arrows(sq_n, forb_norm_all[, 3], sq_n, forb_norm_all[, 3] - sd_forb_norm[, 3], length=0.05, angle=90, col = "brown", lwd = 2, lty = 1 )


lines(sq_n, forb_norm_all[, 4],  lwd = 2, col = "tomato", lty = 1, type = "b")
arrows(sq_n, forb_norm_all[, 4], sq_n, forb_norm_all[, 4] + sd_forb_norm[, 4], length=0.05, angle=90, col = "tomato", lwd = 2, lty = 1 )
arrows(sq_n, forb_norm_all[, 4], sq_n, forb_norm_all[, 4] - sd_forb_norm[, 4], length=0.05, angle=90, col = "tomato", lwd = 2, lty = 1 )


lines(sq_n, forb_norm_all[, 5],  lwd = 2, col = "yellow", lty = 1, type = "b")
arrows(sq_n, forb_norm_all[, 5], sq_n, forb_norm_all[, 5] + sd_forb_norm[, 5], length=0.05, angle=90, col = "yellow", lwd = 2, lty = 1 )
arrows(sq_n, forb_norm_all[, 5], sq_n, forb_norm_all[, 5] - sd_forb_norm[, 5], length=0.05, angle=90, col = "yellow", lwd = 2, lty = 1 )


lines(sq_n, forb_norm_all[, 6],  lwd = 2, col = "slateblue", lty = 1, type = "b")
arrows(sq_n, forb_norm_all[, 6], sq_n, forb_norm_all[, 6] + sd_forb_norm[, 6], length=0.05, angle=90, col = "slateblue", lwd = 2, lty = 1 )
arrows(sq_n, forb_norm_all[, 6], sq_n, forb_norm_all[, 6] - sd_forb_norm[, 6], length=0.05, angle=90, col = "slateblue", lwd = 2, lty = 1 )


lines(sq_n, forb_norm_all[, 7],  lwd = 2, col = "magenta3", lty = 1, type = "b")
arrows(sq_n, forb_norm_all[, 7], sq_n, forb_norm_all[, 7] + sd_forb_norm[, 7], length=0.05, angle=90, col = "magenta3", lwd = 2, lty = 1 )
arrows(sq_n, forb_norm_all[, 7], sq_n, forb_norm_all[, 7] - sd_forb_norm[, 7], length=0.05, angle=90, col = "magenta3", lwd = 2, lty = 1 )


lines(sq_n, forb_norm_all[, 8],  lwd = 2, col = "palegreen", lty = 1, type = "b")
arrows(sq_n, forb_norm_all[, 8], sq_n, forb_norm_all[, 8] + sd_forb_norm[, 8], length=0.05, angle=90, col = "palegreen", lwd = 2, lty = 1 )
arrows(sq_n, forb_norm_all[, 8], sq_n, forb_norm_all[, 8] - sd_forb_norm[, 8], length=0.05, angle=90, col = "palegreen", lwd = 2, lty = 1 )



lines(sq_n, forb_norm_all[, 9],  lwd = 2, col = "brown", lty = 2, type = "b")
arrows(sq_n, forb_norm_all[, 9], sq_n, forb_norm_all[, 9] + sd_forb_norm[, 9], length=0.05, angle=90, col = "brown", lwd = 2, lty = 1 )
arrows(sq_n, forb_norm_all[, 9], sq_n, forb_norm_all[, 9] - sd_forb_norm[, 9], length=0.05, angle=90, col = "brown", lwd = 2, lty = 1 )


lines(sq_n, forb_norm_all[, 10],  lwd = 2, col = "tomato", lty = 2, type = "b")
arrows(sq_n, forb_norm_all[, 10], sq_n, forb_norm_all[, 10] + sd_forb_norm[, 10], length=0.05, angle=90, col = "tomato", lwd = 2, lty = 1 )
arrows(sq_n, forb_norm_all[, 10], sq_n, forb_norm_all[, 10] - sd_forb_norm[, 10], length=0.05, angle=90, col = "tomato", lwd = 2, lty = 1 )


lines(sq_n, forb_norm_all[, 11],  lwd = 2, col = "yellow", lty = 2, type = "b")
arrows(sq_n, forb_norm_all[, 11], sq_n, forb_norm_all[, 11] + sd_forb_norm[, 11], length=0.05, angle=90, col = "yellow", lwd = 2, lty = 1 )
arrows(sq_n, forb_norm_all[, 11], sq_n, forb_norm_all[, 11] - sd_forb_norm[, 11], length=0.05, angle=90, col = "yellow", lwd = 2, lty = 1 )


lines(sq_n, forb_norm_all[, 12],  lwd = 2, col = "slateblue", lty = 2, type = "b")
arrows(sq_n, forb_norm_all[, 12], sq_n, forb_norm_all[, 12] + sd_forb_norm[, 12], length=0.05, angle=90, col = "slateblue", lwd = 2, lty = 1 )
arrows(sq_n, forb_norm_all[, 12], sq_n, forb_norm_all[, 12] - sd_forb_norm[, 12], length=0.05, angle=90, col = "slateblue", lwd = 2, lty = 1 )


lines(sq_n, forb_norm_all[, 13],  lwd = 2, col = "magenta3", lty = 2, type = "b")
arrows(sq_n, forb_norm_all[, 13], sq_n, forb_norm_all[, 13] + sd_forb_norm[, 13], length=0.05, angle=90, col = "magenta3", lwd = 2, lty = 1 )
arrows(sq_n, forb_norm_all[, 13], sq_n, forb_norm_all[, 13] - sd_forb_norm[, 13], length=0.05, angle=90, col = "magenta3", lwd = 2, lty = 1 )


lines(sq_n, forb_norm_all[, 14],  lwd = 2, col = "palegreen", lty = 2, type = "b")
arrows(sq_n, forb_norm_all[, 14], sq_n, forb_norm_all[, 14] + sd_forb_norm[, 14], length=0.05, angle=90, col = "palegreen", lwd = 2, lty = 1 )
arrows(sq_n, forb_norm_all[, 14], sq_n, forb_norm_all[, 14] - sd_forb_norm[, 14], length=0.05, angle=90, col = "palegreen", lwd = 2, lty = 1 )



legend("bottom", legend = names_var,
       col = col_choic, lwd = 2, lty = lin_typ, cex=1.2,
       box.lty=1, box.lwd=2, ncol = 5)

dev.off()


##### Interval ratios
for(n in 1 : length(sq_n))
{
  for( l in 1 : 6)
  {
    
    for(m in 1 : 5)
    {
      interv_ratio[n, l, m] <- mean( ratio_ibs_ebs[n, ,l ,m])
      interv_ratio[n, 6 + l, m] <- mean( ratio_ibs_ebs_ls[n, ,l ,m])
      #interv_ratio[n, 12 + l, m] <- mean( ratio_ebs_ls_ebs[n, ,l ,m])
      Reps <- length(cover_orc)
      sd_interv_ratio[n, l, m] <- sd( ratio_ibs_ebs[n, ,l ,m]) / sqrt(Reps)
      sd_interv_ratio[n, 6 + l, m] <- sd( ratio_ibs_ebs_ls[n, ,l ,m]) / sqrt(Reps)
     # sd_interv_ratio[n, 12 + l, m] <- sd( ratio_ebs_ls_ebs[n, ,l ,m]) / sqrt(Reps)
      
    }
    
  }
}

sq_n <- log10(c(5e4, 1e5, 2e5, 5e5, 8e5,  1e6, 5e6))


col_choic <- c("brown", "tomato", "yellow", "slateblue",  "magenta3", "palegreen" , "brown", "tomato", "yellow", "slateblue",  "magenta3", "palegreen")
lin_typ <- c( rep(1, 6), rep(2, 6))
names_var <- c( expression(paste("c = 0.1, ", beta[1])), expression(paste("c = 0.1, ", beta[2])), expression(paste("c = 0.1, ", beta[3])), expression(paste("c = 1, ", beta[1])),  expression(paste("c = 1, ", beta[2])) , expression(paste("c = 1, ", beta[3])), expression(paste("LS c = 0.1, ", beta[1])), expression(paste("LS c = 0.1, ", beta[2])), expression(paste("LS c = 0.1, ", beta[3])), expression(paste("LS c = 1, ", beta[1])),  expression(paste("LS c = 1, ", beta[2])) , expression(paste("LS c = 1, ", beta[3])))


pdf("C:/Users/Hp/Documents/GitHub/Batch_Means_Online/linear/out/interval_ratio_first_comp_linear_indep.pdf",         # File name
    width = 11, height = 12,   # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk")          # Paper size



plot(sq_n, interv_ratio[,  1, 1],  lwd = 2, col = "brown", lty = 1, type = "b", ylim = c(0.6, 1.8),  xlab = "Log Sample Size", ylab = "Interval Ratio First Comp")
arrows(sq_n, interv_ratio[,  1, 1], sq_n, interv_ratio[,  1, 1] + sd_interv_ratio[,  1, 1], length=0.05, angle=90, col = "brown", lwd = 2, lty = 1 )
arrows(sq_n, interv_ratio[,  1, 1], sq_n, interv_ratio[,  1, 1] - sd_interv_ratio[,  1, 1], length=0.05, angle=90, col = "brown", lwd = 2, lty = 1 )


lines(sq_n, interv_ratio[,  2, 1],  lwd = 2, col = "tomato", lty = 1, type = "b")
arrows(sq_n, interv_ratio[,  2, 1], sq_n, interv_ratio[,  2, 1] + sd_interv_ratio[,  2, 1], length=0.05, angle=90, col = "tomato", lwd = 2, lty = 1 )
arrows(sq_n, interv_ratio[,  2, 1], sq_n, interv_ratio[,  2, 1] - sd_interv_ratio[,  2, 1], length=0.05, angle=90, col = "tomato", lwd = 2, lty = 1 )



lines(sq_n, interv_ratio[,  3, 1],  lwd = 2, col = "yellow", lty = 1, type = "b")
arrows(sq_n, interv_ratio[,  3, 1], sq_n, interv_ratio[,  3, 1] + sd_interv_ratio[,  3, 1], length=0.05, angle=90, col = "yellow", lwd = 2, lty = 1 )
arrows(sq_n, interv_ratio[,  3, 1], sq_n, interv_ratio[,  3, 1] - sd_interv_ratio[,  3, 1], length=0.05, angle=90, col = "yellow", lwd = 2, lty = 1 )


lines(sq_n, interv_ratio[,  4, 1],  lwd = 2, col = "slateblue", lty = 1, type = "b")
arrows(sq_n, interv_ratio[,  4, 1], sq_n, interv_ratio[,  4, 1] + sd_interv_ratio[,  4, 1], length=0.05, angle=90, col = "slateblue", lwd = 2, lty = 1 )
arrows(sq_n, interv_ratio[,  4, 1], sq_n, interv_ratio[,  4, 1] - sd_interv_ratio[,  4, 1], length=0.05, angle=90, col = "slateblue", lwd = 2, lty = 1 )



lines(sq_n, interv_ratio[,  5, 1],  lwd = 2, col = "magenta3", lty = 1, type = "b")
arrows(sq_n, interv_ratio[,  5, 1], sq_n, interv_ratio[,  5, 1] + sd_interv_ratio[,  5, 1], length=0.05, angle=90, col = "magenta3", lwd = 2, lty = 1 )
arrows(sq_n, interv_ratio[,  5, 1], sq_n, interv_ratio[,  5, 1] - sd_interv_ratio[,  5, 1], length=0.05, angle=90, col = "magenta3", lwd = 2, lty = 1 )


lines(sq_n, interv_ratio[,  6, 1],  lwd = 2, col = "palegreen", lty = 1, type = "b")
arrows(sq_n, interv_ratio[,  6, 1], sq_n, interv_ratio[,  6, 1] + sd_interv_ratio[,  6, 1], length=0.05, angle=90, col = "palegreen", lwd = 2, lty = 1 )
arrows(sq_n, interv_ratio[,  6, 1], sq_n, interv_ratio[,  6, 1] - sd_interv_ratio[,  6, 1], length=0.05, angle=90, col = "palegreen", lwd = 2, lty = 1 )


lines(sq_n, interv_ratio[,  7, 1],  lwd = 2, col = "brown", lty = 2, type = "b")
arrows(sq_n, interv_ratio[,  7, 1], sq_n, interv_ratio[,  7, 1] + sd_interv_ratio[,  7, 1], length=0.05, angle=90, col = "brown", lwd = 2, lty = 1 )
arrows(sq_n, interv_ratio[,  7, 1], sq_n, interv_ratio[,  7, 1] - sd_interv_ratio[,  7, 1], length=0.05, angle=90, col = "brown", lwd = 2, lty = 1 )


lines(sq_n, interv_ratio[,  8, 1],  lwd = 2, col = "tomato", lty = 2, type = "b")
arrows(sq_n, interv_ratio[,  8, 1], sq_n, interv_ratio[,  8, 1] + sd_interv_ratio[,  8, 1], length=0.05, angle=90, col = "tomato", lwd = 2, lty = 1 )
arrows(sq_n, interv_ratio[,  8, 1], sq_n, interv_ratio[,  8, 1] - sd_interv_ratio[,  8, 1], length=0.05, angle=90, col = "tomato", lwd = 2, lty = 1 )


lines(sq_n, interv_ratio[, 9, 1],  lwd = 2, col = "yellow", lty = 2, type = "b")
arrows(sq_n, interv_ratio[, 9, 1], sq_n, interv_ratio[, 9, 1] + sd_interv_ratio[, 9, 1], length=0.05, angle=90, col = "yellow", lwd = 2, lty = 1 )
arrows(sq_n, interv_ratio[, 9, 1], sq_n, interv_ratio[, 9, 1] - sd_interv_ratio[, 9, 1], length=0.05, angle=90, col = "yellow", lwd = 2, lty = 1 )



lines(sq_n, interv_ratio[, 10, 1],  lwd = 2, col = "slateblue", lty = 2, type = "b")
arrows(sq_n, interv_ratio[, 10, 1], sq_n, interv_ratio[, 10, 1] + sd_interv_ratio[, 10, 1], length=0.05, angle=90, col = "slateblue", lwd = 2, lty = 1 )
arrows(sq_n, interv_ratio[, 10, 1], sq_n, interv_ratio[, 10, 1] - sd_interv_ratio[, 10, 1], length=0.05, angle=90, col = "slateblue", lwd = 2, lty = 1 )



lines(sq_n, interv_ratio[,  11, 1],  lwd = 2, col = "magenta3", lty = 2, type = "b")
arrows(sq_n, interv_ratio[,  11, 1], sq_n, interv_ratio[,  11, 1] + sd_interv_ratio[,  11, 1], length=0.05, angle=90, col = "magenta3", lwd = 2, lty = 1 )
arrows(sq_n, interv_ratio[,  11, 1], sq_n, interv_ratio[,  11, 1] - sd_interv_ratio[,  11, 1], length=0.05, angle=90, col = "magenta3", lwd = 2, lty = 1 )

lines(sq_n, interv_ratio[,  12, 1],  lwd = 2, col = "palegreen", lty = 2, type = "b")
arrows(sq_n, interv_ratio[,  12, 1], sq_n, interv_ratio[,  12, 1] + sd_interv_ratio[,  12, 1], length=0.05, angle=90, col = "palegreen", lwd = 2, lty = 1 )
arrows(sq_n, interv_ratio[,  12, 1], sq_n, interv_ratio[,  12, 1] - sd_interv_ratio[,  12, 1], length=0.05, angle=90, col = "palegreen", lwd = 2, lty = 1 )


legend("bottom", legend = names_var,
       col = col_choic, lwd = 2, lty = lin_typ, cex=1.2,
       box.lty=1, box.lwd=2, ncol = 5)

dev.off()
