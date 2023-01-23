#############
## load data
#############
load("out/ebs_lugsail_bias.RData")

############################
## bias plots
############################

#########################
## plot for alpha = 0.75
#########################
pdf(file = "ebs_bias_51.pdf")
# EBS bias plot against the sample size
plot(sam_siz, out_ebs[,1], type = 'b', pch = 19, col = 'black', ylim = c(min(out_ebs), 0), 
     ylab="Bias upto a constant multiple", xlab="Sample size", 
     cex.main=1.25, cex.lab=1.5, cex.axis=1.1)

# Add a line for Lugsail-EBS bias
lines(sam_siz, out_ebs[,2], type = 'b', pch = 18, col = "blue", lty=2)

# Add legends
legend("bottomright", legend=c("EBS", "Lugsail-EBS"),
       col=c("black", "blue"), lty=1:2, cex= 1)

dev.off()

#########################
## plot for alpha = 0.75
########################
pdf(file = "ebs_bias_75.pdf")
# EBS bias plot against the sample size
plot(sam_siz, out_ebs[,3], type = 'b', pch = 19, col = 'black', ylim = c(min(out_ebs), 0), 
     ylab="Bias upto a constant multiple", xlab="Sample size", 
     cex.main=1.25, cex.lab=1.5, cex.axis=1.1)

# Add a line for Lugsail-EBS bias
lines(sam_siz, out_ebs[,4], type = 'b', pch = 18, col = "blue", lty=2)

# Add legends
legend("topright", legend=c("EBS", "Lugsail-EBS"),
       col=c("black", "blue"), lty=1:2, cex= 1)

dev.off()

