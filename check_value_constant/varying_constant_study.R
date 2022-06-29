rm(list=ls())
setwd("C:/Users/Hp/Documents/GitHub/Batch_Means_Online/check_value_constant")
source("ibs_lng.R")
source("ebs_lng.R")
alp <- .51
cns_sq <- c(0.01,0.1,1,2,5)
sq_n <- seq(1e4,1e6,by=5e3)
leng_ebs <- matrix(nrow=5,ncol=length(sq_n))
leng_ibs <- matrix(nrow=5,ncol=length(sq_n))


for( k in 1:length(sq_n)){
  n <- sq_n[k]
  for (m in 1 : length(cns_sq)){
 
  leng_ebs[m,k] <- ebs_lng(n,cns=cns_sq[m],alp,bet_typ=1)
  
  leng_ibs[m,k] <- ibs_lng(n,cns=cns_sq[m],alp)
}
  
}

sampl_size_log10 <- log10(sq_n)
pdf("zero_pt_01_beta_1.pdf")
plot(sampl_size_log10,leng_ebs[1,],type="l",ylim=c(20,2500),xlab="Sample Size ",ylab="Number of batches")
lines(sampl_size_log10,leng_ibs[1,],col="green")
legend( x="topleft",legend=c("EBS", "IBS"),
       col=c("black", "green"), lty=c(1,1), cex=0.8)
dev.off()


pdf("zero_pt_1_beta_1.pdf")
plot(sampl_size_log10,leng_ebs[2,],type="l",ylim=c(20,250),xlab="Sample Size ",ylab="Number of batches")
lines(sampl_size_log10,leng_ibs[2,],col="green")
legend( x="topleft",legend=c("EBS", "IBS"),
        col=c("black", "green"), lty=c(1,1), cex=0.8)
dev.off()

pdf("zero_1_beta_1.pdf")
plot(sampl_size_log10,leng_ebs[3,],type="l",ylim=c(2,30),xlab="Sample Size ",ylab="Number of batches")
lines(sampl_size_log10,leng_ibs[3,],col="green")
legend( x="topleft",legend=c("EBS", "IBS"),
        col=c("black", "green"), lty=c(1,1), cex=0.8)

dev.off()


pdf("zero_2_beta_1.pdf")
plot(sampl_size_log10,leng_ebs[4,],type="l",ylim=c(2,25),xlab="Sample Size ",ylab="Number of batches")
lines(sampl_size_log10,leng_ibs[4,],col="green")
legend( x="topleft",legend=c("EBS", "IBS"),
        col=c("black", "green"), lty=c(1,1), cex=0.8)
dev.off()


pdf("zero_5_beta_1.pdf")
plot(sampl_size_log10,leng_ebs[5,],type="l",ylim=c(1,20),xlab="Sample Size ",ylab="Number of batches")
lines(sampl_size_log10,leng_ibs[5,],col="green")
legend( x="topleft",legend=c("EBS", "IBS"),
        col=c("black", "green"), lty=c(1,1), cex=0.8)
dev.off()
############################################################
#New setup for beta
sq_n <- seq(1e4,1e6,by=5e3)
leng_ebs <- matrix(nrow=5,ncol=length(sq_n))
leng_ibs <- matrix(nrow=5,ncol=length(sq_n))


for( k in 1:length(sq_n)){
  n <- sq_n[k]
  for (m in 1 : length(cns_sq)){
    
    leng_ebs[m,k] <- ebs_lng(n,cns=cns_sq[m],alp,bet_typ=2)
    
    leng_ibs[m,k] <- ibs_lng(n,cns=cns_sq[m],alp)
  }
  
}




sampl_size_log10 <- log10(sq_n)
pdf("zero_pt_01_beta_2.pdf")
plot(sampl_size_log10,leng_ebs[1,],type="l",ylim=c(20,8000),xlab="Sample Size ",ylab="Number of batches")
lines(sampl_size_log10,leng_ibs[1,],col="green")
legend( x="topleft",legend=c("EBS", "IBS"),
        col=c("black", "green"), lty=c(1,1), cex=0.8)
dev.off()


pdf("zero_pt_1_beta_2.pdf")
plot(sampl_size_log10,leng_ebs[2,],type="l",ylim=c(20,900),xlab="Sample Size ",ylab="Number of batches")
lines(sampl_size_log10,leng_ibs[2,],col="green")
legend( x="topleft",legend=c("EBS", "IBS"),
        col=c("black", "green"), lty=c(1,1), cex=0.8)
dev.off()

pdf("zero_1_beta_2.pdf")
plot(sampl_size_log10,leng_ebs[3,],type="l",ylim=c(5,80),xlab="Sample Size ",ylab="Number of batches")
lines(sampl_size_log10,leng_ibs[3,],col="green")
legend( x="topleft",legend=c("EBS", "IBS"),
        col=c("black", "green"), lty=c(1,1), cex=0.8)

dev.off()


pdf("zero_2_beta_2.pdf")
plot(sampl_size_log10,leng_ebs[4,],type="l",ylim=c(2,40),xlab="Sample Size ",ylab="Number of batches")
lines(sampl_size_log10,leng_ibs[4,],col="green")
legend( x="topleft",legend=c("EBS", "IBS"),
        col=c("black", "green"), lty=c(1,1), cex=0.8)
dev.off()


pdf("zero_5_beta_2.pdf")
plot(sampl_size_log10,leng_ebs[5,],type="l",ylim=c(1,20),xlab="Sample Size ",ylab="Number of batches")
lines(sampl_size_log10,leng_ibs[5,],col="green")
legend( x="topleft",legend=c("EBS", "IBS"),
        col=c("black", "green"), lty=c(1,1), cex=0.8)
dev.off()
########################################################

sq_n <- seq(1e3,1e6,by=5e3)
leng_ebs <- matrix(nrow=5,ncol=length(sq_n))
leng_ibs <- matrix(nrow=5,ncol=length(sq_n))


for( k in 1:length(sq_n)){
  n <- sq_n[k]
  for (m in 1 : length(cns_sq)){
    
    leng_ebs[m,k] <- ebs_lng(n,cns=cns_sq[m],alp,bet_typ=3)
    
    leng_ibs[m,k] <- ibs_lng(n,cns=cns_sq[m],alp)
  }
  
}




sampl_size_log10 <- log10(sq_n)
pdf("zero_pt_01_beta_3.pdf")
plot(sampl_size_log10,leng_ebs[1,],type="l",ylim=c(2,9000),xlab="Sample Size ",ylab="Number of batches")
lines(sampl_size_log10,leng_ibs[1,],col="green")
legend( x="topleft",legend=c("EBS", "IBS"),
        col=c("black", "green"), lty=c(1,1), cex=0.8)
dev.off()


pdf("zero_pt_1_beta_3.pdf")
plot(sampl_size_log10,leng_ebs[2,],type="l",ylim=c(20,950),xlab="Sample Size ",ylab="Number of batches")
lines(sampl_size_log10,leng_ibs[2,],col="green")
legend( x="topleft",legend=c("EBS", "IBS"),
        col=c("black", "green"), lty=c(1,1), cex=0.8)
dev.off()

pdf("zero_1_beta_3.pdf")
plot(sampl_size_log10,leng_ebs[3,],type="l",ylim=c(2,100),xlab="Sample Size ",ylab="Number of batches")
lines(sampl_size_log10,leng_ibs[3,],col="green")
legend( x="topleft",legend=c("EBS", "IBS"),
        col=c("black", "green"), lty=c(1,1), cex=0.8)

dev.off()


pdf("zero_2_beta_3.pdf")
plot(sampl_size_log10,leng_ebs[4,],type="l",ylim=c(2,50),xlab="Sample Size ",ylab="Number of batches")
lines(sampl_size_log10,leng_ibs[4,],col="green")
legend( x="topleft",legend=c("EBS", "IBS"),
        col=c("black", "green"), lty=c(1,1), cex=0.8)
dev.off()


pdf("zero_5_beta_3.pdf")
plot(sampl_size_log10,leng_ebs[5,],type="l",ylim=c(1,20),xlab="Sample Size ",ylab="Number of batches")
lines(sampl_size_log10,leng_ibs[5,],col="green")
legend( x="topleft",legend=c("EBS", "IBS"),
        col=c("black", "green"), lty=c(1,1), cex=0.8)
dev.off()



