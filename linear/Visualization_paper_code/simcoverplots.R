
##plotting some stuff together
set.seed(573)
source("outputfun.R")
reps = 2000


cover.plot = function(cover.samples,
	reps=2000,overall.conf=.95,colour = c("red4","royalblue1"),leng=.02,...){

	##some required stuff
	cover.est = colMeans(cover.samples)
	cover.cov = cov(cover.samples)/reps
	cover.se = sqrt(diag(cover.cov))

	ints.conf = n.sim.int(cover.cov,overall.conf,cover.est)$ints
	xplace = c(1,1,2,2,3,3)


	##start plotting
	plot(NA,xlim=c(.5,3.5),ylim=c(.5,1),xaxt="n",ylab="Coverage",...)
	axis(1,at=c(1,2,3),labels=c("LB","SI","UB"))
	
	abline(h=c(.8,.9),col=colour,lty=3)

	##plotting the estimates
	points(xplace,cover.est,pch=rep(c(15,16),3),col=colour)


	##plotting the little lines to serve as upper and lower bounds
	segments(xplace-leng,ints.conf[,1],xplace+leng,ints.conf[,1],col=colour)
	segments(xplace-leng,ints.conf[,2],xplace+leng,ints.conf[,2],col=colour)

	##Connecting the upper and lower bounds lines
	segments(x0 = xplace, y0 = ints.conf[,1], 
		x1 = xplace, y1 = ints.conf[,2], col = colour)
}




i500 = read.table("IIDsamples500.txt",header=TRUE)
i1000 = read.table("IIDsamples1k.txt",header=TRUE)
i5000 = read.table("IIDsamples5k.txt",header=TRUE)
i10000 = read.table("IIDsamples10k.txt",header=TRUE)

m2500 = read.table("MCMCsamples2500.txt",header=TRUE)
m5000 = read.table("MCMCsamples5k.txt",header=TRUE)
m10000 = read.table("MCMCsamples10k.txt",header=TRUE)
m25000 = read.table("MCMCsamples25k.txt",header=TRUE)
m50000 = read.table("MCMCsamples50k.txt",header=TRUE)


pdf("simultaneouscoverage.pdf",width=6.4,height=5.5,pointsize=10)
par(mfrow=c(2,4))
par(mar=c(3.5,3,2,.5),mgp=c(2,.5,0),cex.main=1)
cover.plot(i500,overall.conf=.95,xlab="n = 500",main="IID",font.main=1)
cover.plot(i1000,overall.conf=.95,xlab="n = 1000",main="IID",font.main=1)
cover.plot(i5000,overall.conf=.95,xlab="n = 5000",main="IID",font.main=1)
cover.plot(i10000,overall.conf=.95,xlab="n = 10000",main="IID",font.main=1)
##dev.off()

##par(mfrow=c(1,3))
cover.plot(m2500,overall.conf=.95,xlab="n = 2500",main="MCMC",font.main=1)
cover.plot(m5000,overall.conf=.95,xlab="n = 5000",main="MCMC",font.main=1)
cover.plot(m25000,overall.conf=.95,xlab="n = 25000",main="MCMC",font.main=1)
cover.plot(m50000,overall.conf=.95,xlab="n = 50000",main="MCMC",font.main=1)
dev.off()




pdf("simcov_iid_10k.pdf",width=2,height=3.75,pointsize=10)
par(mar=c(3.5,3,2,.5),mgp=c(2,.5,0),cex.main=1,fin=c(2,3.75))
cover.plot(i10000,overall.conf=.95,xlab="n = 10000",main="IID",font.main=1)
dev.off()


pdf("simcov_mcmc_10k.pdf",width=2,height=3.75,pointsize=10)
par(mar=c(3.5,3,2,.5),mgp=c(2,.5,0),cex.main=1,fin=c(2,3.75))
cover.plot(m10000,overall.conf=.95,xlab="n = 10000",main="MCMC",font.main=1)
dev.off()



