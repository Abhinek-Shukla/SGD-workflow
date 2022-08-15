###############################################################################
##generating plots to showcase upper and lower bounds
set.seed(573)
library(mvtnorm)
library(ellipse)


alpha = .1
conf = 1 - alpha
p = 2

##p: perfect correlation
##h: rho = .75
##i: independent
Sigma.pcor = matrix(c(9,6,6,4),nrow=2)
Sigma.half = matrix(c(9,4.5,4.5,4),nrow=2)
Sigma.ind = matrix(c(9,0,0,4),nrow=2)
s.x = 3
s.y = 2


n = 100000
X.p = rmvnorm(n,c(0,0),Sigma.pcor)
X.h = rmvnorm(n,c(0,0),Sigma.half)
X.i = rmvnorm(n,c(0,0),Sigma.ind)

##m: marginal
##b: bonferroni
z.m = qnorm(1-alpha/2)
z.b = qnorm(1 - (alpha/2)/2)

hwu.1 = z.m*s.x
hwo.1 = z.b*s.x
hwu.2 = z.m*s.y
hwo.2 = z.b*s.y

checks = 1000

ratio.12 = hwu.1/hwu.2
xs = seq(hwu.1,hwo.1,length=checks)
ys = xs/ratio.12

bl.u = c(-hwu.1,-hwu.2)
tl.u = c(-hwu.1,hwu.2)
br.u = c(hwu.1,-hwu.2)
tr.u = c(hwu.1,hwu.2)

bl.o = c(-hwo.1,-hwo.2)
tl.o = c(-hwo.1,hwo.2)
br.o = c(hwo.1,-hwo.2)
tr.o = c(hwo.1,hwo.2)

ends.u = rbind(bl.u,br.u,tl.u,tr.u)
ends.o = rbind(bl.o,br.o,tl.o,tr.o)


pdf("figsimconf2.pdf",width=4,height=3.25,pointsize=10)

par(mar=c(3,3,1.5,.5),mgp=c(1.5,.5,0),ps=10,cex.main=1,fin=c(4,3.25))
plot(ends.o,col="red4",pch=16, xlab="Component 1", ylab="Component 2",
	main="Bounds of Simultaneous Confidence Intervals",font.main=1,
	xlim=c(-7,7), ylim=c(-5,5))
points(ends.u, col="blue2",pch=16)
lines(xs,ys, col="black")

lines(rbind(bl.u,br.u),col="blue", lty=c(9))
lines(rbind(bl.u,tl.u),col="blue", lty=c(9))
lines(rbind(br.u,tr.u),col="blue", lty=c(9))
lines(rbind(tl.u,tr.u),col="blue", lty=c(9))


lines(rbind(bl.o,br.o),col="red", lty=c(2))
lines(rbind(bl.o,tl.o),col="red", lty=c(2))
lines(rbind(br.o,tr.o),col="red", lty=c(2))
lines(rbind(tl.o,tr.o),col="red", lty=c(2))


legend("center", legend=c("Lower Bound", "Upper Bound", "Potential Halfwidths",
	"90% Ellipse"),
      col=c("blue", "red", 1, "magenta"),
	lty=c(9,2,1,1), pch = c(16,16,NA,NA), cex=1)


##Sigma.ind, Sigma.pcor, Sigma.half
lines(ellipse(Sigma.half, level = .9), lty=1, col="magenta")
dev.off()

