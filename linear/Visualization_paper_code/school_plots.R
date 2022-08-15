
################################################################################
################################################################################
##code to generate plots used in Gelman's School Data example

##update: Aug 7, 2020
##march 30,2019
##generate our plots for the gelman school example
rm(list=ls())
set.seed(573)
source("outputfun.R")

plot.sim.up = function(x,ests,ints,dens.n=1000,se.den.n=100,ran.range=NULL,
		xlab=NULL,ylab=NULL,ylim=NULL,est.lwd=2,...){
	p = length(ests)
	if(is.null(ran.range)){ran = seq( min(ints)-1, max(ints)+1,
		length=dens.n)}else{ran = seq(min(ran.range),
			max(ran.range),length=dens.n)}
	den.ests = density(x, from=min(ran), to=max(ran), n=length(ran))$y
	if(is.null(xlab)){xlab="x"}
	if(is.null(ylab)){ylab="f(x)"}
	if(is.null(ylim)){ylim=c(0,max(den.ests))}
	plot(NULL,type="l",xlab=xlab,xlim=c(min(ran),max(ran)),
		ylab=ylab,ylim=ylim,yaxs="i",...)
	for(i in 1:p){
		xs = c(ints[i,1],seq(ints[i,1],ints[i,2],length=se.den.n),
			ints[i,2])
		ys=c(0,density(x, from=xs[2], to = xs[length(xs)-1],
			n = (length(xs)-2))$y,0)
		polygon(xs,ys,col="powderblue",border=NA)
	}
	abline(h=0)
	lines(ran,den.ests)
	den.of.est = numeric(p)
	for(i in 1:p){
		den.of.est[i] = density(x, from=ests[i],to=ests[i],n=1)$y
	}
	segments(ests, rep(0,p), ests, den.of.est, col="black", lwd=est.lwd)
}


in.cred.ible = function(X,ests,ints,LWD=2,...){

	p = length(ests)
	p.1 = dim(X)[2]

	low.80 = ests[(p.1+1):(2*p.1)]
	upp.80 = ests[(2*p.1+1):(3*p.1)]
	low.95 = ests[(3*p.1+1):(4*p.1)]
	upp.95 = ests[(4*p.1+1):(5*p.1)]

	xrange = c(min(ints)-1, max(ints)+1)

	index = rep(1:p.1,5)
	qw = .1

	par(bg="white")
	ylabels = sapply(1:p.1, function(i){bquote(theta[.(i)])})

	plot(NULL, xlim=xrange,ylim=c(1,p.1),yaxt="n",xlab="",ylab="",...)
	axis(2, at=1:p.1, labels=c(expression(theta[1]),ylabels[-1]), las=2)

	segments(low.95,index,low.80,index, lwd=LWD)
	segments(upp.80,index,upp.95,index, lwd=LWD)

	for(i in 1:p.1){
		obx = c(low.80[i],low.80[i],upp.80[i],upp.80[i])
		oby = c(index[i]-qw,index[i]+qw,index[i]+qw,index[i]-qw)
		polygon(obx,oby,col="grey90",border=NA)
	}

	cols = c(rep("dodgerblue",p.1),rep("brown1",2*p.1),rep("coral",2*p.1))
	for(i in 1:p){
		xs = c(ints[i,1],ints[i,1],ints[i,2],ints[i,2])
		ys=c(index[i]-qw,index[i]+qw,index[i]+qw,index[i]-qw)
		polygon(xs,ys,col=cols[i],border=NA)
	}

	segments(ests,index-qw,ests,index+qw,col="black", lwd=LWD)

}
##############################################################################

##########################################################################
##sampler functions
theta.update <- function (mu,tau){
	theta.hat <- (mu/tau^2 + y/sigma.y^2)/(1/tau^2 + 1/sigma.y^2)
	V.theta <- 1/(1/tau^2 + 1/sigma.y^2)
	rnorm (8, theta.hat, sqrt(V.theta))
}
mu.update <- function (theta,tau){
	rnorm (1, mean(theta), tau/sqrt(8))
}
tau.update <- function (theta,mu){
	sqrt(sum((theta-mu)^2)/rchisq(1,8-1))
}
#########################################################################









library(R2OpenBUGS)
data(schools)
y = schools[,2]
sigma.y = schools[,3]


##run one sim to generate plots for each setting
##run for n=10000 and n=100000
n = 100000
conf = .9
p.1 = 8

sims = matrix(nrow=n,ncol=p.1+2)
mu = 1
tau = 1
##run the sampler
for (t in 1:n){
	theta <- theta.update(mu,tau)
	mu <- mu.update(theta,tau)
	tau <- tau.update(theta,mu)
	sims[t,] <- c(theta, mu, tau)
}

X = sims[1:n,1:8]


###########################################################################
##our improvment to the stan plot, so mean, 80%CI and 95%CI
e.ind.m4q = rep(1,p.1)
q.ind.m4q = c(rep(.1,p.1),rep(.9,p.1),rep(.025,p.1),rep(.975,p.1))
col.q.m4q = c(rep(1:p.1,4))
mbm.G.m4q = mbm.g(X, e.ind.m4q, q.ind.m4q, col.q.m4q)
cov.G.m4q = mbm.G.m4q$Cov
est.G.m4q = mbm.G.m4q$Est
m.int.m4q = n.sim.int(Sigma=cov.G.m4q,conf=.9,center=est.G.m4q,
	epsilon=.001)$ints
in.cred.ible(X,est.G.m4q,m.int.m4q)

###########################################################################
##now just an 80% cred interval
e.ind.c = rep(0,p.1)
col.q.c = c(sapply(1:p.1,function(i){rep(i,2)}))
q.ind.c = c(rep(c(.1,.9),p.1))
mbm.G.c = mbm.g(X, e.ind.c, q.ind.c, col.q.c)
cov.G.c = mbm.G.c$Cov
est.G.c = mbm.G.c$Est
m.int.c = n.sim.int(Sigma=cov.G.c,conf=.9,center=est.G.c,epsilon=.001)$ints

par(mfrow=c(2,4))
xlabels = sapply(1:p.1, function(i){bquote(theta[.(i)])})
for(i in 1:p.1){
	plot.sim.up(X[,i],est.G.c[which(col.q.c==i)],
		m.int.c[which(col.q.c==i),],xlab=xlabels[i],ylab="Density")
}

par(mfrow=c(1,1))


















##########################################
m = 10000
W = sims[1:m,1:8]


###########################################################################
w.mbm.G.m4q = mbm.g(W, e.ind.m4q, q.ind.m4q, col.q.m4q)
w.cov.G.m4q = w.mbm.G.m4q$Cov
w.est.G.m4q = w.mbm.G.m4q$Est
w.m.int.m4q = n.sim.int(Sigma=w.cov.G.m4q,conf=.9,center=w.est.G.m4q,
	epsilon=.001)$ints
in.cred.ible(W,w.est.G.m4q,w.m.int.m4q)

###########################################################################
##now just an 80% cred interval
w.mbm.G.c = mbm.g(W, e.ind.c, q.ind.c, col.q.c)
w.cov.G.c = w.mbm.G.c$Cov
w.est.G.c = w.mbm.G.c$Est
w.m.int.c = n.sim.int(Sigma=w.cov.G.c,conf=.9,center=w.est.G.c,
	epsilon=.001)$ints










xlb = sapply(1:p.1, function(i){min(w.m.int.c[c((i-1)*(2)+1,i*2),])})-1
xub = sapply(1:p.1, function(i){max(w.m.int.c[c((i-1)*(2)+1,i*2),])})+1

############################################################################
############################################################################
EST.lwd = 1.25




xlabels = sapply(1:p.1, function(i){bquote(theta[.(i)])})
for(i in 1:p.1){
pdf(paste(paste(paste("credTheta",i,sep=""),"100k",sep="_"),"pdf",sep="."),
    width=2.2,height=1.6,pointsize=10)
	par(mfrow=c(1,1),mar=c(3,3,2,1),mgp=c(2,1,0))
	plot.sim.up(X[,i],est.G.c[which(col.q.c==i)],ran.range=c(xlb[i],xub[i]),
		m.int.c[which(col.q.c==i),],xlab=xlabels[i],ylab="density",
	ylim=c(0,.08),main=bquote(theta[.(i)]~" Posterior"),est.lwd=EST.lwd)
dev.off()
}

xlabels = sapply(1:p.1, function(i){bquote(theta[.(i)])})
for(i in 1:p.1){
pdf(paste(paste(paste("credTheta",i,sep=""),"10k",sep="_"),"pdf",sep="."),
    width=2.2,height=1.6,pointsize=10)
	par(mfrow=c(1,1),mar=c(3,3,2,1),mgp=c(2,1,0))
	plot.sim.up(W[,i],w.est.G.c[which(col.q.c==i)],ran.range=c(xlb[i],xub[i]),
		w.m.int.c[which(col.q.c==i),],xlab=xlabels[i],ylab="density",
	ylim=c(0,.08),main=bquote(theta[.(i)]~" Posterior"),est.lwd=EST.lwd)
dev.off()
}






pdf("100kimpstan.pdf",width=3.2,height=3.2,pointsize=8)
par(mar=c(2,2,3,1))
in.cred.ible(X,est.G.m4q,m.int.m4q,main="n = 100,000",font.main=1,LWD=1.25)
dev.off()

pdf("10kimpstan.pdf",width=3.2,height=3.2,pointsize=8)
par(mar=c(2,2,3,1))
in.cred.ible(W,w.est.G.m4q,w.m.int.m4q,main="n = 10,000",font.main=1,LWD=1.25)
dev.off()
















xlabels = sapply(1:p.1, function(i){bquote(theta[.(i)])})
pdf("100k_part1.pdf",pointsize=8,height=2.5,width=6.4)
par(mfrow=c(1,4),mar=c(4,4,2,0))
i=1
plot.sim.up(X[,i],est.G.c[which(col.q.c==i)],ran.range=c(xlb[i],xub[i]),
	m.int.c[which(col.q.c==i),],xlab=xlabels[i],ylab="Density",font.main=1,
	ylim=c(0,.08),main=bquote(theta[.(i)]~" Posterior"),
	axes=FALSE,frame.plot=TRUE,est.lwd=EST.lwd)
axis(side=2,las=1)
axis(side=1,las=1,xlim=c(xlb[i],xub[i]))
par(mar=c(4,0,2,0))
i=2
plot.sim.up(X[,i],est.G.c[which(col.q.c==i)],ran.range=c(xlb[i],xub[i]),
	m.int.c[which(col.q.c==i),],xlab=xlabels[i],ylab="",font.main=1,
	ylim=c(0,.08),main=bquote(theta[.(i)]~" Posterior"),
	axes=FALSE,frame.plot=TRUE,est.lwd=EST.lwd)
axis(side=1,las=1,xlim=c(xlb[i],xub[i]))
par(mar=c(4,0,2,0))
i=3
plot.sim.up(X[,i],est.G.c[which(col.q.c==i)],ran.range=c(xlb[i],xub[i]),
	m.int.c[which(col.q.c==i),],xlab=xlabels[i],ylab="",font.main=1,
	ylim=c(0,.08),main=bquote(theta[.(i)]~" Posterior"),
	axes=FALSE,frame.plot=TRUE,est.lwd=EST.lwd)
axis(side=1,las=1,xlim=c(xlb[i],xub[i]))
par(mar=c(4,0,2,1))
i=4
plot.sim.up(X[,i],est.G.c[which(col.q.c==i)],ran.range=c(xlb[i],xub[i]),
	m.int.c[which(col.q.c==i),],xlab=xlabels[i],ylab="",font.main=1,
	ylim=c(0,.08),main=bquote(theta[.(i)]~" Posterior"),
	axes=FALSE,frame.plot=TRUE,est.lwd=EST.lwd)
axis(side=1,las=1,xlim=c(xlb[i],xub[i]))
dev.off()

pdf("100k_part2.pdf",pointsize=8,height=2.5,width=6.4)
par(mfrow=c(1,4),mar=c(4,4,2,0))
i=5
plot.sim.up(X[,i],est.G.c[which(col.q.c==i)],ran.range=c(xlb[i],xub[i]),
	m.int.c[which(col.q.c==i),],xlab=xlabels[i],ylab="Density",font.main=1,
	ylim=c(0,.08),main=bquote(theta[.(i)]~" Posterior"),
	axes=FALSE,frame.plot=TRUE,est.lwd=EST.lwd)
axis(side=2,las=1)
axis(side=1,las=1,xlim=c(xlb[i],xub[i]))
par(mar=c(4,0,2,0))
i=6
plot.sim.up(X[,i],est.G.c[which(col.q.c==i)],ran.range=c(xlb[i],xub[i]),
	m.int.c[which(col.q.c==i),],xlab=xlabels[i],ylab="",font.main=1,
	ylim=c(0,.08),main=bquote(theta[.(i)]~" Posterior"),
	axes=FALSE,frame.plot=TRUE,est.lwd=EST.lwd)
axis(side=1,las=1,xlim=c(xlb[i],xub[i]))
par(mar=c(4,0,2,0))
i=7
plot.sim.up(X[,i],est.G.c[which(col.q.c==i)],ran.range=c(xlb[i],xub[i]),
	m.int.c[which(col.q.c==i),],xlab=xlabels[i],ylab="",font.main=1,
	ylim=c(0,.08),main=bquote(theta[.(i)]~" Posterior"),
	axes=FALSE,frame.plot=TRUE,est.lwd=EST.lwd)
axis(side=1,las=1,xlim=c(xlb[i],xub[i]))
par(mar=c(4,0,2,1))
i=8
plot.sim.up(X[,i],est.G.c[which(col.q.c==i)],ran.range=c(xlb[i],xub[i]),
	m.int.c[which(col.q.c==i),],xlab=xlabels[i],ylab="",font.main=1,
	ylim=c(0,.08),main=bquote(theta[.(i)]~" Posterior"),
	axes=FALSE,frame.plot=TRUE,est.lwd=EST.lwd)
axis(side=1,las=1,xlim=c(xlb[i],xub[i]))
dev.off()









xlabels = sapply(1:p.1, function(i){bquote(theta[.(i)])})
pdf("10k_part1.pdf",pointsize=8,height=2.5,width=6.4)
par(mfrow=c(1,4),mar=c(4,4,2,0))
i=1
plot.sim.up(W[,i],w.est.G.c[which(col.q.c==i)],ran.range=c(xlb[i],xub[i]),
	w.m.int.c[which(col.q.c==i),],xlab=xlabels[i],ylab="Density",font.main=1,
	ylim=c(0,.08),main=bquote(theta[.(i)]~" Posterior"),
	axes=FALSE,frame.plot=TRUE,est.lwd=EST.lwd)
axis(side=2,las=1)
axis(side=1,las=1,xlim=c(xlb[i],xub[i]))
par(mar=c(4,0,2,0))
i=2
plot.sim.up(W[,i],w.est.G.c[which(col.q.c==i)],ran.range=c(xlb[i],xub[i]),
	w.m.int.c[which(col.q.c==i),],xlab=xlabels[i],ylab="",font.main=1,
	ylim=c(0,.08),main=bquote(theta[.(i)]~" Posterior"),
	axes=FALSE,frame.plot=TRUE,est.lwd=EST.lwd)
axis(side=1,las=1,xlim=c(xlb[i],xub[i]))
par(mar=c(4,0,2,0))
i=3
plot.sim.up(W[,i],w.est.G.c[which(col.q.c==i)],ran.range=c(xlb[i],xub[i]),
	w.m.int.c[which(col.q.c==i),],xlab=xlabels[i],ylab="",font.main=1,
	ylim=c(0,.08),main=bquote(theta[.(i)]~" Posterior"),
	axes=FALSE,frame.plot=TRUE,est.lwd=EST.lwd)
axis(side=1,las=1,xlim=c(xlb[i],xub[i]))
par(mar=c(4,0,2,1))
i=4
plot.sim.up(W[,i],w.est.G.c[which(col.q.c==i)],ran.range=c(xlb[i],xub[i]),
	w.m.int.c[which(col.q.c==i),],xlab=xlabels[i],ylab="",font.main=1,
	ylim=c(0,.08),main=bquote(theta[.(i)]~" Posterior"),
	axes=FALSE,frame.plot=TRUE,est.lwd=EST.lwd)
axis(side=1,las=1,xlim=c(xlb[i],xub[i]))
dev.off()

pdf("10k_part2.pdf",pointsize=8,height=2.5,width=6.4)
par(mfrow=c(1,4),mar=c(4,4,2,0))
i=5
plot.sim.up(W[,i],w.est.G.c[which(col.q.c==i)],ran.range=c(xlb[i],xub[i]),
	w.m.int.c[which(col.q.c==i),],xlab=xlabels[i],ylab="Density",font.main=1,
	ylim=c(0,.08),main=bquote(theta[.(i)]~" Posterior"),
	axes=FALSE,frame.plot=TRUE,est.lwd=EST.lwd)
axis(side=2,las=1)
axis(side=1,las=1,xlim=c(xlb[i],xub[i]))
par(mar=c(4,0,2,0))
i=6
plot.sim.up(W[,i],w.est.G.c[which(col.q.c==i)],ran.range=c(xlb[i],xub[i]),
	w.m.int.c[which(col.q.c==i),],xlab=xlabels[i],ylab="",font.main=1,
	ylim=c(0,.08),main=bquote(theta[.(i)]~" Posterior"),
	axes=FALSE,frame.plot=TRUE,est.lwd=EST.lwd)
axis(side=1,las=1,xlim=c(xlb[i],xub[i]))
par(mar=c(4,0,2,0))
i=7
plot.sim.up(W[,i],w.est.G.c[which(col.q.c==i)],ran.range=c(xlb[i],xub[i]),
	w.m.int.c[which(col.q.c==i),],xlab=xlabels[i],ylab="",font.main=1,
	ylim=c(0,.08),main=bquote(theta[.(i)]~" Posterior"),
	axes=FALSE,frame.plot=TRUE,est.lwd=EST.lwd)
axis(side=1,las=1,xlim=c(xlb[i],xub[i]))
par(mar=c(4,0,2,1))
i=8
plot.sim.up(W[,i],w.est.G.c[which(col.q.c==i)],ran.range=c(xlb[i],xub[i]),
	w.m.int.c[which(col.q.c==i),],xlab=xlabels[i],ylab="",font.main=1,
	ylim=c(0,.08),main=bquote(theta[.(i)]~" Posterior"),
	axes=FALSE,frame.plot=TRUE,est.lwd=EST.lwd)
axis(side=1,las=1,xlim=c(xlb[i],xub[i]))
dev.off()















pdf("post_cred_100k.pdf",width=4,height=5,pointsize=10)
	par(mfrow=c(4,2),mar=c(3,3,2,1),mgp=c(2,1,0))
for(i in 1:p.1){
	plot.sim.up(X[,i],est.G.c[which(col.q.c==i)],ran.range=c(xlb[i],xub[i]),
		m.int.c[which(col.q.c==i),],xlab=xlabels[i],ylab="Density",
	ylim=c(0,.08),main=bquote(theta[.(i)]~" Posterior"),est.lwd=EST.lwd)
}
dev.off()



pdf("post_cred_10k.pdf",width=4,height=5,pointsize=10)
	par(mfrow=c(4,2),mar=c(3,3,2,1),mgp=c(2,1,0))
for(i in 1:p.1){
	plot.sim.up(W[,i],w.est.G.c[which(col.q.c==i)],ran.range=c(xlb[i],xub[i]),
		w.m.int.c[which(col.q.c==i),],xlab=xlabels[i],ylab="Density",
	ylim=c(0,.08),main=bquote(theta[.(i)]~" Posterior"),est.lwd=EST.lwd)
}
dev.off()











