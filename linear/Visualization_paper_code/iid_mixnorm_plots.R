

##March 26, 2019
##up march 27, 2020
##up Aug 7, 2020
##IID sampling mixnorm example to generate plots


rm(list=ls())
set.seed(573)
source("outputfun.R")

###########################################################################
########################################################################
mix.this = c(1,5,11)
vars = c(2.5,4,3)
sds = sqrt(vars)

mix.p = c(.3,.5,.2)

f = function(X){
out = mix.p[1]*dnorm(X,mix.this[1],sds[1])+
	mix.p[2]*dnorm(X,mix.this[2],sds[2])+
	mix.p[3]*dnorm(X,mix.this[3],sds[3])
return(out)
}

rw = function(X.curr){
	##variance closer to dist var gives less cor in chain
	eps = rnorm(1,0,3)
	X.star = X.curr + eps
	ratio = f(X.star)/f(X.curr)
	a = min(1, ratio)
	u = runif(1,0,1)
	out = ifelse(u<ratio,X.star,X.curr)
	return(out)
}

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



####################################################################

e.ind = c(1)
q.ind = c(.1,.9)
col.q = c(1,1)
p = sum(e.ind) + length(q.ind)

conf.level = .9

n = 50000
n1 = 10000
n2 = 1000

labels = sample(1:3,n,prob=mix.p,replace=TRUE)
x = rnorm(n,mix.this[labels],sds[labels])

x1 = x[1:n1]
x2 = x[1:n2]

iid = iid.g(x,e.ind,q.ind,col.q)
iid1 = iid.g(x1,e.ind,q.ind,col.q)
iid2 = iid.g(x2,e.ind,q.ind,col.q)

SI.ints = n.sim.int(iid$Cov,conf.level,iid$Est)$ints
ests = iid$Est

SI.ints1 = n.sim.int(iid1$Cov,conf.level,iid1$Est)$ints
ests1 = iid1$Est
SI.ints2 = n.sim.int(iid2$Cov,conf.level,iid2$Est)$ints
ests2 = iid2$Est

xlims = c(-2,14)
ylims = c(0,.11)




Est.lwd = .75

fw = 2.9
fh = 2.75

pdf("mixnorm_iid_1k.pdf",width=fw,height=fh,pointsize=10)
	par(mar=c(3,3,1.5,.5),mgp=c(1.5,.5,0),ps=10,cex.main=1,fin=c(fw,fh))
plot.sim.up(x2,ests2,SI.ints2,ran.range=xlims,est.lwd=Est.lwd,
	xlab=NULL,ylab="Density",ylim=ylims,main="IID (n = 1,000)",font.main=1)
dev.off()


pdf("mixnorm_iid_10k.pdf",width=fw,height=fh,pointsize=10)
	par(mar=c(3,3,1.5,.5),mgp=c(1.5,.5,0),ps=10,cex.main=1,fin=c(fw,fh))
plot.sim.up(x1,ests1,SI.ints1,ran.range=xlims,est.lwd=Est.lwd,
	xlab=NULL,ylab="Density",ylim=ylims,main="IID (n = 10,000)",font.main=1)
dev.off()


pdf("mixnorm_iid_50k.pdf",width=fw,height=fh,pointsize=10)
	par(mar=c(3,3,1.5,.5),mgp=c(1.5,.5,0),ps=10,cex.main=1,fin=c(fw,fh))
plot.sim.up(x,ests,SI.ints,ran.range=xlims,est.lwd=Est.lwd,
	xlab=NULL,ylab="Density",ylim=ylims,main="IID (n = 50,000)",font.main=1)
dev.off()

