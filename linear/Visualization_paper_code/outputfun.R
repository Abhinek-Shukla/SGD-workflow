
##some new functions used in mcmc output analysis

library(mvtnorm)
library(mcmcse)


n.sim.int = function(Sigma, conf = .9, center=rep(0,dim(Sigma)[1]),
 epsilon = .001, n.star = 100){

	p = dim(Sigma)[1]
	crit.under = qnorm(1-(1-conf)/2)
	crit.over = qnorm(1-(1-conf)/2/p)

	i = 1
	LB.i = crit.under*sqrt(diag(Sigma))
	UB.i = crit.over*sqrt(diag(Sigma))
	X.i = (LB.i+UB.i)/2
	F.X.i = pmvnorm(lower=c(-X.i),upper=c(X.i),mean=rep(0,p),sigma=Sigma)

	while(abs(F.X.i[1] - conf) > epsilon & i <= n.star){
		i = i + 1
		##equality should never happen since then it would terminate?
		if(F.X.i[1] - conf < 0){
			LB.i = X.i
			UB.i = UB.i
			X.i = (X.i + UB.i)/2
		}
		if(F.X.i[1] - conf > 0){
			LB.i = LB.i
			UB.i = X.i
			X.i = (X.i + LB.i)/2
		}
		F.X.i = pmvnorm(lower=c(-X.i),upper=c(X.i),
			mean=rep(0,p),sigma=Sigma)
	}

	out = list("prob_level_error" = F.X.i[1] - conf, "iterations" = i,
		"pmvnorm_error" = attributes(F.X.i)$error,
		"ints" = cbind(center - X.i, center + X.i))
	return(out)
}



mbm.g = function(x, e.ind, q.ind, col.q){

	if(is.matrix(x)==FALSE){x<-matrix(x,nrow=length(x),ncol=1)}
	if(length(e.ind)!=dim(x)[2]){print("Error, check inputs")}
	p.1 = sum(e.ind)
	p.2 = length(q.ind)
	p = p.1 + p.2
	n = dim(x)[1]

	##Find our estimates, call vector Theta.n
	xi.hat = sapply(1:p.2,function(i){quantile(x[,col.q[i]],q.ind[i])})
	

	if(sum(e.ind)==0){mu.n=NULL}else{
		mu.n = apply(x,2,mean)[which(e.ind==1)]
	}

	Theta.n = c(mu.n,xi.hat)


	##generate our columns relating to quantiles
	x.I = sapply(1:p.2,function(i){ifelse(x[,col.q[i]]>xi.hat[i],1,0)})

	X.g = cbind(x[,which(e.ind==1)],x.I)


	sigma.n = mcse.multi(X.g)$cov
	a = floor(n/floor(sqrt(n)))
	df = a
	ESS = multiESS(X.g,covmat=sigma.n)


	den.xi = sapply(1:p.2,function(i){
		density(x[,col.q[i]],from=xi.hat[i],to=xi.hat[i],n=1)$y})

	Delt = diag(c(rep(1,p.1),1/den.xi))
	Omega.n = t(Delt)%*%sigma.n%*%(Delt)
	Cov.n = Omega.n/n
	out = list("Cov" = Cov.n, "Est" = Theta.n, "df" = df, "n" = n,
		"pre.delt.cov" = sigma.n, "ESS" = ESS, "p" = p, "p.1" = p.1,
		"p.2" = p.2, "which.expectations" = e.ind, 
		"quantiles.of.col" = col.q, "quantile.prob" = q.ind)
	return(out)
}



iid.g = function(x, e.ind, q.ind, col.q){

	if(is.matrix(x)==FALSE){x<-matrix(x,nrow=length(x),ncol=1)}
	if(length(e.ind)!=dim(x)[2]){print("Error, check inputs")}
	p.1 = sum(e.ind)
	p.2 = length(q.ind)
	p = p.1 + p.2
	n = dim(x)[1]

	##Find our estimates, call vector Theta.n
	xi.hat = sapply(1:p.2,function(i){quantile(x[,col.q[i]],q.ind[i])})
	

	if(sum(e.ind)==0){mu.n=NULL}else{
		mu.n = apply(x,2,mean)[which(e.ind==1)]
	}

	Theta.n = c(mu.n,xi.hat)


	##generate our columns relating to quantiles
	x.I = sapply(1:p.2,function(i){ifelse(x[,col.q[i]]>xi.hat[i],1,0)})

	X.g = cbind(x[,which(e.ind==1)],x.I)

	sigma.n = cov(X.g)

	den.xi = sapply(1:p.2,function(i){
		density(x[,col.q[i]],from=xi.hat[i],to=xi.hat[i],n=1)$y})

	Delt = diag(c(rep(1,p.1),1/den.xi))
	Omega.n = t(Delt)%*%sigma.n%*%(Delt)
	Cov.n = Omega.n/n
	out = list("Cov" = Cov.n, "Est" = Theta.n, "n" = n,
		"pre.delt.cov" = sigma.n, "p" = p, "p.1" = p.1,
		"p.2" = p.2, "which.expectations" = e.ind, 
		"quantiles.of.col" = col.q, "quantile.prob" = q.ind)
	return(out)
}


#############################################################################
##function to cal coverage of an individual 
cover = function(Theta,ints){
	p = length(Theta)
	is.cover = numeric(p)
	for(i in 1:p){
		is.cover[i] = ifelse(ints[i,1]<Theta[i] & Theta[i]<ints[i,2],1,0)
	}
	return(ifelse(sum(is.cover) == p,1,0))
}
#############################################################################