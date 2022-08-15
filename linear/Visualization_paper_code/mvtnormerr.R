
##11/25/2018
##testing the error in the pmvnorm function in package mvtnorm
library(mvtnorm)

##testing error rates of pmvnorm



test.mvn = function(p){
	sig.half = matrix(rnorm(p^2,0,2),nrow=p)
	sig = sig.half%*%t(sig.half)

	hw = runif(p,1,p^(1.75))
	low = c(-hw)
	upp = c(hw)
	prob.mvn = pmvnorm(low,upp,mean=rep(0,p),sigma=sig)
	mv.err = attributes(prob.mvn)$error
	mv.msg = attributes(prob.mvn)$msg
	mv.prob = prob.mvn[1]
	
	out = list("err" = mv.err, "msg" = mv.msg, "prob" = mv.prob)
	return(out)
}


set.seed(573)
start.time = Sys.time()
mvn.err.2 = t(sapply(rep(2,1000),test.mvn))
mvn.err.3 = t(sapply(rep(3,1000),test.mvn))
mvn.err.4 = t(sapply(rep(4,1000),test.mvn))
mvn.err.5 = t(sapply(rep(5,1000),test.mvn))
mvn.err.10 = t(sapply(rep(10,1000),test.mvn))
mvn.err.20 = t(sapply(rep(20,1000),test.mvn))
mvn.err.50 = t(sapply(rep(50,1000),test.mvn))
mvn.err.100 = t(sapply(rep(100,1000),test.mvn))
end.time = Sys.time()
sim.time = end.time - start.time

errs = matrix(unlist(c(mvn.err.2[,1],mvn.err.3[,1],mvn.err.4[,1],
	mvn.err.5[,1],mvn.err.10[,1],mvn.err.20[,1],mvn.err.50[,1],
	mvn.err.100[,1])),nrow=1000,byrow=FALSE)
msgs = matrix(unlist(c(mvn.err.2[,2],mvn.err.3[,2],mvn.err.4[,2],
	mvn.err.5[,2],mvn.err.10[,2],mvn.err.20[,2],mvn.err.50[,2],
	mvn.err.100[,2])),nrow=1000,byrow=FALSE)
probs = matrix(unlist(c(mvn.err.2[,3],mvn.err.3[,3],mvn.err.4[,3],
	mvn.err.5[,3],mvn.err.10[,3],mvn.err.20[,3],mvn.err.50[,3],
	mvn.err.100[,3])),nrow=1000,byrow=FALSE)

apply(errs,2,mean)
apply(errs,2,var)

##see which entries have problems
msgs[which(msgs!="Normal Completion")]

max(errs)
min(errs)
which(errs==max(errs))

which(msgs!="Normal Completion")
summary(errs)
summary(probs)

par(mfrow=c(4,4))
sapply(1:8,function(i){hist(errs[,i],freq=FALSE)})
sapply(1:8,function(i){hist(probs[,i],freq=FALSE)})


dims = c(2,3,4,5,10,20,50,100)
par(mfrow=c(2,4))
invisible(sapply(1:8,function(i){
	plot(probs[,i],errs[,i], main=bquote("p" == .(dims[i])),
		xlab="probability", ylab = "error")
}))



