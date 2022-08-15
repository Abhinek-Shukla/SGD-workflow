##Aug 7, 2020
##IID sampling mixnorm example with joint conf regions
##repeat for n=500,1000,5000,10000
##conf levels .8,.9,.95

rm(list=ls())
set.seed(573)
source("outputfun.R")

conf.level = c(.8,.9,.95)
cfs = length(conf.level)

n = 10000
reps = 2000

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

#######################################################################

op.it.1 = function(up){(integrate(f,lower=-Inf,upper=up)$value - .1)^2}
op.it.9 = function(up){(integrate(f,lower=-Inf,upper=up)$value - .9)^2}

t.xi.1 = optimize(op.it.1,lower=0,upper=10)$min
t.xi.9 = optimize(op.it.9,lower=10,upper=15)$min
Theta = c(sum(mix.p*mix.this),t.xi.1,t.xi.9)

#######################################################################




e.ind = c(1)
q.ind = c(.1,.9)
col.q = c(1,1)
p = sum(e.ind) + length(q.ind)

ests = matrix(nrow=reps,ncol=p)
SI.cover = matrix(nrow=reps,ncol=cfs)
LB.cover = matrix(nrow=reps,ncol=cfs)
UB.cover = matrix(nrow=reps,ncol=cfs)


start.time = Sys.time()
##now this part is replicated
for(r in 1:reps){
labels = sample(1:3,n,prob=mix.p,replace=TRUE)
x = rnorm(n,mix.this[labels],sds[labels])

iid = iid.g(x,e.ind,q.ind,col.q)

SI.ints = array(data=c(sapply(1:cfs,function(i){
	n.sim.int(iid$Cov,conf.level[i],iid$Est)$ints})),dim=c(p,2,cfs))

crit.under = qnorm(1-(1-conf.level)/2)
crit.over = qnorm(1-(1-conf.level)/2/p)
sd.diag = sqrt(diag(iid$Cov))
LB.ints = array(data=c(sapply(1:cfs,function(i){
	cbind(iid$Est-sd.diag*crit.under[i],iid$Est+sd.diag*crit.under[i])})),
	dim=c(p,2,cfs))
UB.ints = array(data=c(sapply(1:cfs,function(i){
	cbind(iid$Est-sd.diag*crit.over[i],iid$Est+sd.diag*crit.over[i])})),
	dim=c(p,2,cfs))

ests[r,] = iid$Est
SI.cover[r,] = sapply(1:cfs,function(i){cover(Theta,SI.ints[,,i])})
LB.cover[r,] = sapply(1:cfs,function(i){cover(Theta,LB.ints[,,i])})
UB.cover[r,] = sapply(1:cfs,function(i){cover(Theta,UB.ints[,,i])})

if(r%%100==0){print(r)}
}
end.time = Sys.time()
sim.time = end.time - start.time
##########################################################################

cover.rep = cbind(LB.cover,SI.cover,UB.cover)
colnames(cover.rep) = c("LB .8", "LB .9", "LB .95", "SI .8","SI .9","SI .95",
	"UB .8", "UB .9", "UB .95")

##now to look at coverage results
coverages = matrix(colMeans(cover.rep),nrow=cfs,ncol=3,byrow=FALSE)
rownames(coverages) = c(.8,.9,.95)
colnames(coverages) = c("LB","SI","UB")

sim.time
coverages















cover.samples = cover.rep[,c(1,2,4,5,7,8)]


##save restults of simulation
write.table(cover.samples,"IIDsamples10k.txt",col.names=TRUE,row.names=FALSE)






