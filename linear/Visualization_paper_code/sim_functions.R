###############################################################
## This file has some of the functions that I need repeatedly
## in my simulations and analysis. These are NOT all of the 
## functions
###############################################################


####################################################################
## Simulation study to compare lasso and ols estimates 
##
## p = number of covariates
## n = number of observations
## N = Monte Carlo replications
## theta = ar(1) correlation for the covariates
## beta.star = truth
####################################################################
#library(glmnet)
#library(mvtnorm)
library(mcmcse)
sim_study <- function(p = 25, n = 100, N = 1e2, theta = .90, beta.star = NULL)
{

  if(is.null(beta.star)){ beta.star <- rnorm(p, sd = p^(-1/2)) }
  omega <- matrix( ,ncol = p-1, nrow = p-1)
  for(i in 1:(p-1))
  {
  for(j in 1:(p-1))
  {
    omega[i,j] <- theta^(abs(i-j))
  }
  }

  ## Eigenvalue decomposition to find Sigma^(1/2)
  eig <- eigen(omega, symmetric=TRUE)
  omega.root <- eig$vectors%*%diag(eig$values^(1/2))%*%t(eig$vectors)

  lam.vec <- 10^seq(-8,8,by = 0.5 )
  mse <- matrix(0, ncol = 3, nrow = N)
  for(i in 1:N)
  {
    ## Generating data using the fact that X = 0 + Z\Sigma^(1/2)
  X <- cbind(1,matrix(rnorm(n*(p-1)), nrow = n, ncol = (p-1))%*%omega.root)
  y <- X%*%beta.star + rnorm(n, mean=0, sd=sqrt(1))

  lam.lasso <- cv.glmnet(x = X, y = y, lambda = lam.vec)$lambda.min
  lasso <- glmnet(x = X, y = y, lambda = lam.lasso)$beta

  lam.ridge <- cv.glmnet(x = X, y = y, lambda = lam.vec, alpha = 0)$lambda.min
  ridge <- glmnet(x = X, y = y, lambda = lam.ridge, alpha = 0)$beta

  ols <- lm(y ~ X -1)$coef

  mse[i, ] <- c(sum((lasso - beta.star)^2), sum((ridge - beta.star)^2), sum((ols - beta.star)^2))

  }
  colnames(mse) <- c("LASSO", "RIDGE", "OLS")
  return(mse)
}

###############################################################
## Function that estimates the covariance matrix
## x = Monte Carlo samples
## cred = credible intervals needed
## n.star = maximum number of iterations for the algorithm
###############################################################

mc.mqq <- function(x, cred=.80)
{
  if(is.matrix(x)==FALSE){x<-matrix(x,nrow=length(x),ncol=1)}
  p.1 <- dim(x)[2]
  n <- dim(x)[1]
  q.1 <- (1-cred)/2
  q.2 <- 1-q.1
  q.3 <- .5
  p <- 4*p.1

  xi.1.n <- apply(x,2,quantile,q.1)
  xi.2.n <- apply(x,2,quantile,q.2)
  xi.3.n <- apply(x,2,quantile,q.3)
  mu.n <- apply(x,2,mean)
  Theta.n <- c(mu.n,xi.3.n,xi.1.n,xi.2.n)

  # track <- rep(1:p.1, 4)
  # q.track <- rep(c(0,.5, q.1, q.2), each = 2)

  # omega.n <- matrix(0, nrow = p, ncol = p)
  # for(i in (p.1+1):p)
  # {
  #   c1 <- track[i]
  #   for(j in (p.1+1):p)
  #   {
  #     c2 <- track[j]
  #     if(c1 == c2)
  #     {
  #       omega.n[i,j] <- min(q.track[c(i,j)])*(1 - max(q.track[c(i,j)]))
  #     } else{
  #       omega.n[i,j] <- sum(colSums(apply(x[, c(c1,c2)], 1, function(t) t <= Theta.n[c(i, j)])) ==2)/n- prod(q.track[c(i,j)])
  #       print(omega.n[i,j])
  #     }      
  #   }
  # }

# sum(colSums(apply(x[, c(c1,c2)], 1, function(t) t <= Theta.n[c(i, j)])) ==2)/n
  x.q1 <- NULL
  for(j in 1:p.1)
  {
    x.j.q1 <- ifelse(x[,j] <= Theta.n[j+p.1], 1, 0)
    x.q1 <- cbind(x.q1,x.j.q1)
  }
  x.q2 <- NULL
  for(j in 1:p.1){
    x.j.q2 <- ifelse(x[,j] <= Theta.n[j+2*p.1], 1, 0)
    x.q2 <- cbind(x.q2,x.j.q2)
  }
  x.q3 <- NULL
  for(j in 1:p.1){
    x.j.q3 <- ifelse(x[,j] <= Theta.n[j+3*p.1], 1, 0)
    x.q3 <- cbind(x.q3,x.j.q3)
  }

  X.mqq <- cbind(x, x.q3, x.q1, x.q2)
  omega.n <- cov(X.mqq)

  den.q1 <- numeric(p.1)
  den.q2 <- numeric(p.1)
  den.q3 <- numeric(p.1)
  for(j in 1:p.1)
  {
    den.q1[j] <- density(x[,j], from = xi.1.n[j], to = xi.1.n[j], n=1)$y
    den.q2[j] <- density(x[,j], from = xi.2.n[j], to = xi.2.n[j], n=1)$y
    den.q3[j] <- density(x[,j], from = xi.3.n[j], to = xi.3.n[j], n=1)$y
  }
  Delt <- diag(c(rep(1,p.1),-1/den.q3,-1/den.q1,-1/den.q2))
  Sigma.n <- t(Delt)%*%omega.n%*%(Delt)
  Cov.n <- Sigma.n/n
  out = list("Cov" = Cov.n, "Est" = Theta.n, "df" = df, "n" = n,
    "pre.delt.cov" = omega.n, "p" = p, "p.1" = p.1)
  return(out)
}

###############################################################
## Function that estimates the covariance matrix
## x = Monte Carlo samples
## cred = credible intervals needed
## n.star = maximum number of iterations for the algorithm
###############################################################

mbm.mqq <- function(x, cred=.80)
{
  if(is.matrix(x)==FALSE){x<-matrix(x,nrow=length(x),ncol=1)}
  p.1 <- dim(x)[2]
  n <- dim(x)[1]
  q.1 <- (1-cred)/2
  q.2 <- 1-q.1
  q.3 <- .5
  p <- 4*p.1

  xi.1.n <- apply(x,2,quantile,q.1)
  xi.2.n <- apply(x,2,quantile,q.2)
  xi.3.n <- apply(x,2,quantile,q.3)
  mu.n <- apply(x,2,mean)
  Theta.n <- c(mu.n,xi.3.n,xi.1.n,xi.2.n)

  x.q1 <- NULL
  for(j in 1:p.1)
  {
    x.j.q1 <- ifelse(x[,j] <= Theta.n[j+p.1], 1, 0)
    x.q1 <- cbind(x.q1,x.j.q1)
  }
  x.q2 <- NULL
  for(j in 1:p.1){
    x.j.q2 <- ifelse(x[,j] <= Theta.n[j+2*p.1], 1, 0)
    x.q2 <- cbind(x.q2,x.j.q2)
  }
  x.q3 <- NULL
  for(j in 1:p.1){
    x.j.q3 <- ifelse(x[,j] <= Theta.n[j+3*p.1], 1, 0)
    x.q3 <- cbind(x.q3,x.j.q3)
  }
  X.mqq <- cbind(x, x.q3, x.q1, x.q2)

  ## IID sampling so BM is the only one that makes sense (since low variance)
  mcse.obj <- mcse.multi(X.mqq, size = "cuberoot")
  omega.n <- mcse.obj$cov
  a <-  floor(n/mcse.obj$size)
  df <- a
  ESS <- multiESS(X.mqq, covmat=omega.n)

  den.q1 <- numeric(p.1)
  den.q2 <- numeric(p.1)
  den.q3 <- numeric(p.1)
  for(j in 1:p.1)
  {
    den.q1[j] <- density(x[,j], from = xi.1.n[j], to = xi.1.n[j], n=1)$y
    den.q2[j] <- density(x[,j], from = xi.2.n[j], to = xi.2.n[j], n=1)$y
    den.q3[j] <- density(x[,j], from = xi.3.n[j], to = xi.3.n[j], n=1)$y
  }
  Delt <- diag(c(rep(1,p.1),-1/den.q3,-1/den.q1,-1/den.q2))
  Sigma.n <- t(Delt)%*%omega.n%*%(Delt)
  Cov.n <- Sigma.n/n
  out = list("Cov" = Cov.n, "Est" = Theta.n, "df" = df, "n" = n,
    "pre.delt.cov" = omega.n, "ESS" = ESS, "p" = p, "p.1" = p.1)
  return(out)
}



###############################################################
## Function produces simultaenous intervals
## n.sim  = to calculate the degrees of freedom
## Sigma = estimated asymptotic covariance matrix
## conf = confidence level for simultaneous regions
## center = the mean estimate around which to create intervals
## epsilon = tolerance level for SI algorithm
## n.star = maximum number of iterations for the algorithm
###############################################################


new.sim.int <-  function(Sigma, conf = .90, center = rep(0,dim(Sigma)[1]),
 epsilon = .001,  n.star = 100)
{

  p <- dim(Sigma)[1]
  crit.under <- qnorm(1-(1-conf)/2)
  crit.over <- qnorm(1-(1-conf)/2/p)

  i <- 1
  LB.i <- crit.under*sqrt(diag(Sigma))
  UB.i <- crit.over*sqrt(diag(Sigma))
  X.i <- (LB.i+UB.i)/2
  F.X.i <- pmvnorm(lower=c(-X.i),upper=c(X.i),sigma=Sigma)

  while(abs(F.X.i[1] - conf) > epsilon & i <= n.star)
  {
    i = i + 1
    ##equality should never happen since then it would terminate?
    if(F.X.i[1] - conf < 0)
    {
      LB.i <- X.i
      UB.i <- UB.i
      X.i <- (X.i + UB.i)/2
    }
    if(F.X.i[1] - conf > 0)
    {
      LB.i <- LB.i
      UB.i <- X.i
      X.i <- (X.i + LB.i)/2
    }
    F.X.i <- pmvnorm(lower = c(-X.i), upper = c(X.i), sigma = Sigma)
  }

  out <- list("prob_level_error" = F.X.i[1] - conf, "iterations" = i,
    "pmvt_error" = attributes(F.X.i)$error,
    "ints" = cbind(center - X.i, center + X.i))
  return(out)
}


difficult <- function(X, ests, ints, plain = FALSE, title=NULL, xlabels = NULL, range = 1.5, yrange = NULL, ylab = NULL)
{

  p <- length(ests)
  p.1 <- dim(X)[2]

  low.50 <- ests[(p.1+1):(2*p.1)]
  upp.50 <- ests[(2*p.1+1):(3*p.1)]
  # low.95 <- ests[(3*p.1+1):(4*p.1)]
  # upp.95 <- ests[(4*p.1+1):(5*p.1)]


  whisk.length <- range*(upp.50 - low.50)
  if(is.null(yrange))
  {
    yrange <- c(max(0, min(.8*(low.50 - whisk.length))), 1.1*max(upp.50 + whisk.length) )
  }


  mins <- apply(X, 2, min)
  maxs <- apply(X, 2, max)

  whisk.up <- pmin(maxs, upp.50 + whisk.length) 
  whisk.dn <- pmax(mins, low.50 - whisk.length)

  index <- rep(1:p.1,3)
  qw = .2


  par(bg="white")
  plot(NULL, ylim=yrange,xlim=c(.5,p.1+.5),xaxt="n", ann=FALSE)
  title(main = eval(title), cex.lab = 1.2, ylab = ylab, font.main = 1)
  # abline(v=seq(round(xrange[1]),round(xrange[2]),10),col="lightsteelblue")
  # abline(v=1:p.1,col="lightsteelblue")
  axis(1, at=1:p.1, labels= xlabels, las=1)

  if(plain == FALSE)
  {
    cols = adjustcolor(c(rep("lightskyblue1",p.1),rep("lightskyblue1",2*p.1),rep("lightskyblue1",2*p.1)), alpha.f = .6)
    bord.col = adjustcolor("black", alpha.f = .3)
    for(i in 1:p){
      ys = c(ints[i,1],ints[i,1],ints[i,2],ints[i,2])
      xs=c(index[i]-qw,index[i]+qw,index[i]+qw,index[i]-qw)
      polygon(xs,ys,col=cols[i], border = NA)
    }
    segments(index-qw,ests,index+qw,ests,col="black", lwd=1.5)
  }

  # segments(mins,index,low.50,index, lwd=2)
  # segments(maxs,index,upp.50,index, lwd=2)

  segments(index, whisk.dn,index,low.50, lwd=1, lty = 2)
  segments(index, whisk.up,index,upp.50, lwd=1, lty = 2)
  segments(index[-(1:p.1)] - .1, whisk.up, index[-(1:p.1)] + .1, whisk.up)
  segments(index[-(1:p.1)] - .1, whisk.dn, index[-(1:p.1)] + .1, whisk.dn)

  for(i in 1:p.1)
  {
    foo <- X[ ,i] > whisk.up[i]
    points(rep(i, sum(foo)), X[foo, i], cex = .80)
    obx = c(low.50[i],low.50[i],upp.50[i],upp.50[i])
    oby = c(index[i]-qw,index[i]+qw,index[i]+qw,index[i]-qw)
    polygon(oby, obx)
  }

segments(index-qw,ests,index+qw,ests,col="black", lwd=1.5)

}






