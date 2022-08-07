
library(smoothmest)
grad_lad <- function(sg,y,x, tau = 0.5){
  
  (tau - as.numeric( y - t(x) %*% sg < 0 )) %*% x
  
} 
st <- Sys.time()
max_sam = 1e5
Rep = 10
nparm = 20
cns = c(0.1, 1)
        ncores_par = 1
        eta_cns = 1
        sam_siz = 1e5
                    qlev = 0.95
                   
                    burn_in = 10000
               
                    

sam_siz <- sam_siz[sam_siz <= max_sam]
n <- sam_siz[length(sam_siz)] 


parm <- seq(1 / nparm, 1, length.out = nparm)



sg <- matrix(nrow = n + burn_in, ncol = nparm);
sg_ct <- matrix(nrow = n , ncol = nparm)
#Iterates stored
mse <- matrix(rep(0, 4*nparm ), nrow = 4, ncol = nparm)
for( k in 1 : Rep)
{x <- matrix(rnorm((n + burn_in) * nparm), nrow = (n + burn_in), ncol = nparm)

#noisy Observed Data
y <- x %*% parm + rdoublex(n + burn_in, mu = 0,lambda = 1)
#Learning Rate
eta <- numeric(n + burn_in)
sg[1,] <- rep(0, nparm)
alp = 0.51
for(i in 2 : (n + burn_in)){
  
  eta[i] <- i^( - alp)
  sg[i,] <- sg[i - 1,] + eta_cns * eta[i] * grad_lad(sg[i - 1,], y[i], x[i,])
  
}

sg_ct_full <- sg[(burn_in + 1) : (n + burn_in),]

mse[1, ] <-  mse[1, ] + (colMeans(sg_ct_full)-parm)^2


sg[1,] <- rep(0, nparm)
alp = 2/3
for(i in 2 : (n + burn_in)){
  
  eta[i] <- i^( - alp)
  sg[i,] <- sg[i - 1,] + eta_cns * eta[i] * grad_lad(sg[i - 1,], y[i], x[i,])
  
}

sg_ct_full <- sg[(burn_in + 1) : (n + burn_in),]

mse[2,  ] <- mse[2, ] +  (colMeans(sg_ct_full)-parm)^2

sg[1,] <- rep(0, nparm)
alp = 0.75
for(i in 2 : (n + burn_in)){
  
  eta[i] <- i^( - alp)
  sg[i,] <- sg[i - 1,] + eta_cns * eta[i] * grad_lad(sg[i - 1,], y[i], x[i,])
  
}

sg_ct_full <- sg[(burn_in + 1) : (n + burn_in),]

mse[3,  ] <- mse[3, ] +   (colMeans(sg_ct_full)-parm)^2



sg[1,] <- rep(0, nparm)
alp = 0.95
for(i in 2 : (n + burn_in)){
  
  eta[i] <- i^( - alp)
  sg[i,] <- sg[i - 1,] + eta_cns * eta[i] * grad_lad(sg[i - 1,], y[i], x[i,])
  
}

sg_ct_full <- sg[(burn_in + 1) : (n + burn_in),]

mse[4,  ] <- mse[4, ] +  (colMeans(sg_ct_full)-parm)^2

}

mse <- mse/Rep
print(mse)
print(difftime(Sys.time(), st, units = "secs"))

