#######################################
## Background functions for logistic
## regression example
#######################################

# Gradient Function for Logistic Model
gradnt_log <- function(y_this, a, beta)
{
  
  tmp <- as.numeric(( 1 + exp(y_this * sum(a * beta))))
  p_thet <- 1/ tmp
  
  return(-y_this * p_thet)
}

logistic_sgd <- function(n, burn_in, dta, init, alp, eta_cns, epochs)
{
  y <- as.vector(dta[ ,1])
  x <- as.matrix(dta[, -1])
  ### Conversion for applying GLM ###
  y1      <- (y[1:init] + 1) / 2 
  x1      <- x[1:init, ]
  #indx_rm <-	which(apply(x1, 2, var)==0)
  
  model_train <- glm(formula = y1 ~ x1 + 0, family = binomial)
  
  n.mod <- length(y) - init
  sg <- matrix(nrow = n.mod*epochs, ncol = dim(x)[2]) 
  temp   <- as.vector(model_train$coefficients) #MLE
  sg[1, ] <- temp
  sg[1, is.na(temp)]   <-  0
  sg[1, ] <- rep(0, length(temp))
  
  for(i in 2:(n.mod*epochs))
  {  
    ind <- sample(((init+1):n.mod), size = 1)
    sg[i,] <- (sg[(i - 1),] - eta_cns*i^(-alp) *
                 gradnt_log( y[ind], x[ind, ], sg[(i-1), ]) * x[ind,])
  }
  sg_ct_full  <- sg[(burn_in + 1):(n.mod*epochs), ]
  return(sg_ct_full)
}

