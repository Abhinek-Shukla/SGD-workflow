#######################################
## Background functions for logistic
## regression example
#######################################

# Gradient Function for Logistic Model
gradnt_log <- function(y, a, sg)
{
  
  tmp <- as.numeric(( 1 + exp(  y * sum(a * sg))))
  p_thet <- 1/ tmp
  
  return(-y * p_thet)
}

logistic_sgd <- function(n, burn_in, dta, init, alp, eta_cns)
{
	y <- as.vector(dta[1:(n + burn_in + init), 1])
	x <- as.matrix(dta[1:(n + burn_in + init), -1])
	
	sg <- matrix(nrow = n + burn_in, ncol = dim(x)[2]) 
	### Conversion for applying GLM ###
	y1 <- (y[1:init] + 1) / 2 
	x1 <- x[1:init, ]
	model_train <- glm(formula = y1 ~ x1 + 0, family = binomial)
	sg[1,]      <- as.vector(model_train$coefficients) #MLE

	# SGD iterates 
	eta <- numeric(n + burn_in)
	for(i in 2:(n + burn_in))
	{  
	  eta[i] <- eta_cns*i^(-alp)
	  sg[i,] <- sg[(i - 1),] - eta[i] *
	    gradnt_log( y[(init + i)], x[(init + i),], sg[(i-1), ]) * x[(init + i),]
	}
	sg_ct_full  <- sg[(burn_in + 1):(n + burn_in), ]
	return(sg_ct_full)
}

