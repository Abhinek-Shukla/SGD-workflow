ibs_jasa_mean <- function(sgd, alp = 0.51, cns = 1) 
{
  n <- nrow(sgd)
  nparm <- ncol(sgd)

  # JASA Online BM Estimators
  am <- ceiling( cns*(1:1e3)^(2/(1-alp)))
  am <- c(am[am < n], (n+1))

  bm <- diff(am)
  batch_means <- t(sapply(1:(length(am)-1), function(i) colMeans( matrix(sgd[am[i] : (am[i+1]-1), ], ncol = nparm ))))

  batch_means <- bm*scale(batch_means, center = colMeans(sgd), scale = FALSE)
  out <- (t(batch_means) %*% batch_means)/n

  return(out)
}

