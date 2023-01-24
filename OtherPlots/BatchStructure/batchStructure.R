n <- seq(1e2, 1e6, by = 200)
alpha <- .51
beta <- (1 + alpha)/2
c <- .1


bnfun <- function(this.n, beta)
{
  upp <- c*this.n^beta
  bn.vector <- 2^seq(1, 100, by = 1)
  gamma <- which(bn.vector >= upp)[1]
  return(c(2^gamma, this.n/(2^gamma)))
}

out <- sapply(n, bnfun, beta = beta)
bn <- out[1, ]
an <- out[2, ]

pdf("batchsize.pdf", height = 5, width = 5)
plot(n, 2*c*n^(beta), type = 'l', 
     ylab = "Batch Size", xlab = "Iteration Length",
     lty = 3, col = "purple", lwd = 1.1)
lines(n, c*n^(beta), lty = 3, col = "purple", lwd = 1.1)
lines(n, bn, col = "black")
dev.off()

pdf("nbatch.pdf", height = 5, width = 5)
plot(n, n^(1-beta)/c, type = 'l',
     ylab = "Number of Batches", xlab = "Iteration Length",
     lty = 3, col = "purple", lwd = 1.1)
lines(n, n^(1 - beta)/c/2, lty = 3, col = "purple", lwd = 1.1)
lines(n, an, col = "black")
dev.off()