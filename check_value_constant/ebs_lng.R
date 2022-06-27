ebs_lng <- function(n,cns,alp){
bet <- 0.5416 + 0.4671*alp - 0.6930/log10(n)
two_seq <- 2^(seq(10:40))
#Equal Batch Size Configuration
bn <- min(two_seq[two_seq >= cns*n^bet])
#No. of batches
print(bn)
an <- floor(n/bn)
return(an)
}