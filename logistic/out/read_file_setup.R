
eta_typ <- 1
foo <- paste("out/logistic_", nam_matrix, "_n_",max_sam,"_dim_",nparm,"_eta_", eta_typ , ".RData",sep="")
load(foo)


#For dimension 20
cns = seq(0.2, 0.6, by =0.1)
colMeans(cover_ibs)
colMeans(cover_orc)
cov_ebs_matr <- array(dim = c(3,length(cns),5), dimnames = list(  c("beta1", "beta2", "beta3"), cns, 1:5))
for (k in 1: 5)
{ for(l in 1:length(cns) ){
  tmp1 <- (3*l-2)
  tmp2 <- 3*l
  cov_ebs_matr[,l, k  ] <-  colMeans(cover_ebs[k,, tmp1 : tmp2])
  
}
}
cov_ebs_ls_matr <- array(dim = c(3,length(cns),5), dimnames = list(  c("beta1", "beta2", "beta3"), cns, 1:5))
for (k in 1: 5)
{ for(l in 1:length(cns)){
  tmp1 <- (3*l-2)
  tmp2 <- 3*l
  cov_ebs_ls_matr[,l, k  ] <-  colMeans(cover_ebs_ls[k,, tmp1 : tmp2])
  
}
}


