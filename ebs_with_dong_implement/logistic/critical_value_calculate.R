nparm <- 20
an <- 28
crt_val_ebs <- qf(qlev, df1 = nparm, df2 = (an - nparm))
sam_siz <- c( 1e5, 2e5, 5e5, 8e5, 1e6, 5e6,1e7)
updt_crt <- vector()
for ( smpl in 1: length(sam_siz)){
bn <- floor(sam_siz[smpl]/an)
print(bn)
updt_crt[smpl] <- ( sam_siz[smpl] * nparm / ((an - nparm) * bn)) 


}
crt_val_ebs <- crt_val_ebs * updt_crt
print(crt_val_ebs)