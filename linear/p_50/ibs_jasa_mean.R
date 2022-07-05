ibs_jasa_mean <- function(sg_ct,alp) {
n <- nrow(sg_ct)
nparm <- ncol(sg_ct)
# JASA Online BM Estimators
for(m in 1:1000){
  
  am[m] <- floor(0.005*m^(2/(1-alp)))
  
}
am <- c(am[am < n],n)
print((am))
ibs_jasa <- matrix(rep(0,nparm*(length(am)-1)),nrow=(length(am)-1),ncol=nparm)
tot_mean <- rep(0,nparm)

#Equal batch size smart batching (EBS)
for(k in 1 : (length(am)-1)){
  strt_pt <- am[k]
  end_pt <- am[k+1]-1
  if(k==(length(am)-1)){ 
    end_pt <- am[k+1]
  }
  
  ibs_jasa[k,] <-  colSums(sg_ct[strt_pt : end_pt,])
  tot_mean <- tot_mean + ibs_jasa[k,]
}
tot_mean <- tot_mean/(n)
ibs_jasa_mean <- matrix(rep(0,nparm^2),nrow=nparm,ncol=nparm)

for(k in 1 : (length(am)-1)){
  
  tmp_vl <- (ibs_jasa[k,]-(am[k+1]-am[k])*tot_mean)
 
   if(k == (length(am) - 1)){   
  tmp_vl <- (ibs_jasa[k,]-(am[k+1]-am[k]+1)*tot_mean)
    }
 
   ibs_jasa_mean <- ibs_jasa_mean + tmp_vl%*%t(tmp_vl)
  
}
ibs_jasa_mean <- ibs_jasa_mean/(n)
return(ibs_jasa_mean)
}