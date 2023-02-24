set.seed(100)
library(mcmcse)
source("qq_linear.R")


linear_qq(nparm = 5, cns = 0.1, 
	eta_cns = 0.5, sam_siz = 5e4, alp = .51, 
  burn_in = 5000, cns1 = 0.1, Iter = 50)


linear_qq(nparm = 5, cns = 0.1, 
	eta_cns = 0.5, sam_siz = 1e6, alp = .51, 
  burn_in = 5000, cns1 = 0.1, Iter = 50)



# linear_qq(nparm = 20, cns = 0.1, 
# 	eta_cns = 1, sam_siz = 5e4, alp = .51, 
#   burn_in = 5000, cns1 = 0.1, Iter = 50)


# linear_qq(nparm = 20, cns = 0.1, 
# 	eta_cns = 1, sam_siz = 1e6, alp = .51, 
#   burn_in = 5000, cns1 = 0.1, Iter = 50)


