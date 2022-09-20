set.seed(100)
library(smoothmest)
source("qq_lad.R")
lad_qq(nparm = 5, cns = c(0.1), 
	eta_cns = 0.5, sam_siz = 5e4, alp = .51, 
  burn_in = 5000, cns1 = 0.1, Iter = 50, nsampl = 100)


lad_qq(nparm = 5, cns = c(0.1), 
	eta_cns = 0.5, sam_siz = 1e6, alp = .51, 
  burn_in = 5000, cns1 = 0.1, Iter = 50, nsampl = 100)



lad_qq(nparm = 20, cns = c(0.1), 
	eta_cns = 1, sam_siz = 5e4, alp = .51, 
  burn_in = 5000, cns1 = 0.1, Iter = 50, nsampl = 100)


lad_qq(nparm = 20, cns = c(0.1), 
	eta_cns = 1, sam_siz = 1e6, alp = .51, 
  burn_in = 5000, cns1 = 0.1, Iter = 50, nsampl = 100)


