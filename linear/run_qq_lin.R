set.seed(1220)
library(ggplot2)
source("qq_linear.R")
linear_qq(nparm = 5, cns = c(0.1), eta_cns = 0.5, sam_siz = 1e5, alp = .51, 
  burn_in = 5000, cns1 = 0.1, Iter = 50, nsampl = 1000)

