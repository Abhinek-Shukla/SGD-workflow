set.seed(1220)
source("hist_linear.R")
linear_hist(nparm = 20, cns = c(0.1), eta_cns = 0.5, sam_siz = 5e6, qlev = 0.95, alp = .51, 
  burn_in = 10000, nam_matrix = "indep", cns1 = 0.1)

