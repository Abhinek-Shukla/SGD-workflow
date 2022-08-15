
##code will generate output files in working directory

##"outputfun.R" contains some functions for our output analysis that will be used throughout 
## will be called in each script as needed

##Generate Ellipse/upper/lower bounds plot
source("figsimconf.R")

##Generate Mixture normal plots used in intro and example
source("iid_mixnorm_plots.R")
source("mcmc_mixnorm_plots.R")

##Generate Coverage Probabilities example run
##Involves replicated simulation and may take some time
source("IID500.R")
source("IID1k.R")
source("IID5k.R")
source("IID10k.R")
source("MCMC2500.R")
source("MCMC5k.R")
source("MCMC10k.R")
source("MCMC25k.R")
source("MCMC50k.R")
source("simcoverplots.R")

##Generate boxplots for Regression Example
source("sim_functions.R")
source("sim_run.R")

##Generate school dataset plots
source("school_plots.R")

##Testing error rates of pmvnorm
##see "mvtnormerr.R"