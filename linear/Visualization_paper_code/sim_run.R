set.seed(1234)
# source("sim_plots.R")
source("sim_functions.R")

pdf("simstudy100.pdf", height = 8, width = 4)
par(mfrow = c(2,1), mar=c(3, 4.5, 2.5, 0), oma = c(1, 0, 1, .8), mgp = c(3,1.5,0))
yrange <- c(0,5)
n <- 100
p <- 20
sparse <- 10
beta.star <- c(rep(0,sparse), rnorm(p-sparse, sd = 4))

out <- sim_study(p = p, n = n, N = 1e2, beta.star = beta.star)
med <- mc.mqq(x = out, cred = .50)
med.est <- med$Est[4:12]
med.cov <- med$Cov[4:12, 4:12]

intervals <- new.sim.int(Sigma = med.cov, n.sim = 1e2, center = med.est)

difficult(out, ests = med.est, title = "Reps = 100",
  ints = intervals$ints, xlabels = colnames(out), plain = TRUE, yrange = yrange, ylab = "Estimation error")
difficult(out, ests = med.est, title = "Reps = 100",
  ints = intervals$ints, xlabels = colnames(out), yrange = yrange, ylab = "Estimation error")
dev.off()

pdf("simstudy500.pdf", height = 8, width = 4)
par(mfrow = c(2,1), mar=c(3, 4.5, 2.5, 0), oma = c(1, 0, 1, .8), mgp = c(3,1.5,0))
out <- sim_study(p = p, n = n, N = 5e2, beta.star = beta.star)
med <- mbm.mqq(x = out, cred = .50)
med.est <- med$Est[4:12]
med.cov <- med$Cov[4:12, 4:12]

intervals <- new.sim.int(Sigma = med.cov, n.sim = 5e2, center = med.est)

difficult(out, ests = med.est, title = "Reps = 500",
  ints = intervals$ints, xlabels = colnames(out), plain = TRUE, yrange = yrange, ylab = "Estimation error")
difficult(out, ests = med.est, title = "Reps = 500",
  ints = intervals$ints, xlabels = colnames(out), yrange = yrange, ylab = "Estimation error")
dev.off()


pdf("simstudy2000.pdf", height = 8, width = 4)
par(mfrow = c(2,1), mar=c(3, 4.5, 2.5, 0), oma = c(1, 0, 1, .8), mgp = c(3,1.5,0))
out <- sim_study(p = p, n = n, N = 2e3, beta.star = beta.star)
med <- mbm.mqq(x = out, cred = .50)
med.est <- med$Est[4:12]
med.cov <- med$Cov[4:12, 4:12]

intervals <- new.sim.int(Sigma = med.cov, n.sim = 2e3, center = med.est)

difficult(out, ests = med.est, title = "Reps = 2000",
  ints = intervals$ints, xlabels = colnames(out), plain = TRUE, yrange = yrange, ylab = "Estimation error")
difficult(out, ests = med.est, title = "Reps = 2000",
  ints = intervals$ints, xlabels = colnames(out), yrange = yrange, ylab = "Estimation error")
dev.off()