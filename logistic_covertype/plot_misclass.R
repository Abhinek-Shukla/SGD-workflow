plot_misclass <- function(cutoffs, train_percent, misclass, misclass.lb)
  {
  
pdf(file = paste0("train_", train_percent, ".pdf"),   
    width = 4, height = 4)
  
plot(cutoffs, misclass, type = "l", col = "blue", lwd = 1, 
     xlab = "Cutoff", ylab = "misclassification rate", 
     main = paste0(train_percent, "% Train ",
                  100-train_percent, "% Test"), 
                  ylim = c(min(c(misclass, misclass.lb)), 
                           max(c(misclass, misclass.lb))))

# Add the second line to the same plot
lines(cutoffs, misclass.lb, col = "red", lwd = 1)

# Add a legend
legend("topright", legend = c("misclass", "misclass.lb"), 
       col = c("blue", "red"), lwd = 2, cex = 0.5)
dev.off()

}