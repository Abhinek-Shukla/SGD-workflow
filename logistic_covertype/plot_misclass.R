plot_misclass <- function(cutoffs, misclass, misclass.lb)
  {
  
pdf(file = "misclass_covtype.pdf",   
    width = 4, height = 4)
  
plot(cutoffs, misclass, type = "l", col = "black", lwd = 1, 
     xlab = "Softmax Threshold", ylab = "Misclassification Rate", 
                  ylim = c(min(c(misclass, misclass.lb)), 
                           max(c(misclass, misclass.lb))))

# Add the second line to the same plot
lines(cutoffs, misclass.lb, col = "blue", lwd = 1)

# Add a legend
legend("topright", legend = c("No confidence intervals", "With confidence intervals"), 
       col = c("black", "blue"), lwd = 1, cex = 0.5)
dev.off()

}