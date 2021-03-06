############################################################################################
#
# Author:   Luigi Augugliaro
# e-mail:   luigi.augugliaro@unipa.it
# home:     http://dssm.unipa.it/augugliaro/
# data:     27-03-2018

load("Sim_ebic.RData")

measures <- array(0, dim = c(length(n), length(g), 2, 2, 2),
    dimnames = list(n = n, g = g, mes = c("mean", "sd"), type = c("ll", "Qval"),
    summary = c("TPR", "FDR")))

for(k in 1:length(n)){
    for(j in 1:5){
        measures[k, j, "mean", "ll", "FDR"] <- mean(na.omit(results[k, , j, "ll", "FP"] / (results[k, , j, "ll", "FP"] + results[k, , j, "ll", "TP"])))
        measures[k, j, "mean", "ll", "TPR"] <- mean(na.omit(results[k, , j, "ll", "TP"] / (results[k, , j, "ll", "TP"] + results[k, , j, "ll", "FN"])))
        measures[k, j, "mean", "Qval", "FDR"] <- mean(na.omit(results[k, , j, "Qval", "FP"] / (results[k, , j, "Qval", "FP"] + results[k, , j, "Qval", "TP"])))
        measures[k, j, "mean", "Qval", "TPR"] <- mean(na.omit(results[k, , j, "Qval", "TP"] / (results[k, , j, "Qval", "TP"] + results[k, , j, "Qval", "FN"])))
        measures[k, j, "sd", "ll", "FDR"] <- sd(na.omit(results[k, , j, "ll", "FP"] / (results[k, , j, "ll", "FP"] + results[k, , j, "ll", "TP"])))
        measures[k, j, "sd", "ll", "TPR"] <- sd(na.omit(results[k, , j, "ll", "TP"] / (results[k, , j, "ll", "TP"] + results[k, , j, "ll", "FN"])))
        measures[k, j, "sd", "Qval", "FDR"] <- sd(na.omit(results[k, , j, "Qval", "FP"] / (results[k, , j, "Qval", "FP"] + results[k, , j, "Qval", "TP"])))
        measures[k, j, "sd", "Qval", "TPR"] <- sd(na.omit(results[k, , j, "Qval", "TP"] / (results[k, , j, "Qval", "TP"] + results[k, , j, "Qval", "FN"])))
    }
}

######################################################
# First part of the Table 2 in Supplementary Materials
######################################################
ftable(round(measures[, c(1, 3, 5), , , ],2), col.vars = c(4, 5, 2))

##########################################################
# Second column of the Figure 1 in Supplementary Materials
##########################################################
library(latex2exp)

lgd <- c(TeX('$\\bar{BIC}_{0.0}(\\hat{E}^{\\rho})$'),
TeX('$\\bar{BIC}_{0.5}(\\hat{E}^{\\rho})$'),
TeX('$\\bar{BIC}_{1.0}(\\hat{E}^{\\rho})$'))
pr <- paste(c("0", "20", "40", "60", "80", "100"), "%", sep = "")

par(mfrow = c(2, 1), mai = c(1.02, 1.2, 0.82, 0.42), cex.axis = 1.3, cex.main = 1.5, cex.lab = 1.5, mgp = c(3, 0, 0))
matplot(n, measures[, c(1, 3, 5), "mean", "Qval", "TPR"], col = 1, type = "b",
ylim = c(0, 1), pch = 1:3, axes = FALSE, lwd = 2,
ylab = "True Positive Rate", xlab = "Sample size")
legend("topleft", legend = lgd, lty = 1:3, pch = 1:3, box.col = 0, cex = 1.3, lwd = 2)
axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), labels = pr, las = 1)
axis(1, at = n, labels = n)

matplot(n, measures[, c(1, 3, 5), "mean", "Qval", "FDR"], col = 1, type = "b",
ylim = c(-0.1, 1), pch = 1:3, axes = FALSE, lwd = 2,
ylab = "False Discovery Rate", xlab = "Sample size")
legend("topright", legend = lgd, lty = 1:3, pch = 1:3, box.col = 0, cex = 1.3, lwd = 2)
axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), labels = pr, las = 1)
axis(1, at = n, labels = n)

