load(model_1_H25.RData)

#####################################################################
# Output
#####################################################################

prcurve_a <- array(0, dim = c(nsim, nrho, 2, 3), dimnames = list(nsim = NULL, nrho = NULL, type = c("recall", "precision"), estimator = c("rcglasso", "glasso", "mglasso")))
prcurve_data <- array(0, dim = c(30, 2, 3), dimnames = list(rho = 1:30, type = c("recall", "precision"), estimator = c("rcglasso", "glasso", "mglasso")))
auc_tbl <- array(0, dim = c(nsim, 3), dimnames = list(nsim = NULL, estimator = c("rcglasso", "glasso", "mglasso")))
mse_a <- array(0, dim = c(nsim, nrho, 2, 3), dimnames = list(nsim = NULL, nrho = NULL, par = c("muh", "thetah"), estimator = c("rcglasso", "glasso", "mglasso")))
mse_tbl <- array(0, dim = c(2, 2, 3), dimnames = list(measure = c("mu", "theta"), type = c("mean", "sd"), estimator = c("rcglasso", "glasso", "mglasso")))

#####################################################################
# Starting Analysis
#####################################################################

U1 <- upper.tri(ThetaTRUE, diag = TRUE)
for(j in 1:nsim){
    #################################
    # MSE section
    mse_a[j, , "muh", "rcglasso"] <- apply(muh_rcglasso[j, , ], 2, function(muh) crossprod(muh - mu))
    mse_a[j, , "muh", "mglasso"] <- apply(muh_mglasso[j, , ], 2, function(muh) crossprod(muh - mu))
    mse_a[j, , "thetah", "rcglasso"] <- apply(thetah_rcglasso[j, , , ], 3, function(M) sum((M[U1] - ThetaTRUE[U1])^2))
    mse_a[j, , "thetah", "glasso"] <- apply(thetah_glasso[j, , , ], 3, function(M) sum((M[U1] - ThetaTRUE[U1])^2))
    mse_a[j, , "thetah", "mglasso"] <- apply(thetah_mglasso[j, , , ], 3, function(M) sum((M[U1] - ThetaTRUE[U1])^2))
    
    #################################
    # Precizion-Recall Curve
    
    # rcglasso
    out_prcurve <- prcurve(thetah = thetah_rcglasso[j, , , ], thetat = ThetaTRUE, rho = rho_rcglasso[j, ])
    prcurve_a[j, , "recall", "rcglasso"] <- out_prcurve$recall
    prcurve_a[j, , "precision", "rcglasso"] <- out_prcurve$precision
    auc_tbl[j, "rcglasso"] <- out_prcurve$auc
    
    # glasso
    out_prcurve <- prcurve(thetah = thetah_glasso[j, , , ], thetat = ThetaTRUE, rho = rho_rcglasso[j, ])
    prcurve_a[j, , "recall", "glasso"] <- out_prcurve$recall
    prcurve_a[j, , "precision", "glasso"] <- out_prcurve$precision
    auc_tbl[j, "glasso"] <- out_prcurve$auc
    
    # missglasso
    out_prcurve <- prcurve(thetah = thetah_mglasso[j, , , ], thetat = ThetaTRUE, rho = rho_rcglasso[j, ])
    prcurve_a[j, , "recall", "mglasso"] <- out_prcurve$recall
    prcurve_a[j, , "precision", "mglasso"] <- out_prcurve$precision
    auc_tbl[j, "mglasso"] <- out_prcurve$auc
}

prcurve_data[, "recall", "glasso"] <- apply(prcurve_a[, , "recall", "glasso"], 2, mean)
prcurve_data[, "precision", "glasso"] <- apply(prcurve_a[, , "precision", "glasso"], 2, mean)
prcurve_data[, "recall", "rcglasso"] <- apply(prcurve_a[, , "recall", "rcglasso"], 2, mean)
prcurve_data[, "precision", "rcglasso"] <- apply(prcurve_a[, , "precision", "rcglasso"], 2, mean)
prcurve_data[, "recall", "mglasso"] <- apply(prcurve_a[, , "recall", "mglasso"], 2, mean)
prcurve_data[, "precision", "mglasso"] <- apply(prcurve_a[, , "precision", "mglasso"], 2, mean)

#####################################################################
# Table + Figures
#####################################################################

#########################
# Table 2 (Main Document)
mse_tbl["mu", "mean", "rcglasso"] <- mean(apply(mse_a[, , "muh", "rcglasso"], 1, min))
mse_tbl["mu", "sd", "rcglasso"] <- sd(apply(mse_a[, , "muh", "rcglasso"], 1, min))
mse_tbl["mu", "mean", "mglasso"] <- mean(apply(mse_a[, , "muh", "mglasso"], 1, min))
mse_tbl["mu", "sd", "mglasso"] <- sd(apply(mse_a[, , "muh", "mglasso"], 1, min))

mse_tbl["theta", "mean", "rcglasso"] <- mean(apply(mse_a[, , "thetah", "rcglasso"], 1, min))
mse_tbl["theta", "sd", "rcglasso"] <- sd(apply(mse_a[, , "thetah", "rcglasso"], 1, min))
mse_tbl["theta", "mean", "glasso"] <- mean(apply(mse_a[, , "thetah", "glasso"], 1, min))
mse_tbl["theta", "sd", "glasso"] <- sd(apply(mse_a[, , "thetah", "glasso"], 1, min))
mse_tbl["theta", "mean", "mglasso"] <- mean(apply(mse_a[, , "thetah", "mglasso"], 1, min))
mse_tbl["theta", "sd", "mglasso"] <- sd(apply(mse_a[, , "thetah", "mglasso"], 1, min))

tbl <- cbind(as.matrix(ftable(mse_tbl, col.vars = c("measure", "estimator"), row.vars = "type"))[, -2], rbind(apply(auc_tbl, 2, mean), apply(auc_tbl, 2, sd)))
xtable(tbl, digits = 2)

#########################
# Fig. 1 in Main Document and Fig. 2 in Supplementary Materials

# precision-recall curve
pdf(file = "prcurve-model1-H25.pdf")
par(mai = c(1.02, 1.2, 0.82, 0.42), cex.axis = 1.3, cex.main = 1.5, cex.lab = 1.5)
xlim <- c(0, 1)
par(mai = c(1.02, 1.2, 0.82, 0.42), cex.axis = 1.3, cex.main = 1.5, cex.lab = 1.5)
plot(prcurve_data[, , "rcglasso"], xlim = xlim, ylim = c(0, 1),  type = "l", lty = 1, lwd = 2, axes = FALSE, xlab = "Recall", ylab = "Precision")
points(prcurve_data[, , "glasso"], xlim = xlim, ylim = c(0, 1), type = "l", lty = 2, lwd = 2)
points(prcurve_data[, , "mglasso"], xlim = xlim, ylim = c(0, 1),  type = "l", lty = 3, lwd = 2)
axis(1)
axis(2)
legend("topright", legend = c("cglasso", "glasso", "MissGlasso"), lty = 1:3, box.col = 0, cex = 1.3, lwd = 2)
dev.off()

# thetah box-plot
data_boxplot <- cbind(apply(mse_a[, , "thetah", "rcglasso"], 1, min), apply(mse_a[, , "thetah", "glasso"], 1, min), apply(mse_a[, , "thetah", "mglasso"], 1, min))

pdf(file = "box-theta-model1-H25.pdf")
par(mai = c(1.02, 1.2, 0.82, 0.42),cex.axis = 1.3, cex.main = 1.5, cex.lab = 1.5)
boxplot(data_boxplot, outline = FALSE, notch = TRUE, axes = FALSE, col = "gray90", cex = 1.3, ylab = expression(paste(min[rho],"MSE(", hat(Theta^rho),")")), ylim = c(0, 140))
axis(2)
axis(1, col = 0, labels = c("cglasso", "glasso", "MissGlasso"), at = 1:3)
dev.off()

# muh box-plot
data_boxplot <- cbind(apply(mse_a[, , "muh", "rcglasso"], 1, min), apply(mse_a[, , "muh", "mglasso"], 1, min))
pdf(file = "box-mu-model1-H25.pdf")
par(mai = c(1.02, 1.2, 0.82, 0.42),cex.axis = 1.3, cex.main = 1.5, cex.lab = 1.5)
boxplot(data_boxplot, outline = FALSE, notch = TRUE, axes = FALSE, col = "gray90", cex = 1.3, ylab = expression(paste(min[rho],"MSE(", hat(mu^rho),")")))
axis(2)
axis(1, col = 0, labels = c("cglasso", "glasso", "MissGlasso"), at = 1:3)
dev.off()
