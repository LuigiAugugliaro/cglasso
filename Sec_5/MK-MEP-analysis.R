############################################################################################
#
# Author:   Luigi Augugliaro
# e-mail:   luigi.augugliaro@unipa.it
# home:     http://dssm.unipa.it/augugliaro/
# data:     05-07-2018
#
# Description: R code used for the analysis of the real dataset
# describted in "l1-Penalized Censored Gaussian Graphical Model",
# Section 5
#
# Journal: Biostatistics

#########################################################################
# library + R code
library(glasso)
library(igraph)
library(latex2exp)

source("rcglasso_v_3.0.4.R")
source("rcglasso.fitmu.R")
source("rcglasso.fitTheta.R")
source("rcglasso_mle.R")
source("AIC.cglasso.R")
source("loglik.cglasso.R")
source("Q.R")
source("Qmle.R")

######################################
# Reading Normalized Ct-values

NormCt <- read.table("MK-MEP.txt")
dim(NormCt)
# 48 87

# right censoring value
up <- 40

# filtering out unobserved genes
id <- which(apply(NormCt >= up, 2, sum) == dim(NormCt)[1])
subNormCt <- NormCt[, -id]
dim(subNormCt)
# 48 76

######################################
# Fig 2 (a) in Main Document

y <- t(apply(subNormCt, 2, function(x) c(sum(x >= up), sum(x < up))))
colnames(y) <- c("Miss", "Obs")
x <- apply(subNormCt, 2, function(x) mean(x[x < up]))
mis.prop <- y[, 1] / sum(y[1, ])
sort(mis.prop)
out.glm <- glm(y ~ x, family = binomial)

pdf(file = "propmissVSmeandCt.pdf")
par(mai = c(1.02, 1.02, 0.82, 0.42))
plot(x, mis.prop, pch = 20, xlab = "Mean normalized cycle-threshold", ylab = "Proportion of right-censored data", axes = TRUE, cex.lab = 1.7, cex.axis = 1.2)
points(x[order(x)], out.glm$fitted.values[order(x)], type = "l", lwd = 2, lty = 1)
abline(h = 0.85, lty = 2, lwd = 2)
dev.off()

# filtering out genes with a proporsion of censoring about 85%
id <- apply(subNormCt >= up, 2, mean) < 0.85
X <- subNormCt[, id]
X <- as.matrix(X)
dim(X)
# 48 63

######################################
# Fitting rcglasso model

# setting
nrho <- 200
minrho <- 0.85
eps <- 1.0e-5

out_rcglasso <- rcglasso(X = X, k = up, fitmean = TRUE, scale = FALSE, nrho = nrho, minrho = minrho, verbose = TRUE, eps = 1e-5)

# computing EBIC criterion

muh_ini <- thetah_ini <- 0
conv <- npar <- Qval <- vector(length = nrho)

for(m in 1:nrho){
    
    S_ini <- out_rcglasso$S[, , m]
    model <- abs(out_rcglasso$thetah[, , m]) > 1.0e-6
    if(muh_ini[1] == 0) muh_ini <- out_rcglasso$muh[, m]
    if(!is.matrix(thetah_ini)) thetah_ini <- out_rcglasso$thetah[, , m]
    
    out_mle <- rcglasso_mle(X = X, k = up, model = model, S_ini = S_ini, muh_ini = muh_ini, thetah_ini = thetah_ini, fitmean = TRUE, scale = FALSE, verbose = TRUE, eps = 1e-5)
    
    if(out_mle$conv == 0){
        tht_h <- out_mle$thetah
        npar[m] <- sum(abs(tht_h[upper.tri(tht_h, diag = FALSE)]) > 1.0e-6)
        Qval[m] <- Qmle(out_mle)$Qval
        muh_ini <- out_mle$muh
        thetah_ini <- tht_h
    } else {
        thetah_ini <- 0
        conv[m] <- 1
        break
    }
    
    cat("m =", m, "\n")
}

n <- dim(X)[1]  # sample size
p <- dim(X)[2]  # number of variables
g <- 0.5        # parameter of the EBIC criterion

ebic <- -2 * Qval + npar * (log(n) + 4 * g * log(p))

######################################
# Fig 2 (b) in Main Document

subrho <- seq(2, 200, length = 100)

pdf(file = "ebicpath.pdf")
par(mai = c(1.02, 1.02, 0.82, 0.42))
plot(out_rcglasso$rho[subrho], ebic[subrho], type = "l", xlab = expression(rho), ylab = TeX('$\\bar{BIC}_{0.5}(\\hat{E}^{\\rho})$'), lty = 1, lwd = 2, cex.lab = 1.7, cex.axis = 1.2)
abline(v = out_rcglasso$rho[subrho][which.min(ebic[subrho])], lwd = 2, lty = 2)
dev.off()

######################################
# Fig 2 (c) in Main Document

th.opt <- out_rcglasso$thetah[ , , which.min(ebic)]
adj <- 1 * (th.opt != 0)
diag(adj) <- 0
dim(adj)
# 63 63

mean(adj[upper.tri(adj)])
# 0.006144393

# filtering out isolated vertices
adj <- adj[apply(adj, 1, sum) != 0, apply(adj, 2, sum) != 0]
dim(adj)
# 12 12

# renames genes
rownames(adj)[c(2, 4, 5, 6)]
#[1] "CD117/KIT"  "CD42/GP1BA" "CD61/ITGB3" "CD71/TFRC"
rownames(adj)[c(2, 4, 5, 6)] <- colnames(adj)[c(2, 4, 5, 6)] <- c("CD117", "CD42", "CD61", "CD71")

net <- graph.adjacency(adj, mode = "undirected")
net.layout <- layout_with_lgl(net)

pdf(file = "mk_mep_network.pdf")
V(net)$color <- "white"
V(net)$size <- 1
V(net)$frame.color <- NA
V(net)$label.color <- "black"
V(net)$label.cex <- 1.8
E(net)$color <- "gray75"
plot(net, layout = net.layout)
mtext(text = "Megakaryocytic MEP population", cex = 2.5, line = 1)
dev.off()
