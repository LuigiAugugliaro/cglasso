############################################################################################
#
# Author:   Antonino Abbruzzo
# e-mail:   antonino.abbruzzo@unipa.it
# data:     04-10-2018
#
# Description: R code used for the simulation study
# reported in "l1-Penalized Censored Gaussian Graphical Model"
# Section 4.3
#
# Journal: Biostatistics

# Library and R code
library("glasso")
library("MASS")
library("tmvtnorm")
library("imputeLCMD")
library("missForest")
library("VIM")

source("rcglasso_v_3.0.3.R")
source("rcglasso.fitmu.R")
source("rcglasso.fitTheta.R")
source("loglik.glasso.R")
source("AIC.cglasso.R")
source("Q.R")
source("impute_rcglasso.R")
source("missglasso.R")
source("impute_missglasso.R")
source("prcurve.R")

#############################################
# reading dataset
DataObs <- read.table("Arabinobis.txt")
# we work on scaled data
XObs <- scale(DataObs)

#############################################
# output
mse.tbl <- auc.tbl <- matrix(0, nrow = 3, ncol = 5)
colnames(auc.tbl) <- colnames(mse.tbl) <- c("cglasso", "imputedLMCD", "MissGlasso", "missknn", "missForest")

#############################################
# we fit glasso model on the complete dataset
#
# setting
n <- dim(XObs)[1]
p <- dim(XObs)[2]
minrho <- 0.107
maxrho <- 1
nrho <- 200
rholist <- seq(from = minrho, to = maxrho, length.out = nrho)
S <- cov(XObs)

# fitting glasso model
out_glasso <- glassopath(S, penalize.diagonal = FALSE, rholist = rholist)

# in our simulation study the sparse precision matrix estimated using the
# BIC measure is treated as true precision matrix
out_ll <- loglik_func(out_glasso, S, n = nrow(XObs), gamma = 0)
ThtTrue <- out_glasso$wi[, , which.min(out_ll$bic)]

######################################################################
# we simulate three datasets at different
# levels of right censoring, by recording the
# 10%, 20% and 30% of the highest values, respectively.

prob_cens <- c(0.1, 0.2, 0.3) # probability of a censored data

for(m in 1:3){
    #################################
    # preparing datasets
    #################################
    cat("Starting simulation with prob_cens =", prob_cens[m], "...")
    up <- quantile(XObs, probs = 1 - prob_cens[m]) # right censoring values
    Id <- XObs >= up
    Xcens <- XObs
    Xcens[Id] <- up # dataset with censored data
    Xmiss <- XObs
    Xmiss[Id] <- NA # dataset with missing data

    #################################
    # rcglasso section
    #################################
    out_rcglasso <- rcglasso(X = Xcens, k = up, fitmean = TRUE, scale = FALSE, nrho = nrho,
                        minrho = minrho, maxrho = maxrho, verbose = TRUE, eps = 1e-5)
    out_aic <- AIC(out_rcglasso, "Q", k = "bic")
    id_opt <- which.min(out_aic$gof)
    input <- list(X = out_rcglasso$X, k = out_rcglasso$k, R = out_rcglasso$R, muh = out_rcglasso$muh[, id_opt,
                    drop = FALSE], thetah = out_rcglasso$thetah[, , id_opt, drop = FALSE],
                    rho_seq = out_rcglasso$rho_seq[id_opt])
    out_impu_rcglasso <- impute_rcglasso(input)
    mse.tbl[m, "cglasso"] <-  sqrt(sum((out_impu_rcglasso[, , 1][Id] - XObs[Id])^2))
    auc.tbl[m, "cglasso"] <- prcurve(out_rcglasso$thetah, ThtTrue, rho = out_rcglasso$rho_seq)$auc
    
    #################################
    # imputedLMCD section
    #################################
    Xbar <- matrix(0, nrow = n, ncol = p)
    for(i in 1:1000)
        Xbar <- Xbar - impute.MAR.MNAR(-Xmiss, model.selector = rep(0, p))
    Xbar <- Xbar / 1000
    rholist <- seq(from = 0, to = maxrho, length.out = nrho)
    out_glassopath <- glassopath(cov(Xbar), penalize.diagonal = FALSE, rholist = rholist)
    mse.tbl[m, "imputedLMCD"] <- sqrt(sum((Xbar[Id] - XObs[Id])^2))
    auc.tbl[m, "imputedLMCD"] <- prcurve(out_glassopath$wi, ThtTrue, rho = rholist)$auc
    
    #################################
    # MissGlasso section
    #################################
    out_missglasso <- missglasso_path(X = Xcens, k = up, verbose = TRUE, eps = 1e-4, nrho = 1000,
                        minrho = 0.01, maxrho = 0.40, nstp = 800)
    mse_missglasso <- vector(mode = "numeric", length = 1000)
    for(i in 1:1000){
        input <- list(X = out_missglasso$X, k = out_missglasso$k, R = out_missglasso$R,
                        muh = out_missglasso$muh[, i, drop = FALSE],
                        thetah = out_missglasso$thetah[, , i, drop = FALSE],
                        rho_seq = out_missglasso$rho_seq[i])
        out_impu_missglasso <- impute_missglasso(input)
        mse_missglasso[i] <- sqrt(sum((out_impu_missglasso[ , , 1][Id] - XObs[Id])^2))
    }
    mse.tbl[m, "MissGlasso"] <- mse_missglasso[id_opt]
    auc.tbl[m, "MissGlasso"] <- prcurve(out_missglasso$thetah, ThtTrue, rho = out_missglasso$rho_seq)$auc
    
    #################################
    # missknn section
    #################################
    out_knn <- kNN(Xmiss)
    Xknn <- as.matrix(out_knn[, 1:39])
    out_glassopath <- glassopath(cov(Xknn), penalize.diagonal = FALSE, rholist = rholist)
    mse.tbl[m, "missknn"] <- sqrt(sum((Xknn[Id] - XObs[Id])^2))
    auc.tbl[m, "missknn"] <- prcurve(out_glassopath$wi, ThtTrue, rho = rholist)$auc
    
    #################################
    # missForest section
    #################################
    out_missForest <- missForest(Xmiss)
    out_glassopath <- glassopath(cov(out_missForest$ximp), penalize.diagonal = FALSE, rholist = rholist)
    mse.tbl[m, "missForest"] <- sqrt(sum((out_missForest$ximp[Id] - XObs[Id])^2))
    auc.tbl[m, "missForest"] <- prcurve(out_glassopath$wi, ThtTrue, rho = rholist)$auc
    
    cat("completed!\n")
}

###################################
# Table 3 in Main Document
###################################
round(rbind(mse.tbl, auc.tbl), 2)




























