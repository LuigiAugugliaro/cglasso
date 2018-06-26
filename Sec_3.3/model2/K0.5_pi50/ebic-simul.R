############################################################################################
#
# Author:   Luigi Augugliaro
# e-mail:   luigi.augugliaro@unipa.it
# home:     http://dssm.unipa.it/augugliaro/
# data:     27-03-2018
#
# Description: R code used for the simulation study (model 2 with pi = 0.50)
# reported in "l1-Penalized Censored Gaussian Graphical Model"
# Section 3.3 (see also Section 1 in Supplementary Materials)
#
# Journal: Biostatistics


# libraries and R code

library("tmvtnorm")
library("MASS")
library("glasso")
source("rcglasso_v_3.0.3.R")
source("rcglasso.fitmu.R")
source("rcglasso.fitTheta.R")
source("rcglasso_mle.R")
source("Q.R")
source("loglik.R")


# general setting
zero <- 1e-6
nsim <- 100                     # number of simulations
n <- c(100, 200, 300, 400)      # sample size
tht_ij1 <- 0.2                  # value of the non-zero partial correlation coefficients
tht_ij2 <- 0.1                  # value of the non-zero partial correlation coefficients
nrho <- 100                     # number of tuning parameters
rhoratio <- 0.01                # parameter used to define minrho, i.e., minrho = rhoratio * maxrho
g <- c(0, 0.25, 0.5, 0.75, 1)   # parameter of the ebic criterion
K <- 0.5                        # p approx n^K
d <- 0.05                       # % of the p variables are censured
up <- 40                        # right censoreing value
prob_cens <- 0.50               # probability of censoring

# objects used to save the results
results <- array(0, dim = c(length(n), nsim, length(g), 2, 5),
                dimnames = list(n = n, nsim = NULL, g = g, type = c("ll", "Qval"),
                summary = c("npar", "TP", "FP", "TN", "FN")))

Qval <- ll <- npar <- array(0, dim = c(nrho, nsim, length(n)), dimnames = list(nrho = NULL, nsim = NULL, n = n))

conv <- matrix(0, nrow = nsim, ncol = length(n))

set.seed(123)

for(k in 1:length(n)){
    nk <- n[k]
    p <- round(nk^K)
    Tht <- diag(p)
    diag(Tht[1:(p-1), 2:p]) <- tht_ij1
    diag(Tht[2:p, 1:(p-1)]) <- tht_ij1
    diag(Tht[1:(p-2), 3:p]) <- tht_ij2
    diag(Tht[3:p, 1:(p-2)]) <- tht_ij2
    Sgm = solve(Tht)
    U <- upper.tri(Tht, diag = FALSE)
    A <- abs(Tht[U]) > 0
    notA <- !A
    S <- sort(sample(1:p, ceiling(p * d)))
    mu <- rep(34, p)
    tau <- qnorm(1 - prob_cens)
    mu[S] <- up - sqrt(diag(Sgm)[S]) * tau
    for(i in 1:nsim){
        X <- mvrnorm(nk, mu = mu, Sigma = Sgm)
		X[X > up] <- up
        # fitting rcglasso model
        out <- rcglasso(X, up, fitmean = TRUE, scale = FALSE, nrho = nrho, rhoratio = rhoratio, verbose = FALSE)
        if(out$conv == 0){
            for(m in 1:nrho){
                S_ini <- out$S[, , m]
                muh_ini <- out$muh[, m]
                thetah_ini <- out$thetah[, , m]
                model <- abs(thetah_ini) > zero
                out_mle <- rcglasso_mle(X = X, k = up, model = model, S_ini = S_ini,
                                        muh_ini = muh_ini, thetah_ini = thetah_ini,
                                        fitmean = TRUE, scale = FALSE, verbose = FALSE)
                if(out_mle$conv == 0){
                    tht_h <- out_mle$thetah
                    npar[m, i, k] <- sum(abs(tht_h[upper.tri(tht_h, diag = FALSE)]) > zero)
                    ll[m, i, k] <- loglik(out_mle)$loglik
                    Qval[m, i, k] <- Q(out_mle)$Qval
                } else {
                    conv[i, k] <- 1
                    break
                }
            }
            #######################################
            # eBIC section: computed with ll
            for(j in 1:5){
                ebic <- -2 * ll[, i, k] + npar[, i, k] * (log(nk) + 4 * g[j] * log(p))
                best <- which.min(ebic)
                tht_h_fit <- out$thetah[, , best]
                results[k, i, j, "ll", "npar"] <- npar[best, i, k]
                results[k, i, j, "ll", "TP"] <- sum(abs(tht_h_fit[U][A]) > zero)
                results[k, i, j, "ll", "FN"] <- sum(abs(tht_h_fit[U][A]) <= zero)
                results[k, i, j, "ll", "TN"] <- sum(abs(tht_h_fit[U][notA]) <= zero)
                results[k, i, j, "ll", "FP"] <- sum(abs(tht_h_fit[U][notA]) > zero)
            }
            #######################################
            # eBIC section: computed with Qval
            for(j in 1:5){
                ebic <- -2 * Qval[, i, k] + npar[, i, k] * (log(nk) + 4 * g[j] * log(p))
                best <- which.min(ebic)
                tht_h_fit <- out$thetah[, , best]
                results[k, i, j, "Qval", "npar"] <- npar[best, i, k]
                results[k, i, j, "Qval", "TP"] <- sum(abs(tht_h_fit[U][A]) > zero)
                results[k, i, j, "Qval", "FN"] <- sum(abs(tht_h_fit[U][A]) <= zero)
                results[k, i, j, "Qval", "TN"] <- sum(abs(tht_h_fit[U][notA]) <= zero)
                results[k, i, j, "Qval", "FP"] <- sum(abs(tht_h_fit[U][notA]) > zero)
            }
        } else conv[i, k] <- 1
        cat("simulation ", i, "with n =", nk, "completed\n")
        if(is.element(i, 10 * 1:10)) save.image("Sim_ebic_tmp.RData")
    }
}

save.image("Sim_ebic.RData")






