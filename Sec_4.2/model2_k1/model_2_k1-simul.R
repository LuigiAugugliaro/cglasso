############################################################################################
#
# Author:   Luigi Augugliaro and Antonino Abbruzzo
# e-mail:   luigi.augugliaro@unipa.it
# e-mail:   antonino.abbruzzo@unipa.it
# home:     http://dssm.unipa.it/augugliaro/
# data:     04-09-2018
#
# Description: R code used for the simulation study
# reported in "l1-Penalized Censored Gaussian Graphical Model"
# Section 4.2. Model 2 with k = 1
#
# Journal: Biostatistics


###########################################
# libraries and R code
library("MASS")
library("glasso")
library("tmvtnorm")
library("norm")
library("huge")
library("xtable")

source("rcglasso_path_v_1.1.1.R")
source("rT_maker_approx.R")
source("rT_maker_exact.R")
source("missglasso.R")
source("prcurve.R")
###########################################

# setting
nsim <- 100         # number of simulations
up <- 40            # right censoring value
n <- 100            # sample size
p <- 50             # number of variables
k <- 1              # parameter used to control the parsity in Theta
prob_link <- k/p    # probability that a link will be present
H <- 30             # maximum number of variable that may be censored
nrho <- 30          # number of rho-values

# output
conv <- rep(TRUE, nsim)
rho_rcglasso <- matrix(0, nsim, nrho)
rho_glasso <- matrix(0, nsim, nrho)
rho_mglasso <- matrix(0, nsim, nrho)
muh_rcglasso <- array(0, dim = c(nsim, p, nrho))
muh_mglasso <- array(0, dim = c(nsim, p, nrho))
thetah_rcglasso <- array(0, dim = c(nsim, p, p, nrho))
thetah_glasso <- array(0, dim = c(nsim, p, p, nrho))
thetah_mglasso <- array(0, dim = c(nsim, p, p, nrho))

################################################
# Starting Simulation study

set.seed(123)

mu = c(rep(40, H), runif(p - H, 10, 35))
sim <- huge.generator(n = n, d = p, graph = "random", prob = prob_link)
ThetaTRUE <- round(sim$omega, 5)
SigmaTRUE <- round(sim$sigma, 5)
rm(sim)

lambda.min.ratio <- 0.05
lambda.max <-  max(abs(SigmaTRUE[upper.tri(SigmaTRUE)]))
lambda.min = lambda.min.ratio * lambda.max

for(j in 1:nsim){
    X <- mvrnorm(n, mu = mu, Sigma = SigmaTRUE)
    X[X > up] <- up
  
    #############################################################
    # right censored glasso
    #############################################################
  
    out_rcglasso <- rcglasso_path(X = X, k = up, verbose = FALSE, nrho = nrho, minrho = lambda.min, maxrho = lambda.max, nstp = 1000, eps = 1e-4)
  
    conv[j] <- sum(out_rcglasso$conv) == 0
    muh_rcglasso[j, , ] <- out_rcglasso$muh
    thetah_rcglasso[j, , , ] <- out_rcglasso$thetah
    rho_rcglasso[j, ] <- out_rcglasso$rho

    #############################################################
    # glasso
    #############################################################
    
    S <- cov(X)
    rho_seq <- seq(from = 0.02, to = max(abs(S[upper.tri(S)])), length = nrho)
    out_glasso <- glassopath(S, rholist = rho_seq, penalize.diagonal = FALSE, trace = 0)
    thetah_glasso[j, , , ] <- out_glasso$wi[, , nrho:1]
    rho_glasso[j, ] <- rho_seq

    #############################################################
    # missglasso
    #############################################################
    
    out_missglasso <- missglasso_path(X = X, k = up, verbose = FALSE, nrho = nrho, minrho = 0.005, eps = 1e-4, maxrho = lambda.max, nstp = 1000)
    conv[j] <- conv[j] & sum(out_missglasso$conv) == 0
    
    muh_mglasso[j, , ] <- out_missglasso$muh
    thetah_mglasso[j, , , ] <- out_missglasso$thetah
    rho_mglasso[j, ] <- out_missglasso$rho

    cat("Simulation", j, "completed\n")
    if(j %in% (10 * (1:10))) save.image("model_2_k1_tmp.RData")
}

save.image("model_2_k1.RData")
