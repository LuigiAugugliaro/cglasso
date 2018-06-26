############################################################################################
#
# Author:   Luigi Augugliaro
# e-mail:   luigi.augugliaro@unipa.it
# home:     http://dssm.unipa.it/augugliaro/
# data:     11-30-2017
#
# Description: R code used for the simulation study
# reported in "l1-Penalized Censored Gaussian Graphical Model"
# Section 4.1
#
# Journal: Biostatistics

library("tmvtnorm")
library("MASS")
library("glasso")
library("huge")
source("rcglasso_path_v_1.1.0.R")
source("rT_maker_approx.R")
source("rT_maker_exact.R")

# general setting
nsim <- 100         # number of simulations
n <- 100            # sample size
p <- 10             # number of variables
d <- 2:8            # number of censored variables
k <- 40             # right censoring value
prob_cens <- 0.25   # probability of censoring
nrho <- 30          # number of point of the solution curves
minrho <- 0.001     # smallest rho-value
zero <- 1e-7        # threshold used to define a zero value

# objects used to save the results
time_table <- array(data = 0, dim = c(nsim, 2, length(d)), dimnames = list(nsim = NULL, method = c("exact", "approx"), d = d))
muh_f <- array(data = 0, dim = c(nrho, nsim, length(d)), dimnames = list(rho = NULL, nsim = NULL, d = d))
thetah_f <- array(data = 0, dim = c(nrho, nsim, length(d)), dimnames = list(rho = NULL, nsim = NULL, d = d))

# starting simulation study
for(h in 1:length(d)){
    
    out_huge <- huge.generator(n = n, d = p, graph = "random", prob = 0.10)
    Sgm <- out_huge$sigma
    Tht <- solve(Sgm)
    U <- upper.tri(Tht, diag = TRUE)
    A <- sort(sample(1:p, d[h]))
    mu <- rep(34, p)
    tau <- qnorm(1 - prob_cens)
    mu[A] <- k - sqrt(diag(Sgm)[A]) * tau
    
    for(i in 1:nsim){
        X <- mvrnorm(n, mu = mu, Sigma = Sgm)
        X[X > k] <- k
        
        # exact method
        time_table[i, "exact", h] <- system.time(out_exact <- rcglasso_path(X, k = k, method = "exact", verbose = FALSE, minrho = minrho, nrho = nrho))[3]
        
        # approximated method method
        time_table[i, "approx", h] <- system.time(out_approx <- rcglasso_path(X, k = k, method = "approx", verbose = FALSE, minrho = minrho, nrho = nrho))[3]
        
        muh_f[, i, h] <- apply(out_exact$muh - out_approx$muh, 2, crossprod)
        thetah_f[, i, h] <- apply(thetah_e - thetah_a, 2, crossprod)
    
        if(is.element(i, 10 * 1:10))
            save.image("CompAsp_tmp.RData")
            
        cat("simulation ", i, "with d =", d[h], "completed.\n")
        
    }
    save.image("CompAsp_tmp.RData")
}
save.image("CompAsp.RData")
