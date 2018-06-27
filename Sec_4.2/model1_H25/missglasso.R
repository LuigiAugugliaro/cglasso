################################################################################################
# Author: Luigi Augugliaro
# data: 10-09-2017
# version: 1.0.0
#
# Description: 'missglasso_path' implements the algorithm proposed in Stadler & Buhlmann (2012) to fit
#               a glasso model with missing-at-random values.
#
# Arguments
# X: input matrix, of dimension n x p; each row is an observation vector.
# k: right censoring value.
# nrho: number of rho-values
# minrho: smallest rho-value
# nstp: number of iteration steps
# eps: number used of the convergence of the algorithm
# verbose: flag for printing out information as iterations proceed.

missglasso_path <- function(X, k, nrho = 50, minrho = 1e-7, nstp = 100, eps = 1e-4, 
                            maxrho = 1, verbose = FALSE, ...){
    if(missing(k)) stop("k is missing")
    p <- dim(X)[2]
    if(length(k) == 1) k <- rep(k, p)
    R <- t(apply(X, 1, function(x) x >= k))
    idsub <- apply(R, 1, all)
    if(any(idsub)){
        idsub <- !idsub
        Xsub <- X[idsub, , drop = FALSE]
        R <- R[idsub, , drop = FALSE]
    } else Xsub <- X
    n <- dim(Xsub)[1]
    row_mis <- which(apply(R, 1, any))
    n_mis <- length(row_mis)
    col_mis <- apply(R[row_mis, , drop = FALSE], 1, which)
    Y <- Xsub
    Y[R] <- 0
    T1_obs <- colSums(Y)
    T2_obs <- crossprod(Y)
    muh <- matrix(0, p, nrho)
    thetah <- array(0, dim = c(p, p, nrho))
    niter <- vector(mode = "numeric", length = nrho)
    conv <- vector(mode = "numeric", length = nrho)

    for(m in 1:p){
        rm <- R[, m]
        if(any(rm)){
            nm <- sum(rm)
            xobs <- Xsub[!rm, m, drop = FALSE]
            muh_old <- mean(xobs)
            sigmah_old <- mean(xobs^2) - muh_old^2
            if(verbose) cat("muh", "dmuh", "sigmah", "dsigmah", "\n")
            
            for(i in 1:nstp){
                tmean <- muh_old
                tvar <- sigmah_old
                Y[rm, m] <- tmean
                T1_bar <- T1_obs[m] + nm * tmean
                T2_bar <- T2_obs[m, m] + nm * (tvar + tmean^2)
                muh_new <- T1_bar / n
                sigmah_new <- T2_bar / n - muh_new^2
                dmuh <- abs(muh_new - muh_old)
                dsigmah <- abs(sigmah_new - sigmah_old)
                if(verbose) cat(c(muh_new, dmuh, sigmah_new, dsigmah), "\n")
                if(max(dmuh, dsigmah) < eps) break
                else {
                    muh_old <- muh_new
                    sigmah_old <- sigmah_new
                }
            }
            
            if(i == nstp) stop("maximum number of iterations has been reached at the starting point (increase ''nstp'')")
            if(verbose) cat("n. iter =", i, "\n\n")
            
            muh[m, 1] <- muh_new
            thetah[m, m, 1] <- 1 / sigmah_new
        } else {
            muh[m, 1] <- T1_obs[m] / n
            thetah[m, m, 1] <- 1 / (T2_obs[m, m] / n - muh[m, 1]^2)
        }
        
    }
    
    S_bar <- crossprod(Y) / n - outer(muh[, 1], muh[, 1])
    if(minrho >= maxrho) stop("ridurre minrho; massimo valore di rho = ", maxrho)
    rho_seq <- seq(from = maxrho, to = minrho, length = nrho)
    
    for(h in 2:nrho){
        rho <- rho_seq[h]
        muh_old <- muh[, h - 1]
        thetah_old <- thetah[, , h - 1]

        if(verbose) {
            cat("Model", h, "rho = ", round(rho, 4), "\n")
            cat(c("dmuh", "dthetah", "\n"))
        }
        for(i in 1:nstp){
            T_bar <- Tmaker_missglasso(T1_obs, T2_obs, Xsub, n_mis, row_mis, col_mis, muh_old, thetah_old)
            muh_new <- T_bar$T1_bar / n
            S_bar <- T_bar$T2_bar / n - outer(muh_new, muh_new)
            out <- glasso(S_bar, rho = rho, penalize.diagonal = FALSE, wi.init = thetah_old, ...)
            thetah_new <- (out$wi + t(out$wi)) / 2
            dmuh <- max(abs(muh_new - muh_old))
            dthetah <- max(abs(thetah_new - thetah_old))
            if(verbose) cat(c(dmuh, dthetah), "\n")
            # controllo convergenza
            if(max(dmuh, dthetah) < eps) break
            else{
                muh_old <- muh_new
                thetah_old <- thetah_new
            }
        }
        
        if(i == nstp) {
            warning("model number ", h, ": maximum number of iterations has been reached (increase ''nstp'')")
            conv[h] <- 1
        }
        
        if(verbose) cat("n. iter =", i, "conv =", conv[h], "\n\n")
        
        muh[, h] <- muh_new
        thetah[, , h] <- thetah_new
        niter[h] <- i
        
    }
    
    out <- list(muh = muh, thetah = thetah, rho_seq = rho_seq, X = X, k = k, R = R, niter = niter, conv = conv)
    out
    
}

Tmaker_missglasso <- function(T1_obs, T2_obs, X, n_mis, row_mis, col_mis, muh, Thetah){
    T1_bar <- T1_obs
    T2_bar <- T2_obs
    for(h in 1:n_mis){
        i <- row_mis[h]
        ri <- col_mis[[h]]
        x_obs <- X[i, -ri]
        muh_obs <- muh[-ri]
        muh_mis <- muh[ri]
        Thetah_obs <- Thetah[-ri, -ri, drop = FALSE]
        Thetah_mis <- Thetah[ri, ri, drop = FALSE]
        Thetah_misobs <- Thetah[ri, -ri, drop = FALSE]
        zeta <- drop(Thetah_misobs %*% (x_obs - muh_obs))
        Sigmah_mis_given_obs <- drop(solve(Thetah_mis))
        muh_mis_given_obs <- muh_mis - drop(Sigmah_mis_given_obs %*% zeta)
        T1_bar[ri] <- T1_bar[ri] + muh_mis_given_obs
        T2_bar[-ri, ri] <- T2_bar[-ri, ri, drop = FALSE] + outer(x_obs, muh_mis_given_obs)
        T2_bar[ri, -ri] <- t(T2_bar[-ri, ri, drop = FALSE])
        T2_bar[ri, ri] <- T2_bar[ri, ri] + Sigmah_mis_given_obs + outer(muh_mis_given_obs, muh_mis_given_obs)
    }
    out <- list(T1_bar = T1_bar, T2_bar = T2_bar)
    out
}
