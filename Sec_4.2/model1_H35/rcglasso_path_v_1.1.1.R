################################################################################################
# Author: Luigi Augugliaro
# data: 05-11-2017
# version: 1.1.1

# Arguments
# X = input matrix, of dimension n x p; each row is an observation vector.
# k = right censoring value; it can be a vector of length p specifying the right censoring value for each variable
# method = string used to specify how to compute the matrix S (see section 3.2). If method = "exact", the matrix S is
#           copmuted exact method, otherwise (method = "approximated") S is approximated using expression (3.8) in page 11
# nrho = number of rho-values
# minrho = smallest rho-value
# maxrho = largest rho-value
# nstp = number of iteration steps
# eps = argument used for convergence
# verbose = flag for printing out information as iterations proceed.

rcglasso_path <- function(X, k, method = c("approximated", "exact"), nrho = 50, minrho = 1e-7, maxrho = 1, nstp = 100, eps = 1e-4, verbose = FALSE, ...){
    method <- match.arg(method)
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
            out_mle <- mle.tmvnorm(xobs, upper = k[m], start = list(mu = mean(xobs), sigma = var(xobs)))
            muh_old <- out_mle@coef[1]
            sigmah_old <- out_mle@coef[2]
            if(verbose) cat("muh", "dmuh", "sigmah", "dsigmah", "\n")
        
            for(i in 1:nstp){
                out <- mtmvnorm(mean = muh_old, sigma = sigmah_old, lower = k[m])
                Y[rm, m] <- out$tmean
                T1_bar <- T1_obs[m] + nm * out$tmean
                T2_bar <- T2_obs[m, m] + nm * drop(out$tvar + out$tmean^2)
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
            
            if(i == nstp) 
            	stop("maximum number of iterations has been reached at the starting point (increase ''nstp'')")
            if(verbose) 
            	cat("n. iter =", i, "\n\n")
            
            muh[m, 1] <- muh_new
            thetah[m, m, 1] <- 1 / sigmah_new
        } else {
            muh[m, 1] <- T1_obs[m] / n
            thetah[m, m, 1] <- 1 / (T2_obs[m, m] / n - muh[m, 1]^2)
        }
        
    }
    
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
            T_bar <- if(method == "exact")
                        rT_maker_exact(T1_obs, T2_obs, k, Xsub, n_mis, row_mis, col_mis, muh_old, thetah_old)
                    else
                        rT_maker_approx(T1_obs, T2_obs, k, Xsub, n_mis, row_mis, col_mis, muh_old, thetah_old)
            muh_new <- T_bar$T1_bar / n
            S_bar <- T_bar$T2_bar / n - outer(muh_new, muh_new)
            out <- glasso(S_bar, rho = rho, penalize.diagonal = FALSE, wi.init = thetah_old, ...)
            thetah_new <- (out$wi + t(out$wi)) / 2
            dmuh <- max(abs(muh_new - muh_old))
            dthetah <- max(abs(thetah_new - thetah_old))
            if(verbose) cat(c(dmuh, dthetah), "\n")
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
