################################################################################################
# Author: Luigi Augugliaro
# data: 05-11-2017
# version: 3.0.3

rcglasso <- function(X, k, fitmean = FALSE, scale = TRUE, nrho = 50, maxrho, minrho = 1e-3, nstp = 500, eps = 1e-6, tol = 1e-5, verbose = FALSE){
    this.call <- match.call()
    if(missing(k)) stop("k is missing")
    if(nrho <= 0) stop("argument ", sQuote("nrho"), " is not a posotive value")
    zero <- 1.0e-13
    if(is.null(colnames(X)))
    colnames(X) <- paste("X", 1:dim(X)[2], sep = "")
    if(length(k) == 1) k <- rep(k, dim(X)[2])
    R <- matrix(FALSE, dim(X)[1], dim(X)[2])
    for(m in 1:dim(X)[2])   R[, m] <- X[, m] >= k[m]
    if(sum(R) == 0) stop("dataset without censored data: use glasso to fit the model")
    xv <- vector(mode = "numeric", length = dim(X)[2])
    xm <- vector(mode = "numeric", length = dim(X)[2])
    for(m in 1:dim(X)[2]) {
        xm[m] <- mean(X[!R[, m], m])
        xv[m] <- mean(X[!R[, m], m]^2) - xm[m]^2
    }
    id <- which(xv <= zero)
    if(length(id) > 0){
        msg <- paste("column", sQuote(id), "is removed because the standard deviations is zero", collapse = "\n")
        warning(msg)
        X <- X[, -id]
        R <- R[, -id]
        k <- k[-id]
        xm <- xm[-id]
        xv <- xv[-id]
    }
    p <- dim(X)[2]
    idsub <- apply(R, 1, all)
    if(any(idsub)){
        idsub <- !idsub
        Xsub <- X[idsub, , drop = FALSE]
        R <- R[idsub, , drop = FALSE]
    } else Xsub <- X
    
    if(scale){
        Xsub <- scale(Xsub, center = FALSE, scale = sqrt(xv))
        k <- k / sqrt(xv)
        xm <- xm / sqrt(xv)
    }
    
    nobs <- dim(Xsub)[1]
    row_mis <- which(apply(R, 1, any))
    n_mis <- length(row_mis)
    col_mis <- apply(R[row_mis, , drop = FALSE], 1, which)
    Y <- Xsub
    Y[R] <- 0
    T1_obs <- colSums(Y)
    T2_obs <- crossprod(Y)
    U <- upper.tri(T2_obs)
    
    cat("Fitting marginal models\n")
    muh_old <- vector(mode = "numeric", length = p)
    thetah_old <- matrix(0, p, p)
    if(!fitmean) dxm <- vector(mode = "numeric", length = dim(X)[2])
    for(m in 1:p){
        if(verbose) cat("\tFitting marginal model", m, "\n")
        rm <- R[, m]
        if(any(rm)){
            nm <- sum(rm)
            xobs <- Xsub[!rm, m, drop = FALSE]
            mumh_old <- mean(xobs)
            varmh_old <- mean(xobs^2) - mumh_old^2
            if(verbose) {
                cat("\tmuh", "dmuh", "sigmah", "dsigmah", "\n")
                cat(c("\t", mumh_old, 0, varmh_old, 0, "\n"))
            }
            for(i in 1:nstp){
                z <- (k[m] - mumh_old) / sqrt(varmh_old)
                ptail <- pnorm(z, lower.tail = FALSE)
                if(ptail >= tol) {
                    ratio <- dnorm(z) / ptail
                    tmean <- mumh_old + sqrt(varmh_old) * ratio
                    tvar <- varmh_old * (1 + z * ratio - ratio^2)
                } else {
                    tmean <- mumh_old
                    tvar <- varmh_old
                }
                T1_bar <- T1_obs[m] + nm * tmean
                T2_bar <- T2_obs[m, m] + nm * (tvar + tmean^2)
                Y[rm, m] <- tmean
                if(fitmean){
                    mumh_new <- T1_bar / nobs
                    varmh_new <- T2_bar / nobs - mumh_new^2
                } else {
                    mumh_new <- mumh_old
                    dxm[m] <- xm[m]^2 - 2 * xm[m] * T1_bar / nobs
                    varmh_new <- T2_bar / nobs + dxm[m]
                }
                dmu <- abs(mumh_new - mumh_old)
                dvar <- abs(varmh_new - varmh_old)
                if(verbose) cat(c(mumh_new, dmu, varmh_new, dvar), "\n")
                if(max(dmu, dvar) < eps) break
                else {
                    mumh_old <- mumh_new
                    varmh_old <- varmh_new
                }
            }
            if(i == nstp)
            stop("model number 1: maximum number of iterations has been reached (increase argument ", sQuote("nstp"), ")")
            if(verbose) cat("n. iter =", i, "\n\n")
            muh_old[m] <- mumh_new
            thetah_old[m, m] <- 1 / varmh_new
        } else {
        	muh_old[m] <- T1_obs[m] / nobs
        	thetah_old[m, m] <- 1 / (T2_obs[m, m] / nobs - muh_old[m]^2)
        }
    }
    
    if(fitmean) S_bar <- crossprod(Y) / nobs - outer(muh_old, muh_old)
    else {
        Y <- scale (Y, center = xm, scale = FALSE)
        S_bar <- crossprod(Y) / nobs
    }
    diag(S_bar) <- 1 / diag(thetah_old)
    if(missing(maxrho)) maxrho <- max(abs(S_bar[U]))
    
    nmodels <- nrho
    if(minrho > maxrho){
        warning("argument ", sQuote("minrho"), " is larger than ", sQuote("maxrho"), "(", round(maxrho, 3), "); solution path contains only one solution.")
        maxrho <- minrho
        nmodels <- 1
    }
    muh <- matrix(0, p, nmodels, dimnames = list(colnames(X), NULL))
    thetah <- array(0, dim = c(p, p, nmodels), dimnames = list(colnames(X), colnames(X), NULL))
    S <- array(0, dim = c(p, p, nmodels), dimnames = list(colnames(X), colnames(X), NULL))
    npar <- rep(2 * p, nmodels)
    niter <- vector(mode = "numeric", length = nmodels)
    conv <- 0

    rho_seq <- exp(seq(from = log(maxrho), to = log(minrho), length = nrho))
    for(h in 1:nmodels){
        rho <- rho_seq[h]
        cat("Fitting rcglasso model number", h,"with rho =", rho, "\n")
        for(i in 1:nstp){
            ########################
            if(verbose) cat("\tfitting theta\n")
            out <- rcglasso.fitTheta(R, nobs, Xsub, k, n_mis, row_mis, col_mis, rho, muh_old, thetah_old, nstp, tol, eps, verbose)
            if(is.element(out$conv, 1:2)) break
            S_bar <- out$S_bar
            niter[h] <- niter[h] + out$niter
            thetah_new <- out$thetah
            if(fitmean) {
                ########################
                if(verbose) cat("\tfitting mu\n")
                out <- rcglasso.fitmu(T1_obs, nobs, Xsub, k, n_mis, row_mis, col_mis, rho, muh_old, thetah_new, nstp, tol, eps, verbose)
                if(is.element(out$conv, 1:2)) break
                muh_new <- out$muh_new
                niter[h] <- niter[h] + out$niter
                dmuh <- max(abs(muh_new - muh_old))
                dthetah <- max(abs(thetah_new - thetah_old))
                if(max(dmuh, dthetah) < eps) break
                else {
                    muh_old <- muh_new
                    thetah_old <- thetah_new
                }
            } else {
                muh_new <- xm
                break
            }
        }
        if(i == nstp) {
            warning("maximum number of iterations has been reached (increase argument ", sQuote("nstp"), ")\n")
            conv <- 1
        }
        if(is.element(out$conv, 1:2)) {
            if(nmodels == 1) stop("Impossible to fit the rcglasso model at rho = ", rho)
            muh[, 1:(h - 1), drop = FALSE]
            thetah[, , 1:(h - 1), drop = FALSE]
            S <- S[, , 1:(h - 1), drop = FALSE]
            npar <- npar[1:(h - 1)]
            rho_seq <- rho_seq[1:(h - 1)]
            niter <- niter[1:(h - 1)]
            break
        }
        muh_old <- muh[, h] <- muh_new
        thetah_old <- thetah[, , h] <- thetah_new
        S[, , h] <- S_bar
        npar[h] <- npar[h] + sum(abs(thetah[, , h][U]) > zero)
        if(fitmean) niter[h] <- niter[h] + i
        if(verbose) cat("\tn. iter =", niter[h], "conv =", conv, "\n\n")
    }

    out <- list(muh = muh, thetah = thetah, S = S, npar = npar, rho_seq = rho_seq, nobs = nobs, X = X, k = k, R = R, niter = niter, conv = conv, call = this.call)
    class(out) <- "cglasso"
    out
}




































