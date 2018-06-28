################################################################################################
# Author: Luigi Augugliaro

rcglasso_mle <- function(X, k, model, S_ini, muh_ini, thetah_ini, fitmean = FALSE, scale = TRUE, nstp = 500, eps = 1e-6, tol = 1e-5, verbose = FALSE){
    zero <- 1.0e-13
    
    this.call <- match.call()
    if(missing(S_ini) | missing(muh_ini) | missing(thetah_ini)) stop("muh_ini and thetah_ini are needed to fit the model")
    if(!is.logical(model)) stop("model is not a logical matrix")
    if(missing(k)) stop("k is missing")
    if(is.null(colnames(X)))
    colnames(X) <- paste("X", 1:dim(X)[2], sep = "")
    if(length(k) == 1) k <- rep(k, dim(X)[2])
    R <- matrix(FALSE, dim(X)[1], dim(X)[2])
    for(m in 1:dim(X)[2])   R[, m] <- X[, m] >= k[m]
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
    niter <- 0
    conv <- 0
    U <- upper.tri(thetah_ini, diag = FALSE)

    if(all(!model[upper.tri(model, diag = FALSE)])){
        i <- 1
        S_bar <- S_ini
        muh_new <- muh_ini
        thetah_new <- thetah_ini
    } else {
        rho <- model
        rho[model] <- 1e-3
        rho[!model] <- 1e+13
        muh_old <- muh_ini
        thetah_old <- thetah_ini
        for(i in 1:nstp){
            out <- rcglasso.fitTheta(R, nobs, Xsub, k, n_mis, row_mis, col_mis, rho, muh_old, thetah_old, nstp, tol, eps, verbose)
            if(is.element(out$conv, 1:2)){
                conv <- out$conv
                break
            }
            S_bar <- out$S_bar
            niter <- niter + out$niter
            thetah_new <- out$thetah
            if(fitmean) {
                ########################
                if(verbose) cat("\tfitting mu\n")
                out <- rcglasso.fitmu(T1_obs, nobs, Xsub, k, n_mis, row_mis, col_mis, rho, muh_old, thetah_new, nstp, tol, eps, verbose)
                if(is.element(out$conv, 1:2)){
                    conv <- out$conv
                    break
                }
                muh_new <- out$muh_new
                niter <- niter + out$niter
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
		if(is.element(out$conv, 1:2)) warning("rcglasso_mle does not converge")
    }
    if(i == nstp) {
        warning("maximum number of iterations has been reached (increase argument ", sQuote("nstp"), ")\n")
        conv <- 1
        break
    }
    
    muh <- muh_new
    thetah <- thetah_new
    S <- S_bar
    npar <- 2 * p + sum(abs(thetah[U]) > zero)
    if(fitmean) niter <- niter + i
    if(verbose) cat("\tn. iter =", niter, "conv =", conv, "\n\n")
    
    out <- list(muh = muh, thetah = thetah, S = S, npar = npar, nobs = nobs, X = X, k = k, model = model,
        R = R, niter = niter, conv = conv, call = this.call)
    class(out) <- "cglasso"
    out
}




































