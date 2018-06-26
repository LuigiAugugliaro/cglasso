loglik <- function(object){
    #    require(mvtnorm)
	X <- object$X
	k <- object$k
	R <- object$R
    nobs <- object$nobs
	idsub <- apply(R, 1, all)
    if(any(idsub)){
        idsub <- !idsub
        Xsub <- X[idsub, , drop = FALSE]
        R <- R[idsub, , drop = FALSE]
    } else Xsub <- X
    row_mis <- which(apply(R, 1, any))
    row_obs <- setdiff(1:nobs, row_mis)
    if(length(row_obs) > 0) {
        oXsub <- X[row_obs, , drop = FALSE]
        n_obs <- length(row_obs)
    }
    n_mis <- length(row_mis)
    col_mis <- apply(R[row_mis, , drop = FALSE], 1, which)
    
	muh <- object$muh
	thetah <- object$thetah
    nmodels <- dim(muh)[2]
	p <- dim(muh)[1]
	ll <- vector(mode = "numeric", length = nmodels)
    
	for(kk in 1:nmodels){
        sigmah <- solve(thetah[, , kk])
        if(length(row_obs) > 0){
            Z <- scale(oXsub, center = muh[, kk], scale = FALSE)
            S <- crossprod(Z) / n_obs
            ll[kk] <- 0.5 * n_obs * (determinant(thetah[, , kk])$modulus - sum(S * thetah[, , kk]) - p * log(2 * pi))
        }
		for(i in 1:n_mis){
			ii <- row_mis[i]
			xo <- Xsub[ii, -col_mis[[i]]]
			muh_o <- muh[-col_mis[[i]], kk]
			muh_m <- muh[col_mis[[i]], kk]
			sigmah_oo <- sigmah[-col_mis[[i]], -col_mis[[i]], drop = FALSE]
			sigmah_mm <- sigmah[col_mis[[i]], col_mis[[i]], drop = FALSE]
			sigmah_mo <- sigmah[col_mis[[i]], -col_mis[[i]], drop = FALSE]
            ll[kk] <- ll[kk] + dmvnorm(xo, mean = muh_o, sigma = sigmah_oo, log = TRUE)
            M <- sigmah_mo %*% solve(sigmah_oo)
            muh_m_given_o <- drop(muh_m + M %*% (xo - muh_o))
            sigmah_m_given_o <- sigmah_mm - M %*% t(sigmah_mo)
            ptail <- pmvnorm(lower = k[col_mis[[i]]], mean = muh_m_given_o, sigma = sigmah_m_given_o)
            if(ptail > 0) ll[kk] <- ll[kk] + log(ptail)
		}
	}
    
	out <- list(loglik = ll, rho_seq = object$rho_seq)
	out
	
}













