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
	p <- length(muh)
    sigmah <- solve(thetah)
    ll <- 0
    
    if(length(row_obs) > 0){
        Z <- scale(oXsub, center = muh, scale = FALSE)
        S <- crossprod(Z) / n_obs
        ll <- 0.5 * n_obs * (determinant(thetah)$modulus - sum(S * thetah) - p * log(2 * pi))
    }
        
    for(i in 1:n_mis){
        ii <- row_mis[i]
        xo <- Xsub[ii, -col_mis[[i]]]
        muh_o <- muh[-col_mis[[i]]]
        muh_m <- muh[col_mis[[i]]]
        sigmah_oo <- sigmah[-col_mis[[i]], -col_mis[[i]], drop = FALSE]
        sigmah_mm <- sigmah[col_mis[[i]], col_mis[[i]], drop = FALSE]
        sigmah_mo <- sigmah[col_mis[[i]], -col_mis[[i]], drop = FALSE]
            
        ll <- ll + dmvnorm(xo, mean = muh_o, sigma = sigmah_oo, log = TRUE)
        M <- sigmah_mo %*% solve(sigmah_oo)
        muh_m_given_o <- drop(muh_m + M %*% (xo - muh_o))
        sigmah_m_given_o <- sigmah_mm - M %*% t(sigmah_mo)
        ptail <- pmvnorm(lower = k[col_mis[[i]]], mean = muh_m_given_o, sigma = sigmah_m_given_o)
        if(ptail > 0) ll <- ll + log(ptail)
    }
    
	out <- list(loglik = ll)
	out
	
}













