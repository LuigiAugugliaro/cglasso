impute_rcglasso <- function(object){
    X <- object$X
    k <- object$k
    R <- object$R
    idsub <- apply(R, 1, all)
    if(any(idsub)){
        idsub <- !idsub
        Xsub <- X[idsub, , drop = FALSE]
        R <- R[idsub, , drop = FALSE]
    } else Xsub <- X
    nobs <- dim(Xsub)[1]
    row_mis <- which(apply(R, 1, any))
    row_obs <- setdiff(1:nobs, row_mis)
    if(length(row_obs) > 0) Xsub_obs <- X[row_obs, , drop = FALSE]
    n_mis <- length(row_mis)
    col_mis <- apply(R[row_mis, , drop = FALSE], 1, which)
    
    nrho <- length(object$rho_seq)
    Xnew <- array(0, dim = c(nobs, dim(Xsub)[2], nrho))
    
    for(kk in 1:nrho){
        Xnew[, , kk] <- Xsub
        muh <- object$muh[, kk]
        Thetah <- object$thetah[, , kk]
        for(h in 1:n_mis){
            i <- row_mis[h]
            ri <- col_mis[[h]]
            x_o <- Xsub[i, -ri]
            muh_o <- muh[-ri]
            muh_m <- muh[ri]
            Thetah_oo <- Thetah[-ri, -ri, drop = FALSE]
            Thetah_mm <- Thetah[ri, ri, drop = FALSE]
            Thetah_mo <- Thetah[ri, -ri, drop = FALSE]
            zeta <- drop(Thetah_mo %*% (x_o - muh_o))
            Sigmah_m_given_o <- drop(solve(Thetah_mm))
            muh_m_given_o <- muh_m - drop(Sigmah_m_given_o %*% zeta)
            out <- mtmvnorm(mean = muh_m_given_o, sigma = Sigmah_m_given_o, lower = k[ri], doComputeVariance = FALSE)
            Xnew[i, ri, kk] <- out$tmean
        }
    }
    Xnew
}
