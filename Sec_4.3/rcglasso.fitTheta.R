rcglasso.fitTheta<- function(R, nobs, X, k, n_mis, row_mis, col_mis, rho, muh, thetah, nstp, tol, eps, verbose){
    Z <- scale(X, center = muh, scale = FALSE)
    tau <- k - muh
    Z[R] <- 0
    T2_obs <- crossprod(Z)
    thetah_old <- thetah
    conv <- 0
    if(verbose) cat(c("\t\ti", "dthetah", "\n"))
    for(i in 1:nstp){
        T2_bar <- T2_obs
        for(h in 1:n_mis){
            ri <- col_mis[[h]]
            zobs <- Z[row_mis[h], -ri]
            thetah_obs <- thetah_old[-ri, -ri, drop = FALSE]
            thetah_mis <- thetah_old[ri, ri, drop = FALSE]
            thetah_misobs <- thetah_old[ri, -ri, drop = FALSE]
            zeta <- drop(thetah_misobs %*% zobs)
            sigmah_mis_given_obs <- try(solve(thetah_mis), silent = TRUE)
            if(class(sigmah_mis_given_obs) == "try-error"){
                warning("The conditional covariance matrix can not be computed at rho = ", rho)
                out <- list(conv = 2)
                return(out)
            }
            muh_mis_given_obs <- - drop(sigmah_mis_given_obs %*% zeta)
            tmean <- vector(mode  = "numeric", length = length(ri))
            tvar <- diag(length(ri))
            for(m in 1:length(ri)){
                zi <- (tau[ri][m] - muh_mis_given_obs[m]) / sqrt(sigmah_mis_given_obs[m, m])
                ptail <- pnorm(zi, lower.tail = FALSE)
                if(ptail > tol){
                    ratio <- dnorm(zi) / ptail
                    tmean[m] <- muh_mis_given_obs[m] + sqrt(sigmah_mis_given_obs[m, m]) * ratio
                    tvar[m, m] <- sigmah_mis_given_obs[m, m] * (1 + zi * ratio - ratio^2)
                } else {
                    tmean[m] <- muh_mis_given_obs[m]
                    tvar[m, m] <- sigmah_mis_given_obs[m, m]
                }
            }
            T2_bar[-ri, ri] <- T2_bar[-ri, ri, drop = FALSE] + outer(zobs, tmean)
            T2_bar[ri, -ri] <- t(T2_bar[-ri, ri, drop = FALSE])
            T2_bar[ri, ri] <- T2_bar[ri, ri] + tvar + outer(tmean, tmean)
        }
        if(i == nstp){
            warning("maximum number of iterations has been reached (increase argument ", sQuote("nstp"), ")\n")
            out <- list(conv = 1)
            return(out)
        }
        S_bar <- T2_bar / nobs
        out <- glasso(S_bar, rho = rho, penalize.diagonal = FALSE, wi.init = thetah_old)
        thetah_new <- (out$wi + t(out$wi)) / 2
        dthetah <- max(abs(thetah_new - thetah_old))
        if(verbose) cat(c("\t\t", i, dthetah), "\n")
        if(dthetah < eps) break
        else thetah_old <- thetah_new
    }
    out <- list(S_bar = S_bar, thetah = thetah_new, niter = i, conv = conv)
    out
}
