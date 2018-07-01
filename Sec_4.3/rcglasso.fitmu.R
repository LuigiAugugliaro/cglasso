rcglasso.fitmu <- function(T1_obs, nobs, X, k, n_mis, row_mis, col_mis, rho, muh, thetah, nstp, tol, eps, verbose){
    conv <- 0
    muh_old <- muh
    if(verbose) cat(c("\t\ti", "dmu", "\n"))
    for(i in 1:nstp){
        T1_bar <- T1_obs
        for(h in 1:n_mis){
            ri <- col_mis[[h]]
            x_obs <- X[row_mis[h], -ri]
            muh_obs <- muh_old[-ri]
            muh_mis <- muh_old[ri]
            thetah_obs <- thetah[-ri, -ri, drop = FALSE]
            thetah_mis <- thetah[ri, ri, drop = FALSE]
            thetah_misobs <- thetah[ri, -ri, drop = FALSE]
            zeta <- drop(thetah_misobs %*% (x_obs - muh_obs))
            sigmah_mis_given_obs <- try(solve(thetah_mis), silent = TRUE)
            if(class(sigmah_mis_given_obs) == "try-error"){
                warning("The conditional covariance matrix can not be computed at rho = ", rho)
                out <- list(conv = 2)
                return(out)
            }
            muh_mis_given_obs <- muh_mis - drop(sigmah_mis_given_obs %*% zeta)
            tmean <- vector(mode  = "numeric", length = length(ri))
            for(m in 1:length(ri)){
                zi <- (k[ri][m] - muh_mis_given_obs[m]) / sqrt(sigmah_mis_given_obs[m, m])
                ptail <- pnorm(zi, lower.tail = FALSE)
                if(ptail > tol){
                    ratio <- dnorm(zi) / ptail
                    tmean[m] <- muh_mis_given_obs[m] + sqrt(sigmah_mis_given_obs[m, m]) * ratio
                } else tmean[m] <- muh_mis_given_obs[m]
            }
            T1_bar[ri] <- T1_bar[ri] + tmean
        }
        muh_new <- T1_bar / nobs
        dmuh <- max(abs(muh_new - muh_old))
        if(verbose) cat(c("\t\t", i, dmuh), "\n")
        if(dmuh < eps) break
        else muh_old <- muh_new
    }
    if(i == nstp){
        warning("maximum number of iterations has been reached (increase argument ", sQuote("nstp"), ")\n")
        out <- list(conv = 1)
        return(out)
    }
    out <- list(muh_new = muh_new, niter = i, conv = conv)
    out
}
