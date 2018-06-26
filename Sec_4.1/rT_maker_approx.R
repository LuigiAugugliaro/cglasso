rT_maker_approx <- function(T1_obs, T2_obs, k, X, n_mis, row_mis, col_mis, muh, Thetah){
    # inizializzazione statistiche T1 e T2
    T1_bar <- T1_obs
    T2_bar <- T2_obs
    for(h in 1:n_mis){
        i <- row_mis[h]     # indice della riga con almeno un missing
        ri <- col_mis[[h]]  # indice delle colonne con missing nell'i-esima riga
        x_obs <- X[i, -ri]  # sub-vettore osservato nell'i-esima riga
        muh_obs <- muh[-ri]
        muh_mis <- muh[ri]
        Thetah_obs <- Thetah[-ri, -ri, drop = FALSE]
        Thetah_mis <- Thetah[ri, ri, drop = FALSE]
        Thetah_misobs <- Thetah[ri, -ri, drop = FALSE]
        zeta <- drop(Thetah_misobs %*% (x_obs - muh_obs))
        Sigmah_mis_given_obs <- solve(Thetah_mis)
        muh_mis_given_obs <- muh_mis - drop(Sigmah_mis_given_obs %*% zeta)
        tmean <- vector(mode  = "numeric", length = length(ri))
        tvar <- diag(length(ri))
        for(m in 1:length(ri)){
            out <- mtmvnorm(mean = muh_mis_given_obs[m], sigma = Sigmah_mis_given_obs[m, m], lower = k[ri][m])
            tmean[m] <- out$tmean
            tvar[m, m] <- out$tvar
        }
        #update delle statistiche di sintesi
        T1_bar[ri] <- T1_bar[ri] + tmean
        T2_bar[-ri, ri] <- T2_bar[-ri, ri, drop = FALSE] + outer(x_obs, tmean)
        T2_bar[ri, -ri] <- t(T2_bar[-ri, ri, drop = FALSE])
        T2_bar[ri, ri] <- T2_bar[ri, ri] + tvar + outer(tmean, tmean)
    }
    out <- list(T1_bar = T1_bar, T2_bar = T2_bar)
    out
}
