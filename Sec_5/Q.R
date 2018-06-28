Q <- function(object){
    nobs <- object$nobs
    S <- object$S
    p <- dim(S)[1]
    thetah <- object$thetah
    nmodels <- length(object$rho_seq)
    Qval <- vector(mode = "numeric", length = nmodels)
    for(k in 1:nmodels)
        Qval[k] <- 0.5 * nobs * (determinant(thetah[, , k])$modulus - sum(S[, , k] * thetah[, , k]) - p * log(2 * pi))
    out <- list(Qval = Qval, rho_seq = object$rho_seq)
    out
}
