Q <- function(object){
    nobs <- object$nobs
    S <- object$S
    p <- dim(S)[1]
    thetah <- object$thetah
    Qval <- 0.5 * nobs * (determinant(thetah)$modulus - sum(S * thetah) - p * log(2 * pi))
    out <- list(Qval = Qval)
    out
}
