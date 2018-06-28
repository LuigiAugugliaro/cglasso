AIC.cglasso <- function(object, mof = c("Q", "loglik"), k = c("bic", "aic")){
    mof <- match.arg(mof)
    nobs <- object$nobs
    npar <- object$npar
    if(is.character(k)){
        k <- match.arg(k)
        k <- ifelse(k == "aic", 2, log(nobs))
    }
    if(k <= 0)
    stop("argument ", sQuote("k"), " can not be a negative value")
    fit <- if(mof == "loglik") loglik(object)$loglik
            else Q(object)$Q
    out <- -2 * fit + k * npar
    out <- list(gof = out, rho = object$rho_seq)
    out
}
