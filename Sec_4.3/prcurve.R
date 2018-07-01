prcurve <- function(thetah, thetat, rho){
	U <- upper.tri(thetat, diag = FALSE)
	thetat_v <- thetat[U]
	A <- which(abs(thetat_v) > 0)
	thetah_m <- apply(thetah, 3, function(M) M[U])
    
    nrho <- length(rho)
    recall <- apply(abs(thetah_m[A, ]) > 0, 2, mean)
    precision <- vector(mode = "numeric", length = nrho)
    for(i in 1:nrho){
        id <- which(abs(thetah_m[, i]) > 0)
        precision[i] <- ifelse(length(id) == 0, NA, mean(abs(thetat_v[id]) > 0))
    }
    id <- order(recall)
    precision <- precision[id]
    recall <- recall[id]
    
    recall[nrho] <- 1
    min_prec <- length(A) / length(thetat_v)
    precision[nrho] <- min_prec

    dprecision <- c(diff(precision), 0)
    drecall <- c(diff(recall), 0)
    auc <- sum(na.omit(precision * drecall)) + sum(na.omit(dprecision * drecall)) / 2
    
    out <- list(recall = recall, precision = precision, auc = auc, rho = rho[id])
	out
}
