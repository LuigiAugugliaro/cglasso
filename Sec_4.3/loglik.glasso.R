loglik_func <- function(est, s, n, gamma = 0){

loglik = rep(0,length(est$rholist))
df <- vector()
for (i in 1:length(est$rholist)){
	loglik[i] <- determinant(est$wi[,,i],logarithm=T)$modulus[[1]]-
	  sum(as.matrix(est$wi[,,i])*s)
	df[i]  <- sum((est$wi[,,i] != 0)+0)
	}

est$loglik = loglik
est$bic = -n*est$loglik + log(n)*df + 4*gamma*log(nrow(s))*df

return(est)
}
