load("CompAsp.RData")

##################################
# Table 1 Part A: Average CPU time
##################################

dimnames(time_table)
m_time <- apply(time_table, 3, function(M) apply(M, 2, mean))
sd_time <- apply(time_table, 3, function(M) apply(M, 2, sd))
round(rbind(m_time[1, ], sd_time[1, ], m_time[2, ], sd_time[2, ]), 2)

#######################################################
# Table 1 Part B: Difference between the two estimators
#######################################################

tbl <- array(0, dim = c(length(d), 2, 2), dimnames = list(d = d, dist = c("mu", "theta"), sum = c("mean", "sd")))

for(h in 1:length(d)){
    tbl[h, "mu", "mean"]    <- mean(apply(muh_f[, , h], 2, max))
    tbl[h, "mu", "sd"]      <- sd(apply(muh_f[, , h], 2, max))
    tbl[h, "theta", "mean"] <- mean(apply(thetah_f[, , h], 2, max))
    tbl[h, "theta", "sd"]   <- sd(apply(thetah_f[, , h], 2, max))
}

ftable(tbl, col.vars = 1, row.vars = c(2, 3))



























