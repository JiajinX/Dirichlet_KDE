# If we reduce the covariance matrix of the Gaussian mixture,
# the optimal h of the first derivative will be larger than ones in density.

#1.
mix0 <- list(mu = rbind(c(0, 0)),
             Sigma = ks::invvech(c(1/4, 0, 1/2)),
             pr = 1)
# [r=0] h=0.0408
# [r=1] h=0.0527

#2.
mix0 <- list(mu = rbind(c(1, 0), c(-1, 0)) * 3/4,
            Sigma = rbind(diag(1/5, 2), diag(1/5, 2)),
            pr = c(1, 1) / 2)
# [r=0] h=0.0213
# [r=1] h=0.0247

#3.
mix0 <- list(mu = rbind(c(1, -1), c(-1, 1)) * 2/5,
            Sigma = rbind(ks::invvech(c(1/4, 3/20, 1/4)), diag(1/4, 2)),
            pr = c(1, 1) / 2)
class(mix0) <- "normix"
plot(mix0, euc = T)
plot(mix0)

#4....
