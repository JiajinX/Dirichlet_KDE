# LSCV for derivative density when r = 1, closed form version.

invisible(lapply(dir("R", full.names = TRUE), source))

mix1 <- dirmix(a = c(8, 5, 15))
mix2 <- dirmix(a = c(10, 14, 7))

D <- ncol(mix1$alpha)
beta <- mix1$alpha + mix2$alpha - 1

# Left-hand side
a10 <- sum(mix1$alpha)
a20 <- sum(mix2$alpha)
t <- a10 + a20 - D - 1
t1 <- t * (t - 1)

t2 <- exp(attr(mix1$alpha, "logmbeta") + attr(mix2$alpha, "logmbeta") - Cdir(beta, log = T))

a1_1 <- mix1$alpha - 1
a2_1 <- mix2$alpha - 1
a1a2_2 <- a1_1 + a2_1
subt1 <- sum(a1_1 * a2_1 / (a1a2_2 * (a1a2_2 - 1)))
subt2 <- sum(a1_1 * a2_1 / a1a2_2^2)
subt3 <- sum(a1_1 / a1a2_2) * sum(a2_1 / a1a2_2)
t3 <- (D - 1) / D * subt1 + subt2 / D - subt3 / D
t1 * t2 * t3

# Right-hand side
beta0 <- sum(beta)
K <- t2 * (beta0 - 1) * (beta0 - 2)
M <- outer(c(a1_1), c(a1_1)) / outer(c(beta - 1), c(beta - 1))
diag(M) <- (a1_1) * (a1_1 - 1) / ((beta - 1) * (beta - 2))
KM <- K * M
P <- diag(D) - 1/D
-c(P) %*% c(KM)

# Test for Dirichlet mixture ===================================================

n <- 500
h <- 0.1
c <- 1.5
mix <- dirmix(a = rbind(c(8, 12, 5), c(12, 7, 9), c(20, 5, 5)),
              pr = c(0.2, 0.1, 0.7))
data <- rdirmix(n, mix)
kn <- dirkn(data, h, c = c)

# Left-hand side
a0 <- rowSums(mix$alpha)
D <- ncol(mix$alpha)
outCay <- outer(Cdir(kn$alpha, log = T), attr(mix$alpha, "logmbeta"), "+")
Cay <- Cdir_plus(kn$alpha, mix$alpha, log = TRUE)
Cd <- exp(outCay - Cay)
Cs <- D1S(kn$alpha, mix$alpha)
logC <- prefactor(1/h + c * D, a0, D)
r <- Cd * Cs * exp(logC[rep(1, n), ])
sum(mix$pr * t(r)) / n

# Right-hand side
a1 <- kn$alpha
a2 <- mix$alpha
a1n <- nrow(a1)
a2n <- nrow(a2)
beta <- a1[rep(seq_len(a1n), a2n), ] + a2[rep(seq_len(a2n), each = a1n), ] - 1
beta0 <- rowSums(beta)
dim(beta0) <- c(a1n, a2n)
K <- Cd * (beta0 - 1) * (beta0 - 2)

a1out <- t(apply(a1 - 1, 1, function(x) as.vector(outer(x, x))))
betaout <- t(apply(beta - 1, 1, function(x) as.vector(outer(x, x))))
combM <- a1out[rep(seq_len(a1n), a2n), ] / betaout
# --
a1diag <- (a1 - 1) * (a1 - 2)
betadiag <- (beta - 1) * (beta - 2)
diagM <- a1diag[rep(seq_len(a1n), a2n), ] / betadiag
# --
combM[, seq(1, D*D, by = D + 1)] <- diagM

re <- colSums(c(K) * combM * rep(mix$pr, each = a1n) / a1n)
P <- diag(D) - 1/D
-c(P) %*% re

# By Monte Carlo
douter <- ddir(kn$x, kn$alpha, outer = TRUE)
diag(douter) <- 0
rx <- kn$x[rep(1:n, each = n), , drop = FALSE]
ra <- kn$alpha[rep(1:n, n), , drop = FALSE]
log_ra_1 <- log(ra - 1)
log_rx <- log(rx)
gm <- exp(log_ra_1 - log_rx)
g2m <- t(apply(gm, 1, function(g) kronecker(g, g)))
Dg <- exp(log_ra_1 - 2 * log_rx)
Dg <- -t(apply(Dg, 1, function(r) c(diag(r))))
D2 <- g2m + Dg # g2m + Dg == t(D2) (in previous function)
r2 <- c(t(douter)) * D2
P <- diag(D) - 1/D
-c(P) %*% colSums(r2) / ( n * (n - 1))
