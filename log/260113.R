# LSCV for derivative density when r = 1

invisible(lapply(dir("R", full.names = TRUE), source))

c <- 1.5
h <- 0.1
n <- m <- 100
mix <- dirmix(a = c(8, 5, 15)) # f
x <- rdirmix(n, mix)
kn <- dirkn(x, h = h, c) # f_hat

# By therom
d <- ddir(x, kn$alpha, outer = T)
diag(d) <- 0
rx <- x[rep(1:n, each = m), , drop = FALSE]
ra <- kn$alpha[rep(1:m, n), , drop = FALSE]
gm <- exp(log(ra - 1) - log(rx))
g2m <- t(apply(gm, 1, function(g) kronecker(g, g)))
Dg <- exp(log(ra - 1) - 2 * log(rx))
Dg <- -t(apply(Dg, 1, function(r) c(diag(r))))
D2 <- g2m + Dg # g2m + Dg == t(D2) (in previous function)
r <- c(t(d)) * D2

# Projected
# Eq1
D <- ncol(x)
P <- diag(D) - 1/D
-c(P) %*% colSums(r)/(n * (n - 1))
# Eq2
P2 <- kronecker(P, P)
Pr <- r %*% t(P2)
-sum(c(diag(D)) * colSums(Pr) / (n * (n - 1)))

# MC method
nsp <- 1e4
unidir <- dirmix(a = rep(1,3))
samp <- rdirmix(nsp, unidir)
proj = T

# term2
Dkn <- pddir(samp, obj = kn, proj = proj)
Dfm <- pddir(samp, obj = mix, proj = proj)
mean(rowSums(Dkn * Dfm)) / 2

# term1
t1 <- pddir(samp, obj = kn, proj = proj)
mean(rowSums(t1 * t1))/2

# ========================================================
mix <- Dir_mix(8)
x <- rdirmix(100, mix)
c <- 1.5
vh <- seq(0.5, 0.05, len = 10)
(rth <- sapply(vh, function(h) DLSCV(h, x, c)))
(rmc <- sapply(vh, function(h) DLSCVmc(h, x, c, mix)))

plot.new()
plot.window(xlim = rev(range(vh)), ylim = range(rth, rmc))
axis(1); axis(2); box()
lines(vh, rth, type = "o")
lines(vh, rmc, type = "o", col = 2)
