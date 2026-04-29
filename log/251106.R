invisible(lapply(dir("R", full.names = TRUE), source))

# x <- rep(1, 3) / 3
x <- c(0.2, 0.3, 0.5)
# x <- c(0, 0.5, 0.5)
vh <- c(1, 0.5, 0.1, 0.05, 0.01)
vh <- seq(1, 0.2, by = -0.2)
oldpar <- par(no.readonly = TRUE)
par(mfcol = c(3, 5))
for (h in vh) {
  ait <- kn_ait(x, h)
  cc1 <- kn_cc1(x, h)
  cc2 <- kn_cc2(x, h)
  plot(ait)
  plot(cc1)
  plot(cc2)
}
par(oldpar)

mx <- rbind(rep(1, 3) / 3, c(0.3, 0.2, 0.5), c(0.4, 0.1, 0.5),
            c(0.45, 0.01, 0.54), c(0.45, 0, 0.55))
h <- 0.2
par(mfcol = c(3, 5))
for (i in seq_len(nrow(mx))) {
  ait <- kn_ait(mx[i, ], h)
  cc1 <- kn_cc1(mx[i, ], h)
  cc2 <- kn_cc2(mx[i, ], h)
  plot(ait)
  plot(cc1)
  plot(cc2)
}
par(oldpar)

h <- 0.2
x <- c(0.3, 0.1, 0.6)
(ait <- kn_ait(x, h))
(cc1 <- kn_cc1(x, h))
(cc2 <- kn_cc2(x, h))
data <- c(0.2, 0.7, 0.1)
dkn(data, ait)
dkn(data, cc1)

# Test for perturbation
xg <- simplex_grid(axl = 10)
x <- xg$grid[xg$valid, ]
points_simplex(x)
# X <- c(0, 0, 1)
X <- c(1, 0 ,0)
nx <- perturbe(x, 1/X)
points_simplex(nx, add = T, col = "blue")
points_simplex(X, col = "black", pch = 19, add = T, cex = 1)

# MISE
mix0 <- dirmix(rep(1, 3)/3 / 0.05 + 1) # 1
mix0 <- dirmix(
  rbind(c(0.1, 0.35, 0.55)/0.12,
        c(0.78, 0.13, 0.09)/0.34,
        c(0.34, 0.14, 0.52)/0.06) + 1,
  pr = c(0.24, 0.47, 0.29)) # 2
mix0 <- dirmix(
  rbind(c(0, 0.55, 0.45)/0.3, c(0, 0, 1)/0.1) + 1,
  pr = c(0.8, 0.2)) # 3
mix0 <- dirmix(
  rbind(c(0.12, 0.22, 0.43, 0.23)/0.24, c(0, 0.41, 0.04, 0.55)/0.08) + 1,
  pr = c(0.68, 0.32)) # 4
mix0 <- dirmix(
  rbind(c(0.87, 0.13)/0.24, c(0.24, 0.76)/0.08) + 1,
  pr = c(0.3, 0.7)) # 4

# plot(mix0)
h <- 0.1
n <- 5
B <- 1000

# 1.
mix1 <- samp_comp(mix0, mn = 1e5)
(mise1 <- (MISEcf(h, mix1, n)))

# 2.
mix2 <- samp_mix(mix0, n = n, B = B)
(mise2 <- MISEmc(h, mix2))

# 3.
mix3 <- grid_mix(mix0, n = 500)
(mise3 <- ISEdt(h, mix3))

# h vs MISE
(vh <- seq(0.05, 0.005, by = -0.002))
vh <- 0.05
B <- 10000
# 1.
vmise1 <- sapply(vh, function(h) {
  cat(h, " ")
  MISEcf(h, mix1, n)
})
# 2.
vmise2 <- sapply(vh, function(h) {
  cat(h, " ")
  MISEmc(h, mix2)
})
# 3.
vmise3 <- sapply(vh, function(h) {
  cat(h, " ")
  ISEdt(h, mix3)
})

pdf("MISE.pdf")
plot.new()
plot.window(xlim = rev(range(vh)), ylim = range(vmise1) + c(-0.004, 0.001))
box(); axis(1); axis(2)
lines(vh, vmise1)
lines(vh, vmise2, col = "red")
lines(vh, vmise3, col = "blue")
hat1 <- vh[which.min(vmise1)]
hat2 <- vh[which.min(vmise2)]
hat3 <- vh[which.min(vmise3)]
abline(v = hat1, lty = 2)
abline(v = hat2, lty = 2, col = "red")
# abline(v = hat3, lty = 2, col = "blue")
text(labels = hat1, x = hat1 + 0.002, y = mean(range(vmise1)))
text(labels = hat2, x = hat2 - 0.002, y = mean(range(vmise1)), col = "red")
text(labels = hat3, x = hat3 + 0.002, y = mean(range(vmise1)), col = "blue")
dev.off()

# LSCV
data <- rdirmix(n, mix0)
# h vs LSCV
lscv <- sapply(vh, function(h) LSCV(h, data))
plot(vh, lscv, type = 'l', xlim = rev(range(vh)))
# Compile for (11) -> first term only need upper triangle, check with MC LSCV
