invisible(lapply(dir("R", full.names = TRUE), source))
mix <- simumix(9, c = 1.5)
mix <- gene_compsamp(mix, 5)
n <- 5
vh <- seq(1, 0.1, len = 10)
re <- sapply(vh, function(h) {
  cat(round(h, 4), " ")
  DMISEf(h, mix, n, c = 1.5)
})
plot(vh, re, xlim = rev(range(vh)), type = "o")

D <- 3
P <- diag(D) - matrix(1/D, D, D)

# =======================================
mix <- Dir_mix(1)
x <- rdirmix(1e5, dirmix(a = rep(1, 3)))
d <- ddirpred(x, mix)
mean(d^2) / gamma(3)
Psi0(mix)

# =======================================
# Check DPsi0
mix <- simumix(9, c = 1.5)
x <- rdirmix(1e4, dirmix(a = rep(1, 3)))
pd <- pddir(x, mix, proj = T)
mean(rowSums(pd^2)) / gamma(3)
t(mix$pr) %*% DPsi0(mix) %*% mix$pr

# =======================================
# Check DPsi1
proj = T
c <- 1.5
mix <- simumix(9, c = c)
x <- rdirmix(1e4, dirmix(a = rep(1, 3)))
pd <- pddir(x, mix, proj = proj)

h <- 0.1
data <- rdirmix(500, mix)
kn <- dirkn(data, h = h, c = c)
epd <- pddir(x, kn, proj = proj)
mean(rowSums(pd*epd)) / gamma(3)

mix <- gene_compsamp(mix)
samp <- attr(mix, "samp")
attr(h, "y") <- lapply(samp, function(xi) {
  ya <- xi / h + c
  Cya <- Cdir(ya, log = TRUE)
  list(a = ya, Ca = Cya)
})
t(mix$pr) %*% DPsi1(h, mix, c = c) %*% mix$pr

# =======================================
# Check DPsi2
proj = T
c <- 1.5
mix <- simumix(9, c = c)
x <- rdirmix(1e5, dirmix(a = rep(1, 3)))

h <- 0.1; n <- 10; B <- 50
rterm1 <- rep(NA, B)
for (i in 1:B) {
  print(i)
  data <- rdirmix(n, mix)
  kn <- dirkn(data, h = h, c = c)
  epd <- pddir(x, kn, proj = proj)
  rterm1[i] <- mean(rowSums(epd^2)) / gamma(3)
}
mean(rterm1)

mix <- gene_compsamp(mix)
samp <- attr(mix, "samp")
attr(h, "y") <- lapply(samp, function(xi) {
  ya <- xi / h + c
  Cya <- Cdir(ya, log = TRUE)
  list(a = ya, Ca = Cya)
})
(dpsi2 <- t(mix$pr) %*% ((1 - 1/n) * DPsi2(h, mix, c = c)) %*% mix$pr + DCdh(h, mix, c) / n)
t.test(rterm1, mu = c(dpsi2))

# =======================================
# Check MISE
n <- 10
h <- 1

mix <- simumix(9, c = 1.5)
mix <- gene_compsamp(mix)
mix <- gene_data(n = n, mix, B = 100)
mix <- gene_eval(mix, r = 1, neva = 1e5)

(mcr <- MISEmc(h, mix, r = 1, c = 1.5, verbose = T))
(thor <- DMISEf(h, mix, n = n, c = 1.5))
t.test(attr(mcr, "ISE"), mu = thor)

vh <- seq(0.45, 0.2, len = 20)
re <- DMISEfh(vh, mix, n, c = 1.5)
plot_mise(re)
