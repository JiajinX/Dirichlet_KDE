# Select optimal bandwidth for KDE by finite mixture

invisible(lapply(dir("R", full.names = TRUE), source))

mix0 <- Dir_mix(4)
n <- 500
data <- rdirmix(n, mix0)
est <- optfm(data, verbose = 1)
plot(mix0, prop = F)
plot(est$mix, prop = F)

scale <- 1e-3
vfplot(mix0, scale = scale)
vfplot(est$mix, scale = scale, col = 2, add = T)

mix0 <- gene_compsamp(mix0)

vh <- seq(0.01, 0.005, len = 10)
estmix <- gene_compsamp(est$mix)

h0 <- optMISEh(vh = vh, mix0, n, verbose = TRUE)
kn0 <- dirkn(x = data, h = h0$h)
plot(kn0, prop = F)

hest <- optMISEh(vh = vh, estmix, n, verbose = TRUE)
knest <- dirkn(x = data, h = hest$h)
plot(knest, prop = F)

# =======================================================

c <- 1.5
dvh <- seq(0.03, 0.005, len = 10)
dh0 <- optDMISEh(vh = dvh, mix0, n, c, verbose = TRUE)
dkn0 <- dirkn(x = data, h = dh0$h)
vfplot(dkn0, scale = scale)

dhest <- optDMISEh(vh = dvh, estmix, n, c, verbose = TRUE)
dknest <- dirkn(x = data, h = dhest$h)
vfplot(dknest, scale = scale)
