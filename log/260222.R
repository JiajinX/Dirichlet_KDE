# Limit the boundary for the finite mixture (1.5 for first derivative of MISE)
invisible(lapply(dir("R", full.names = TRUE), source))

# Create a Dirichlet mixture closed to the boundary
mix0 <- Dir_mix(1)
mix0 <- dirmix(a = c(1.6, 7, 12))
mix0 <- dirmix(a = rbind(c(2, 6, 8), c(4, 10, 1.6), c(8, 8, 3)),
               pr = c(0.3, 0.5, 0.2))
mix0 <- dirmix(a = rbind(c(1.6, 6, 8), c(4, 10, 1.6), c(8, 8, 3)),
               pr = c(0.3, 0.5, 0.2))
plot(mix0, prop = F)
# Generate data from mix0
n <- 1000
data <- rdirmix(n, mix0)


# For the density estimation #
## Find optimal bandwidth by MISE ##
(fn <- optfm(data))
plot(fn$mix, prop = F)
mix <- gene_compsamp(fn$mix)
vh = seq(0.01, 0.001, len = 10)
(miseh <- optMISEh(vh, mix, n, verbose = TRUE))
knmise <- dirkn(x = data, h = miseh$h)
plot(knmise, prop = F)
## Find optimal bandwidth by LSCV ##
(lscvh <- optLSCVh(vh, data, verbose = TRUE))
knlscv <- dirkn(x = data, h = lscvh$h)
plot(knlscv, prop = F)


# For the first derivative of density estimation #
## Find optimal bandwidth by DMISE ##
(fn <- optfm(data, bound = 1.5))
plot(fn$mix, prop = F)
mix <- gene_compsamp(fn$mix)
vh = seq(0.03, 0.005, len = 10)
(dmiseh <- optMISEh(vh, mix, n, r = 1, verbose = TRUE))
kndmise <- dirkn(x = data, h = dmiseh$h)
vfplot(mix0, scale = 2e-3)
vfplot(kndmise, scale = 2e-3, add = T, col = 2)
## Find optimal bandwidth by DLSCV ##
(dlscvh <- optLSCVh(vh, data, r = 1, verbose = TRUE))
kndlscv <- dirkn(x = data, h = dlscvh$h)
vfplot(mix0, scale = 2e-3)
vfplot(kndlscv, scale = 2e-3, add = T, col = 2)
