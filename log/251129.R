invisible(lapply(dir("R", full.names = TRUE), source))

mix <- Dir_mix(4)
plot(mix, prop = F)
mixmode(mix, add = TRUE, plot = 1)
n <- 500
hlist <- readRDS("data/251126/Dirmix4_df")
h <- hlist$h[hlist$n == n]
X <- rdirmix(n, mix)
c <- 1
knplot(X, h, c = c, prop = F)
knmode(X, h, c = 1.5, st = attr(mix, "mode"), add = TRUE, plot = 1)

scale <- 1e-3
vfmixplot(mix, scale = scale, col = 2)
vfknplot(X, h = h, c = c, col = "grey50", add = TRUE, scale = scale)
vfdiffplot(mix, X, h = h, c = c, add = T, scale = scale)

# Mode =========================================================================
mix <- Dir_mix(10)
# mix <- dirmix(a = rbind(c(0, 0.2, 0.8) / 0.05, c(0, 0.1, 0.9) / 0.06) + 1)
plot(mix, prop = FALSE)
vfmixplot(mix, add = F, scale = 1e-3)
(tmd <- mixmode(mix, plot = 2, add = T))
