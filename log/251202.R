invisible(lapply(dir("R", full.names = TRUE), source))

# MISE
n <- 500
mix <- Dir_mix(4)
mix <- gene_compsamp(mix)
mix <- gene_eval(mix, neva = 1e4)
# mix <- gene_grid(mix, split = 50)
mix <- gene_data(n, mix, B = 50)

vh <- seq(0.015, 0.008, len = 20)

mf <- MISEfh(vh, mix, n)
plot_mise(mf)
mmc <- MISEh(vh, mix, fun = MISEmc, verbose = T)
# mdt <- MISEh(vh, mix, fun = MISEdt)

plot.new()
plot.window(xlim = rev(range(vh)), ylim = range(mf, mmc$MISE))
box(); axis(1); axis(2)
plot_mise(mf, add = T)
plot_mise(mmc$MISE, add = T, col = 2)
# plot_mise(mdt$MISE, add = T, col = 3)
MISEtest(mmc$ISE, mf)
# MISEtest(mdt$ISE, mf)

plot(mix, prop = F)
dirmode(mix, add = T, plot = 1)
vfplot(mix, add = T)
hopt <- vh[which.min(mf)]
kn <- dirkn(rdirmix(n, mix), h = hopt)
plot(kn, prop = F)
dirmode(kn, st = attr(mix, "mode"), add = T, plot = 1)
vfplot(kn, add = T)
# dirmode(kn, st = rdirmix(20, mix), plot = 2, add = T)

scale <- 1e-3
vfplot(mix, scale = scale, col = 1)
vfplot(kn, add = T, scale = scale, col = 2)
vfdiffplot(mix, kn, add = T, scale = scale)



