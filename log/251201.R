invisible(lapply(dir("R", full.names = TRUE), source))

# First test: one comp, one x.
x <- c(0.2, 0.8, 0)
mix <- dirmix(a = c(3, 4, 1))
plot(mix, prop = FALSE)
points_simplex(x, col = "grey30", add = TRUE)
(pd <- pddir(x, mix))
y <- gdvend(x, pd, scale = 1e-3)
arrows_simplex(x, y, add = TRUE, col = "blue")
dirmode(mix, st = x, plot = 2, col = "grey50", add = TRUE)

# Second test: one comp, grid x.
mix <- Dir_mix(3)
vfplot(mix, col = "blue", axl = 20)

x <- c(0.06, 0.02, 0.92)
(pd <- pddir(x, mix))
y <- gdvend(x, pd, scale = 1e-3)
arrows_simplex(x, y, add = TRUE, col = "red")
dirmode(mix, st = x, plot = 2, col = "grey50")
