invisible(lapply(dir("R", full.names = TRUE), source))
library(rgl)
library(compositions)

Dg <- function(x) {
  if (is.vector(x)) x <- t(x)
  n <- nrow(x)
  D <- ncol(x)
  lapply(1:n, function(i) {
    diag(D) - matrix(x[i, ], D, D)
  })
}

pddirbase <- function(x, mix) {
  if (is.vector(x)) x <- t(x)
  n <- nrow(x)
  m <- length(mix$pr)
  logpr <- rep(log(mix$pr), each = n)
  logd <- ddir(x, mix$alpha, log = TRUE, outer = TRUE) + logpr
  rx <- x[rep(1:n, each = m), , drop = FALSE]
  ra <- mix$alpha[rep(1:m, n), , drop = FALSE]
  logg <- log(ra - 1) - log(rx)
  r <- logg + c(t(logd))
  rowsum(exp(r), rep(1:n, each = m))
}

gdvend <- function(x, pd, scale) {
  if (is.vector(x)) x <- t(x)
  if (is.vector(pd)) pd <- t(pd)
  n <- nrow(x)
  D <- ncol(x)
  xend <- x + pd * scale
}

simplexpanel3D <- function(alpha = 0.2, view = NULL) {
  if (is.null(view))
    view <- rbind(c(0.013, -0.702, 0.712, 0),
                  c(0.926, -0.260, -0.274, 0),
                  c(0.377, 0.663, 0.647, 0),
                  c(0, 0, 0 , 1))
  V <- diag(3)
  open3d()
  triangles3d(V, alpha = alpha)
  axes3d()
  view3d(userMatrix = view)
}

simplex3D_grid <- function(x = c(0, 1), margin = 1e-3, axl = 100) {
  ax <- seq(x[1] + margin, x[2] - margin, l = axl)
  x <- as.matrix(expand.grid(ax, ax))
  x <- cbind(x, 1 - rowSums(x))
  x <- x / rowSums(x)
  i <- rowSums(x >= margin) == 3
  list(grid = x, valid = i, ax = ax)
}

# 2D Plot
mix <- Dir_mix(4)
plot(mix, prop = F)
vfplot(mix, add = T)
dirmode(mix, plot = 2, add = T)

n <- 500
h <- 0.009758621
# simumix_1: 0.007689655, 0.03775862
# Dir_mix_4: 0.009758621, 0.04586207
c <- 1
data <- rdirmix(n, mix)
kn <- dirkn(data, h, c = c)
plot(kn, prop = F)
vfplot(kn, add = T)
dirmode(kn, st = attr(mix, "mode"), plot = 2, add = T)

vfplot(mix); vfplot(kn, col = 2, add = T); vfdiffplot(mix, kn, add = F)

# 3D Plot
gd <- simplex3D_grid(axl = 30)
x <- gd$grid[gd$valid, ]
v <- pddir(x, mix, proj = F)
ve <- pddir(x, kn, proj = F)

y <- gdvend(x, v, scale = 1e-3)
ye <- gdvend(x, ve, scale = 1e-3)
simplexpanel3D()
arrows3D(x, y, angle = 15, length = 0.15)
arrows3D(x, ye, angle = 15, length = 0.15, col = 2)

v1 <- v - rowMeans(v)
ve1 <- ve - rowMeans(ve)
y1 <- gdvend(x, v1, scale = 1e-3)
ye1 <- gdvend(x, ve1, scale = 1e-3)
simplexpanel3D()
arrows3D(x, y1, angle = 15, length = 0.15)
arrows3D(x, ye1, angle = 15, length = 0.15, col = 2)

Dgx <- Dg(x)
v2 <- t(mapply(function(i) v[i, ] %*% Dgx[[i]], 1:nrow(v)))
ve2 <- t(mapply(function(i) ve[i, ] %*% Dgx[[i]], 1:nrow(v)))
y2 <- gdvend(x, v2, scale = 1e-3)
ye2 <- gdvend(x, ve2, scale = 1e-3)
simplexpanel3D()
arrows3D(x, y2, angle = 15, length = 0.15)
arrows3D(x, ye2, angle = 15, length = 0.15, col = 2)

v3 <- v2 - rowMeans(v2)
ve3 <- ve2 - rowMeans(ve2)
y3 <- gdvend(x, v3, scale = 1e-3)
ye3 <- gdvend(x, ve3, scale = 1e-3)
simplexpanel3D()
arrows3D(x, y3, angle = 15, length = 0.15)
arrows3D(x, ye3, angle = 15, length = 0.15, col = 2)

all(round(v1, 10) == round(v3, 10))
