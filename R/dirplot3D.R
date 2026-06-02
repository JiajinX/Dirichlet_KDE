simplexpanel3D <- function(view = NULL, zoom = 1, axes = TRUE, alpha = 0.2) {
  if (is.null(view))
    view <- rbind(c(0.013, -0.702, 0.712, 0),
                  c(0.926, -0.260, -0.274, 0),
                  c(0.377, 0.663, 0.647, 0),
                  c(0, 0, 0, 1))
  V <- diag(3)
  rgl::open3d()
  rgl::triangles3d(V, alpha = alpha)
  if (axes) rgl::axes3d()
  rgl::view3d(userMatrix = view, zoom = zoom)
}

simplex3D_grid <- function(x = c(0, 1), margin = 1e-3, axl = 100) {
  ax <- seq(x[1] + margin, x[2] - margin, l = axl)
  x <- as.matrix(expand.grid(ax, ax))
  x <- cbind(x, 1 - rowSums(x))
  x <- x / rowSums(x)
  i <- rowSums(x >= margin) == 3
  list(grid = x, valid = i, ax = ax)
}

vfplot3D <- function(obj, add = FALSE, col = "blue", scale = 1e-3, axl = 30,
                     view = NULL, zoom = 1, axes = TRUE) {
  gd <- simplex3D_grid(axl = axl)
  x <- gd$grid[gd$valid, ]
  pd <- pddir(x, obj, proj = FALSE)
  y <- gdvend(x, pd, scale = scale)
  if (!add) simplexpanel3D(view, zoom, axes)
  compositions::arrows3D(x, y, col = col, length = 0.1)
}

# Barycentric coordinate -> Cartesian coordinates
# data: simplex data with 4-part components
b2c <- function(data) {
  if (is.vector(data)) data <- t(data)
  p <- matrix(c(sqrt(6) / 4, 3 / 4 * sqrt(2), 0,
                0, 0, 0,
                sqrt(6) / 2, 0, 0,
                sqrt(6) / 4, 1 / 4 * sqrt(2), 1),
              ncol = 3, byrow = TRUE)
  data %*% p
}

# Cartesian coordinates -> Barycentric coordinate
c2b <- function(data) {
  if (is.vector(data)) data <- t(data)
  p <- matrix(c(sqrt(6) / 4, 3 / 4 * sqrt(2), 0,
                0, 0, 0,
                sqrt(6) / 2, 0, 0,
                sqrt(6) / 4, 1 / 4 * sqrt(2), 1),
              ncol = 3, byrow = TRUE)
  inv_p <- solve(t(cbind(1, p)))
  tcrossprod(cbind(1, data), inv_p)
}

# Plot regular tetrahedron
tetrahedron <- function(add = FALSE) {
  p <- matrix(c(sqrt(6) / 4, 3 / 4 * sqrt(2), 0,
                0, 0, 0,
                sqrt(6) / 2, 0, 0,
                sqrt(6) / 4, 1 / 4 * sqrt(2), 1),
              ncol = 3, byrow = TRUE)
  if (add) open3d()
  rgl::segments3d(c(p[1, ], p[2, ]))
  rgl::segments3d(c(p[1, ], p[3, ]))
  rgl::segments3d(c(p[1, ], p[4, ]))
  rgl::segments3d(c(p[2, ], p[3, ]))
  rgl::segments3d(c(p[2, ], p[4, ]))
  rgl::segments3d(c(p[3, ], p[4, ]))

  # Add lab
  gap <- 0.05
  cex <- 1.5
  rgl::text3d(p[1, 1], p[1, 2] + gap, p[1, 3], text = expression(x[1]), cex = cex)
  rgl::text3d(p[2, 1] - gap, p[2, 2], p[2, 3], text = expression(x[2]), cex = cex)
  rgl::text3d(p[3, 1] + gap, p[3, 2], p[3, 3], text = expression(x[3]), cex = cex)
  rgl::text3d(p[4, 1], p[4, 2], p[4, 3] + gap, text = expression(x[4]), cex = cex)
}

# Points in regular tetrahedron
# data: simplex data with 4-part components
points_tetrahedron <- function(data, add = FALSE, col = 2) {
  coord <- b2c(data)
  if (!add) tetrahedron()
  rgl::points3d(coord[, 1], coord[, 2], coord[, 3], col = col)
}

contour3D <- function(obj, data = NULL, n = 50, lv = NULL,
                      col = NULL, alo = 0.1, ahi = 0.5, add = FALSE) {
  x <- seq(0, sqrt(6) / 2, len = n)
  y <- seq(0, 3 / 4 * sqrt(2), len = n)
  z <- seq(0, 1, len = n)
  k <- length(lv)
  if (is.null(col)) col <- viridis::viridis(k)
  al <- seq(alo, ahi, len = k)

  plot3D <- function(x, y, z) {
    coord <- c2b(cbind(x, y, z))
    coord <- coord / rowSums(coord)
    d <- rep(0, length = nrow(coord))
    idx <- which(rowSums(coord > 0) == 4)
    d[idx] <- dspx(coord[idx, ], obj)
    d
  }
  f <- function(x, y, z) plot3D(x, y, z)
  misc3d::contour3d(f, lv, x, y, z, color = col, alpha = al, add = add)
  if (!add) tetrahedron()
  # view3d(theta = 0, phi = -3, zoom = 0.2)
  rgl::view3d(theta = 7, phi = 25, zoom = 0.3)
}

# Draw the Decision boundary of two estimation on simplex
decbondary3D <- function(est1, est2, prior, threshold = NULL,
                         x1 = c(0, sqrt(6)/2), x2 = c(0, 3/4 * sqrt(2)),
                         x3 = c(0, 1), axl = 50, level = 0, alpha = 0.5,
                         add = FALSE, ...) {
  x <- seq(x1[1], x1[2], len = axl)
  y <- seq(x2[1], x2[2], len = axl)
  z <- seq(x3[1], x3[2], len = axl)
  plot3D <- function(x, y, z) {
    coord <- c2b(cbind(x, y, z))
    coord <- coord / rowSums(coord)
    idx <- which(rowSums(coord > 0) == 4)
    d <- rep(-Inf, length = nrow(coord))
    d1 <- dspx(coord[idx, ], est1)
    d2 <- dspx(coord[idx, ], est2)
    if (!is.null(threshold)) d1[d1 < threshold] <- -Inf
    d[idx] <- d1 * prior[1] - d2 * prior[2]
    d
  }
  f <- function(x, y, z) plot3D(x, y, z)
  misc3d::contour3d(f, level, x, y, z, alpha = alpha, add = add, ...)
  if (!add) tetrahedron()
  # view3d(theta = 0, phi = -3, zoom = 0.2)
  rgl::view3d(theta = 7, phi = 25, zoom = 0.3)
}
