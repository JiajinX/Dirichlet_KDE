# Create a simplex panel
simplexpanel <- function(lab = list(expression(x[1]),
                                    expression(x[2]),
                                    expression(x[3]))) {
  par(mar = rep(1, 4))
  liml <- 0.1
  labl <- 0.06
  sidel <- 2 / sqrt(3)
  tcex <- 0.8
  margin <- 0.05
  plot.new()
  plot.window(xlim = c(-liml, sidel + liml), ylim = c(-liml, 1 + liml),
              xaxs = "i", yaxs = "i", asp = 1)

  l <- matrix(c(0, 0, sidel, 0,
                0, 0, sidel / 2, 1,
                sidel / 2, 1, sidel, 0), nrow = 4)
  gap <- c(10, 25, 50, 75, 90) / 100
  temp <- gap / sqrt(3)
  text(sidel / 2, 1 + labl, lab[[1]])
  text(-labl / sidel, -labl / 2, lab[[2]], adj = 1)
  text(sidel + labl / sidel, -labl / 2, lab[[3]], adj = 0)
  segments(temp, gap, sidel - temp, gap, lty = 3, col = "grey")
  segments(temp, gap, 2 * temp, 0, lty = 3, col = "grey")
  segments(sidel - temp, gap, sidel - 2 * temp, 0, lty = 3, col = "grey")
  text(temp - margin, gap, gap * 100, cex = tcex)
  text(sidel - temp + margin, gap, (1 - gap) * 100, cex = tcex)
  text(2 * temp, -0.05, (1 - gap) * 100, cex = tcex)
  # Draw equilateral triangle, labs, axis value and auxiliary line
  segments(l[1, ], l[2, ], l[3, ], l[4, ], col = "grey50")
}

points_simplex <- function(x, add = FALSE, col = 2, pch = 20, cex = 0.8, ...) {
  if (!add) simplexpanel()
  if (is.vector(x)) x <- t(x)
  x <- cbind((2 - 2 * x[, 2] - x[, 1]) / sqrt(3), x[, 1])
  points(x, col = col, pch = pch, cex = cex, ...)
}

arrows_simplex <- function(x, y, add = FALSE, tol = 1e-6,
                           col = "grey60", length = 0.04, angle = 15, ...) {
  if (!add) simplexpanel()
  if (is.vector(x)) x <- t(x)
  if (is.vector(y)) y <- t(y)
  x <- cbind((2 - 2 * x[, 2] - x[, 1]) / sqrt(3), x[, 1])
  y <- cbind((2 - 2 * y[, 2] - y[, 1]) / sqrt(3), y[, 1])
  col <- rep(col, len = nrow(x))
  idx <- rowSums((y - x)^2) > tol
  arrows(x[idx, 1], x[idx, 2], y[idx, 1], y[idx, 2],
         col = col[idx], length = length, angle = angle, ...)
  points(x[!idx, ], col = col[idx], pch = 20, cex = 0.01)
}

simplex_grid <- function(x1 = c(0, 1), x2 = c(0, 2 / sqrt(3)), axl = 100,
                         boundary = TRUE) {
  ax1 <- seq(x1[1], x1[2], l = axl)
  ax2 <- seq(x2[1], x2[2], l = axl)
  x <- as.matrix(expand.grid(ax2, ax1))
  x <- cbind(x[, 2], 1 - (sqrt(3) * x[, 1] + x[, 2]) / 2)
  x <- cbind(x, 1 - rowSums(x))
  x[x < 1e-15 & x > 0] <- 0
  x <- x / rowSums(x)
  px <- if (boundary) x >= 0 else x > 0
  i <- rowSums(px) == 3
  list(grid = x, valid = i, ax1 = ax1, ax2 = ax2)
}

plot.dirkn <- function(kn, data = NULL, col = NULL, add = FALSE,
                       prop = TRUE, lv = NULL, plv = seq(0.5, 0.95, len = 4),
                       x1 = c(0, 1), x2 = c(0, 2 / sqrt(3)), axl = 100) {
  if (is.data.frame(data)) data <- as.matrix(data)
  x <- simplex_grid(x1, x2, axl) # Grid for calculating density
  d <- rep(NA, length = axl^2)
  d[x$valid] <- ddirkn(x$grid[x$valid, ], kn)
  if (prop) {
    fhat <- dirmix(kn$alpha)
    if (is.null(data)) data <- rdirmix(1000, fhat) # Generate 1000 data randomly
    ddata <- ddirkn(data, kn)
    d[x$valid] <- colMeans(outer(ddata, d[x$valid], ">="))
  }
  if (is.null(lv)) { # Density levels
    if (prop) {
      lv <- seq(30, 90, by = 20) / 100
    } else lv <- quantile(d[x$valid], probs = plv)
  }
  lb <- if (prop) paste0(lv * 100, "%") else round(lv, 2) # Labels of levels
  dim(d) <- c(axl, axl)
  if (is.null(col)) col <- paletteer::paletteer_c("viridis::viridis", length(lv))
  if (prop) col <- rev(col)
  if (!add) simplexpanel()
  contour(x$ax2, x$ax1, d, levels = lv, labels = lb,
          col = col, labcex = 0.8, add = TRUE)
}

plot.dirmix <- function(mix, data = NULL, col = NULL, add = FALSE, prop = TRUE,
                        lv = NULL, plv = seq(0.5, 0.95, len = 4),
                        x1 = c(0, 1), x2 = c(0, 2 / sqrt(3)), axl = 100) {
  if (is.data.frame(data)) data <- as.matrix(data)
  x <- simplex_grid(x1, x2, axl) # Grid for calculating density
  d <- rep(NA, length = axl^2)
  d[x$valid] <- ddirmix(x$grid[x$valid, ], mix)
  if (prop) {
    if (is.null(data)) data <- rdirmix(1000, mix) # Generate 1000 data randomly
    ddata <- ddirmix(data, mix)
    d[x$valid] <- colMeans(outer(ddata, d[x$valid], ">="))
  }
  if (is.null(lv)) { # Density levels
    if (prop) {
      lv <- seq(30, 90, by = 20) / 100
    } else lv <- quantile(d[x$valid], probs = plv)
  }
  lb <- if (prop) paste0(lv * 100, "%") else round(lv, 2) # Labels of levels
  dim(d) <- c(axl, axl)
  if (is.null(col)) col <- paletteer::paletteer_c("viridis::viridis", length(lv))
  if (prop) col <- rev(col)
  if (!add) simplexpanel()
  contour(x$ax2, x$ax1, d, levels = lv, labels = lb,
          col = col, labcex = 0.8, add = TRUE)
  points_simplex(attr(mix, "mode"), add = TRUE)
}

# Find the end points of the x based on gradient vector
gdvend <- function(x, pd, scale, unit = FALSE) {
  if (is.vector(x)) x <- t(x)
  if (is.vector(pd)) pd <- t(pd)
  vlen <- sqrt(rowSums(pd^2))
  # color
  ngp <- 25
  breaks <- quantile(vlen, seq(0, 1, length.out = ngp + 1), na.rm = TRUE)
  breaks <- unique(breaks)
  gp <- cut(vlen, breaks, labels = FALSE)
  ngp <- length(unique(gp))
  col <- paletteer::paletteer_c("grDevices::Grays", ngp, direction = -1) # grDevices::Oslo
  col <- col[gp]
  if (unit) pd <- pd / vlen
  xend <- x + pd * scale
  attr(xend, "col") <- col
  xend
}

# Vector field of the partial derivatives of Dirichlet mixture/kernel
# 1. Unit vector with magnitude by color?   2. Keep arrow inside the boundary?
vfplot <- function(obj, scale = NULL, unit = TRUE, add = FALSE,
                   x1 = c(0, 1), x2 = c(0, 2 / sqrt(3)), axl = 30) {
  if (is.null(scale)) scale <- ifelse(unit, 3e-2, 1e-3)
  x <- simplex_grid(x1, x2, axl) # Grid for calculating density
  xval <- x$grid[x$valid, ]
  pd <- pdspx(xval, obj)
  y <- gdvend(xval, pd, scale, unit)
  col <- attr(y, "col")
  if (unit) {
    idx <- !is.nan(y[, 1])
    col <- col[idx]
    xval <- xval[idx, , drop = FALSE]
    y <- y[idx, , drop = FALSE]
  }
  if (!add) simplexpanel()
  arrows_simplex(xval, y, col = col, add = TRUE)
}

vfdiffplot <- function(tf, est, scale = 1e-3, add = FALSE,
                       x1 = c(0, 1), x2 = c(0, 2 / sqrt(3)), axl = 30) {
  x <- simplex_grid(x1, x2, axl)
  xval <- x$grid[x$valid, ]
  pdt <- pddir(xval, tf)
  pde <- pddir(xval, est)
  pdif <- pdt - pde

  cols <- paletteer::paletteer_c("grDevices::Viridis", n = 30)
  lvec <- sqrt(rowSums(pdif^2))
  mapcol <- cols[as.numeric(cut(lvec, 30))]

  if (!add) simplexpanel()
  yt <- gdvend(xval, pdt, scale)
  ye <- gdvend(xval, pde, scale)
  arrows_simplex(yt, ye, col = mapcol, add = TRUE)
  # y <- gdvend(xval, pdif, scale)
  # arrows_simplex(xval, y, col = mapcol, add = TRUE)
}

dirmode <- function(obj, st = NULL, steplen = 1e-3, tol = 1e-10, maxit = 1e3,
                    add = FALSE, plot = 0, col = "blue") {
  if (is.null(st)) {
    # If no start points provided, then the default start points are the modes
    # of each mixture components
    if (class(obj) == "dirmix") st <- attr(mix, "mode")
    else stop("No initial points set")
  }
  if (is.vector(st)) st <- t(st)
  if (!add && plot != 0) simplexpanel()
  if (plot > 1) points_simplex(st, add = TRUE, col = col, pch = 1)
  for (i in seq_len(maxit)) {
    Df <- pddir(st, obj)
    ed <- stepup(st, steplen * Df)
    if (plot == 3)
      arrows_simplex(st, ed, col = col, add = TRUE, length = 0, tol = tol)
    if (max(ed - st) < tol) break
    st <- ed
  }
  if (plot > 0) points_simplex(ed, add = TRUE, col = col)
  ed[!duplicated(round(ed, 9)), , drop = FALSE]
}

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
