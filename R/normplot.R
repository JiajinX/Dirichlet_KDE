contour_spx <- function(mix, data = NULL, col = NULL, add = FALSE, prop = TRUE,
                        lv = NULL, plv = seq(0.5, 0.95, len = 4),
                        x1 = c(0, 1), x2 = c(0, 2 / sqrt(3)), axl = 100) {
  x <- simplex_grid(x1, x2, axl, boundary = FALSE) # Grid for calculating density
  d <- rep(NA, length = axl^2)
  xts <- ilrts(x$grid[x$valid, ])
  d[x$valid] <- dSnormix(xts, mix)
  if (prop) {
    if (is.null(data)) data <- rnormix(1000, mix) # Generate 1000 data randomly
    ddata <- dSnormix(data, mix)
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
  points_simplex(ilrts(mix$mu)$sx, add = TRUE)
}

contour_euc <- function(mix, data = NULL, col = NULL, add = FALSE, prop = TRUE,
                        lv = NULL, plv = seq(0.5, 0.95, len = 4),
                        x1 = c(-3, 3), x2 = c(-3, 3), axl = 100) {
  if (!add) {
    plot.new()
    plot.window(xlim = x1, ylim = x2)
    axis(1)
    axis(2)
    box()
  }
  x1 <- seq(x1[1], x1[2], l = axl)
  x2 <- seq(x2[1], x2[2], l = axl)
  x <- as.matrix(expand.grid(x1, x2))
  d <- dnormix(x, mix)
  if (prop) {
    if (is.null(data)) data <- rnormix(1000, mix) # Generate 1000 data randomly
    ddata <- dnormix(data, mix)
    d <- colMeans(outer(ddata, d, ">="))
  }
  if (is.null(lv)) { # Density levels
    if (prop) {
      lv <- seq(30, 90, by = 20) / 100
    } else lv <- quantile(d, probs = plv)
  }
  lb <- if (prop) paste0(lv * 100, "%") else round(lv, 2) # Labels of levels
  dim(d) <- c(axl, axl)
  if (is.null(col)) col <- paletteer::paletteer_c("viridis::viridis", length(lv))
  if (prop) col <- rev(col)
  contour(x1, x2, d, levels = lv, labels = lb,
          col = col, labcex = 0.8, add = TRUE)
  points(mix$mu, col = 2, pch = 20, cex = 0.8)
}

plot.normix <- function(mix, euc = FALSE, data = NULL,
                        col = NULL, add = FALSE, prop = TRUE,
                        lv = NULL, plv = seq(0.5, 0.95, len = 4),
                        x1 = NULL, x2 = NULL, axl = 100) {
  if (euc) {
    if (is.null(x1)) x1 = c(-3, 3)
    if (is.null(x2)) x2 = c(-3, 3)
    fun <- contour_euc
  } else {
    if (is.null(x1)) x1 = c(0, 1)
    if (is.null(x2)) x2 = c(0, 2 / sqrt(3))
    fun <- contour_spx
  }
  fun(mix, data, col, add, prop, lv, plv, x1, x2, axl)
}
