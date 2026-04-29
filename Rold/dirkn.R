perturbe <- function(x, y) {
  if (is.vector(x)) x <- t(x)
  r <- t(t(x) * y)
  r[is.nan(r)] <- 1 # NaN -> 1
  isinf <- is.infinite(rowSums(r))
  rinf <- r[isinf, , drop = FALSE]
  r[isinf, ] <- ifelse(is.finite(rinf), 0, 1)
  compositions::clo(r)
}

power <- function(x, scale) {
  if (is.vector(x)) x <- t(x)
  compositions::clo(exp(scale * log(x)))
}

# Kernel from Aitchison
kn_ait <- function(x, h) {
  a <- x / h + 1
  kn <- list(alpha = a, h = h, X = x)
  class(kn) <- c("Aitchison", "kn")
  kn
}

# Kernel from Chacon(1)
kn_cc1 <- function(x, h) {
  D <- length(x)
  a <- rep(1 / (h * D), D)
  kn <- list(alpha = a, h = h, X = x)
  class(kn) <- c("Chacon1", "kn")
  kn
}

kn_cc2 <- function(x, h) {
  D <- length(x)
  a <- 1 + rep(1 / D, D)
  kn <- list(alpha = a, h = h, X = x)
  class(kn) <- c("Chacon2", "kn")
  kn
}

dkn <- function(x, kn) {
  type <- class(kn)[1]
  alpha <- kn$alpha
  X <- kn$X
  h <- kn$h
  D <- length(X)
  d <- switch(type,
    "Aitchison" = ddir(x, alpha),
    "Chacon1" = ddir(perturbe(x, 1/X), alpha),
    "Chacon2" = ddir(power(perturbe(x, 1/X), 1/h), alpha) / h^(D - 1)
  )
  d
}

rkn <- function(n, kn) {
  type <- class(kn)[1]
  data <- gtools::rdirichlet(n, kn$alpha)
  h <- kn$h
  X <- kn$X
  switch(type,
         "Aitchison" = data,
         "Chacon1" = perturbe(data, X),
         "Chacon2" = perturbe(power(data, h), X)
  )
}

plot.kn <- function(kn, data = NULL, col = NULL, add = FALSE, prop = TRUE,
                    lv = NULL, plv = seq(0.5, 0.95, len = 4),
                    x1 = c(0, 1), x2 = c(0, 2 / sqrt(3)), axl = 100) {
  x <- simplex_grid(x1, x2, axl) # Grid for calculating density
  d <- rep(NA, length = axl^2)
  d[x$valid] <- dkn(x$grid[x$valid, ], kn)
  if (prop) {
    if (is.null(data)) data <- rkn(1000, kn) # Generate 1000 data randomly
    ddata <- dkn(data, kn)
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
  points_simplex(kn$X, add = TRUE)
}

