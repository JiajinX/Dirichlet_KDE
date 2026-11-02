plot.meanshift <- function(clu, st, pth = FALSE, obj = NULL,
                           txt = NULL, contour = FALSE,
                           col = NULL, lab = NULL, txt.cex = 0.5) {
  m <- max(clu$cluster)
  if (is.null(col)) col <- paletteer::paletteer_d("ggthemes::gdoc", m)
  simplexpanel(lab, labcex = 1.2)
  if (contour) plot(obj, data = st, add = TRUE, col = "grey60",
                    lv = seq(30, 90, len = 4) / 100, axl = 200)
  for (i in seq_len(m)) {
    sti <- st[clu$cluster == i, ]
    if (pth) MEM(sti, obj, pth = TRUE, add = TRUE, col = col[i])
    if (is.null(txt)) {
      points_simplex(sti, col = col[i], add = TRUE)
    } else {
      txti <- txt[clu$cluster == i]
      text_simplex(sti, txti, col = col[i], add = TRUE, cex = txt.cex)
    }
    points_simplex(clu$centers, add = TRUE, pch = 3, cex = 1.2, col = 1)
  }
}

simplex_trigrid <- function(n) {
  m <- 2^n
  # Generate all non-negative integer triples (i,j,k)
  # satisfying i + j + k = s
  triples <- function(s) {
    len <- (s + 1):1
    i <- rep.int(0:s, times = len)
    j <- sequence(len) - 1
    k <- s - i - j
    cbind(i, j, k)
  }
  # Triangles with the same orientation as the simplex
  up <- triples(m - 1)
  cup <- (up + 1 / 3) / m
  # Inverted triangles
  down <- triples(m - 2)
  cdown <- (down + 2 / 3) / m
  centers <- rbind(cup, cdown)
  ort <- c(rep(1, nrow(cup)), rep(-1, nrow(cdown)))
  list(centers = centers, scale = 1 / (sqrt(3) * m), ort = ort)
}

# Triangle points of grid centers.
simplex_tripoint <- function(centers, scale, ort, tol = 1e-10) {
  if (is.vector(centers)) centers <- t(centers)
  n <- nrow(centers)
  if (length(ort) == 1) ort <- rep(ort, n)
  if (length(ort) != n) stop("'ort' has the wrong length.")
  B <- rbind(c( 2, -1, -1), c(-1,  2, -1), c(-1, -1,  2))
  B <- B * scale / sqrt(3)
  tri <- outer(ort, B, "*")
  tri <- sweep(outer(ort, B, "*"), c(1, 3), centers, "+")
  tri[abs(tri) < tol] <- 0
  tri[abs(tri - 1) < tol] <- 1
  tri
}

polygon_simplex <- function(x, ...) {
  xy <- cbind((2 - 2 * x[, 2] - x[, 1]) / sqrt(3), x[, 1])
  polygon(xy[, 1], xy[, 2], ...)
}

plot_cluster_region <- function(centers, cluster, scale, ort, kn,
                                alpha = c(0, 1), gamma = 0.5, qmax = 1,
                                lab = NULL, cols = NULL) {
  if (is.vector(centers)) centers <- t(centers)
  n <- nrow(centers)
  if (length(cluster) != n) stop("Incorrect cluster number.")
  tri <- simplex_tripoint(centers, scale, ort)
  m <- max(cluster)
  d <- ddirkn(centers, kn)
  dmax <- quantile(d, qmax, names = FALSE)
  r <- pmin(d / dmax, 1)^gamma
  a <- alpha[1] + (alpha[2] - alpha[1]) * r
  if (is.null(cols)) cols <- hcl.colors(m, "Set 2")
  rgb <- col2rgb(cols[cluster]) / 255
  fill <- rgb(rgb[1, ], rgb[2, ], rgb[3, ], alpha = a)
  simplexpanel(lab, 1.2)
  # Fill each small triangle
  for (i in seq_len(n)) {
    polygon_simplex(tri[i, , ], col = fill[i], border = NA)
  }
}
