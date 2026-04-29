# Plot for real data
repoint <- function(data, x1, x2, group = NULL) {
  codata <- cbind(data[x1], data[x2], 1 - data[x1] - data[x2])
  col <- paletteer::paletteer_d("ggthemes::excel_Droplet")
  simplexpanel()
  if (!is.null(group)) {
    unig <- unique(group)
    for (i in seq_along(unig)) {
      points_simplex(codata[group == unig[i], ], col = col[i], add = TRUE)
    }
  } else {
    points_simplex(codata, add = TRUE)
  }
}

par_repoint <- function(data, group = NULL, nr = 4) {
  comb <- t(combn(colnames(data), m = 2))
  l <- nrow(comb)
  op <- par(no.readonly = TRUE)
  par(mfrow = c(ceiling(l / nr), nr))
  for (i in seq_len(l)) repoint(data, x1 = comb[i, 1], x2 = comb[i, 2], group)
  par(op)
}

reden <- function(obj, x1, x2) {
  vx1 <- obj$x[, x1]
  vx2 <- obj$x[, x2]
  vx3 <- 1 - vx1 - vx2
  x <- cbind(vx1, vx2, vx3)
  obj <- dirkn(x, h = obj$h, c = obj$c)
  plot(obj)
}

par_reden <- function(kn, nr = 4) {
  comb <- t(combn(colnames(kn$x), m = 2))
  l <- nrow(comb)
  op <- par(no.readonly = TRUE)
  par(mfrow = c(ceiling(l / nr), nr))
  for (i in seq_len(l)) reden(kn, x1 = comb[i, 1], x2 = comb[i, 2])
  par(op)
}
