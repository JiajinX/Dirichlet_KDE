Cdir_plus <- function(x, y, log = FALSE) {
  nx <- nrow(x)
  ny <- nrow(y)
  xpy <- x[rep(seq_len(nx), ny), , drop = FALSE] +
    y[rep(seq_len(ny), each = nx), , drop = FALSE]
  Cxpy <- Cdir(xpy - 1, log = log)
  dim(Cxpy) <- c(nx, ny)
  Cxpy
}

gene_compsamp <- function(mix, mn = 1e5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  m <- length(mix$pr)
  samp <- replicate(2, do.call(rbind, rdircomp(mn, mix)), simplify = FALSE)
  attr(samp, "size") <- c(m, mn)
  attr(mix, "samp") <- samp
  mix
}

Cdh <- function(h, mix) {
  samp <- attr(mix, "samp")
  size <- attr(samp, "size")
  Cy <- attr(h, "Cy")
  if (is.list(Cy)) {
    Cy <- Cy[[1]]
    samp <- samp[[1]]
  }

  C2y <- Cdir(2 * samp / h + 1, log = TRUE)
  r <- exp(2 * Cy - C2y)
  dim(r) <- rev(size)
  sum(colMeans(r) * mix$pr)
}

Psi0 <- function(mix) {
  Ca <- attr(mix$alpha, "logmbeta")
  outCaa <- outer(Ca, Ca, "+")
  Caa <- Cdir_plus(mix$alpha, mix$alpha, log = TRUE)
  exp(outCaa - Caa)
}

Psi1 <- function(h, mix) {
  samp <- attr(mix, "samp")
  size <- attr(samp, "size")
  Ca <- attr(mix$alpha, "logmbeta")
  Cy <- attr(h, "Cy")
  if (is.list(Cy)) {
    Cy <- Cy[[1]]
    samp <- samp[[1]]
  }

  outCay <- outer(Cy, Ca, "+")
  Cya <- Cdir_plus(samp / h + 1, mix$alpha, log = TRUE)
  r <- exp(outCay - Cya)
  K <- kronecker(diag(size[1]), matrix(1/size[2], 1, size[2]))
  K %*% r
}

Psi2 <- function(h, mix) {
  samp <- attr(mix, "samp")
  size <- attr(samp, "size")
  Cy <- attr(h, "Cy")

  v <- seq_len(prod(size))
  i <- rep(v, size[1])
  mtx <- matrix(v, nrow = size[2])
  j <- c(mtx[rep(seq_len(size[2]), size[1]), ])

  up <- Cy[[1]][i] + Cy[[2]][j]
  down <- Cdir((samp[[1]][i, ] + samp[[2]][j, ]) / h + 1, log = TRUE)
  r <- exp(up - down)
  dim(r) <- c(prod(size), size[1])
  K <- kronecker(diag(size[1]), matrix(1/size[2], 1, size[2]))
  K %*% r
}

MISEf <- function(h, mix, n, verbose = FALSE) {
  samp <- attr(mix, "samp")
  if (is.null(samp)) stop("No samples generated.")
  attr(h, "Cy") <- lapply(samp, function(xi) Cdir(xi / h + 1, log = TRUE))

  cdh <- Cdh(h, mix)
  psi0 <- Psi0(mix)
  psi1 <- Psi1(h, mix)
  psi2 <- Psi2(h, mix)

  re <- c(cdh / n + t(mix$pr) %*% ((1 - 1/n) * psi2 - 2 * psi1 + psi0) %*% mix$pr)
  if (verbose) print(paste0(round(h, 4), ": ", round(re, 2)))
  re
}
