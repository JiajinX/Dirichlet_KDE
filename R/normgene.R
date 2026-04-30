# Generate random data from Gaussian mixture
rnormix <- function(n, mix) {
  ks::rmvnorm.mixt(n, mix$mu, mix$Sigma, props = mix$pr)
}

# Density of Gaussian mixture on Euclidean space
dnormix <- function(x, mix, log = FALSE) {
  if (is.vector(x)) x <- t(x)
  d <- ks::dmvnorm.mixt(x, mix$mu, mix$Sigma, props = mix$pr)
  if (log) d <- log(d)
  d
}

# inv = FALSE, Simplex -> R
# inv = TRUE, R -> Simplex
ilrts <- function(x) {
  if (is.vector(x)) x <- t(x)
  inv <- !all(round(rowSums(x), 10) == 1)
  if (inv) {
    rx <- x
    sx <- compositions::ilrInv(rx)
    attributes(sx) <- NULL
    dim(sx) <- dim(rx) + c(0, 1)
  } else {
    sx <- x
    if (any(sx <= 0)) stop("All components of x must be positive.")
    rx <- compositions::ilr(sx)
    attributes(rx) <- NULL
    dim(rx) <- dim(sx) - c(0, 1)
  }
  jacb <- exp(rowSums(log(sx)) + 1 / 2 * log(ncol(sx)))
  x <- list(sx = sx, rx = rx, jacb = jacb)
  class(x) <- "ilrx"
  x
}

# Extract x by index
select_ilrx <- function(x, idx) {
  if (!inherits(x, "ilrx")) stop("The class of x should be ilrx")
  x$sx <- x$sx[idx, , drop = FALSE]
  x$rx <- x$rx[idx, , drop = FALSE]
  x$jacb <- x$jacb[idx]
  x
}

# Density of Gaussian mixture on Simplex
dSnormix <- function(x, mix, log = FALSE) {
  if (!inherits(x, "ilrx")) x <- ilrts(x)
  logd <- dnormix(x$rx, mix, log = TRUE)
  ds <- logd - log(x$jacb)
  if (!log) ds <- exp(ds)
  ds
}

# Partial derivative of Gaussian mixture on Simplex
pdnormix <- function(x, mix, proj = TRUE) {
  if (is.vector(x)) x <- t(x)
  n <- nrow(x)
  D <- ncol(x)
  m <- nrow(mix$mu)
  d <- D - 1

  ilrx <- ilrts(x)
  y <- ilrx$rx
  jacb <- ilrx$jacb
  invx <- 1/x
  V <- compositions::ilrBase(D = D)

  tvz <- matrix(nrow = n * m, ncol = D)
  dy <- matrix(nrow = m, ncol = n)
  idxr <- rep(seq_len(m), n)
  for (i in seq_len(m)) {
    mu <- mix$mu[i, ]
    sig <- mix$Sigma[((i - 1) * d + 1):(i * d), , drop = FALSE]
    dy[i, ] <- mvtnorm::dmvnorm(y, mean = mu, sigma = sig, log = TRUE)
    z <- solve(sig, t(sweep(y, 2, mu, "-")))
    tvz[idxr == i, ] <- t(V %*% z)
  }
  dSx <- exp(t(dy) - log(jacb))
  prdS <- t(dSx) * mix$pr
  pd <- -invx * rowsum(c(prdS) * (tvz + 1), group = rep(seq_len(n), each = m))
  if (proj) pd <- pd - rowMeans(pd, na.rm = TRUE)
  pd
}
