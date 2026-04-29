# Partial derivative of Dirichlet mixture
# x=0, a>1: f(x)=0, (a-1)/x=Inf => pd=0 (all dim)
# x>0, a=1: f(x)=., (a-1)/x=0 => pd=0 (some dim)
# x=0, a=1: f(x)=., (a-1)/x=NaN => pd=0 (some dim)
pddirmix <- function(x, mix, proj = TRUE) {
  if (is.vector(x)) x <- t(x)
  n <- nrow(x)
  m <- length(mix$pr)
  logpr <- rep(log(mix$pr), each = n)
  logd <- ddir(x, mix$alpha, log = TRUE, outer = TRUE) + logpr
  rx <- x[rep(1:n, each = m), , drop = FALSE]
  ra <- mix$alpha[rep(1:m, n), , drop = FALSE]
  logg <- log(ra - 1) - log(rx)
  logg[is.nan(logg)] <- NA
  r <- logg + c(t(logd))
  pd <- rowsum(exp(r), rep(1:n, each = m))
  if (proj) pd <- pd - rowMeans(pd, na.rm = TRUE)
  pd[is.na(pd)] <- 0
  pd
}

# Partial derivative of Dirichlet kernel
pddirkn <- function(x, kn, proj = TRUE) {
  if (is.vector(x)) x <- t(x)
  n <- nrow(x)
  m <- nrow(kn$alpha)
  logd <- ddir(x, kn$alpha, log = TRUE, outer = TRUE)
  rx <- x[rep(1:n, each = m), , drop = FALSE]
  ra <- kn$alpha[rep(1:m, n), , drop = FALSE]
  logg <- log(ra - 1) - log(rx)
  logg[is.nan(logg)] <- NA
  r <- logg + c(t(logd))
  pd <- rowsum(exp(r), rep(1:n, each = m)) / m
  if (proj) pd <- pd - rowMeans(pd, na.rm = TRUE)
  pd[is.na(pd)] <- 0
  pd
}

# Partial derivative
pddir <- function(x, obj, proj = TRUE) {
  if (class(obj) == "dirkn") pddirkn(x, obj, proj) else pddirmix(x, obj, proj)
}

# Second partial derivative of Dirichlet kernel
spddirkn0 <- function(x, kn, proj = TRUE) {
  if (is.vector(x)) x <- t(x)
  D <- ncol(x)
  n <- nrow(x)
  m <- nrow(kn$alpha)
  d <- ddir(x, kn$alpha, outer = TRUE)
  rx <- x[rep(1:n, each = m), , drop = FALSE]
  ra <- kn$alpha[rep(1:m, n), , drop = FALSE]
  outra <- apply(ra, 1, function(r) {
    mat <- outer(r - 1, r - 1)
    diag(mat) <- (r - 1) * (r - 2)
    mat
  })
  outx <- apply(x, 1, function(r) {
    outer(r, r)
  })[, rep(1:n, each = m), drop = FALSE]
  D2 <- outra / outx

  P <- diag(D) - 1/D
  P2D2 <- t(kronecker(P, P) %*% D2)
  P2D2 * c(d)
}

# Second partial derivative of Dirichlet kernel
spddirkn <- function(x, kn, proj = TRUE) {
  if (is.vector(x)) x <- t(x)
  D <- ncol(x)
  n <- nrow(x)
  m <- nrow(kn$alpha)
  d <- ddir(x, kn$alpha, outer = TRUE)
  rx <- x[rep(1:n, each = m), , drop = FALSE]
  ra <- kn$alpha[rep(1:m, n), , drop = FALSE]
  gm <- exp(log(ra - 1) - log(rx))
  g2m <- t(apply(gm, 1, function(g) kronecker(g, g)))
  Dg <- exp(log(ra - 1) - 2 * log(rx))
  Dg <- -t(apply(Dg, 1, function(r) c(diag(r))))
  D2 <- g2m + Dg # g2m + Dg == t(D2) (in previous function)
  # x_1 a_1, ..., x_1 a_D, ..., x_D a_D, x_D a_D.
  if (proj) {
    P <- diag(D) - 1/D
    P2 <- kronecker(P, P)
    D2 <- t(P2 %*% t(D2))
  }
  r <- c(t(d)) * D2
  rowsum(r, rep(1:n, each = m)) / m
}
