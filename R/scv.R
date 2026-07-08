SCV_2nd <- function(h, x, mix) {
  n <- nrow(x)
  ya <- x/h + 1
  Cx <- Cdir(ya, log = TRUE)
  outCay <- outer(Cx, attr(mix$alpha, "logmbeta"), "+")
  Cay <- Cdir_plus(ya, mix$alpha, log = TRUE)
  Cd <- exp(outCay - Cay)
  sum(mix$pr * t(Cd)) / n
}

# The second integral has close-form result if we know true Dirichlet mixture
SCVf <- function(h, x, mix, verbose = FALSE) {
  if (any(!is.finite(x))) stop("x must be finite.")
  if (any(x < 0)) stop("x must be non-negative.")
  term1 <- CV_1st(h, x)
  term2 <- SCV_2nd(h, x, mix)
  re <- term1 - 2 * term2
  if (verbose) print(paste0(round(h, 4), ": ", round(re, 2)))
  re
}

DSCV_2nd <- function(h, x, mix, c = 1.5) {
  n <- nrow(x)
  D <- ncol(x)
  ya <- x/h + c
  Cx <- Cdir(ya, log = TRUE)
  a0 <- rowSums(mix$alpha)
  outCay <- outer(Cx, attr(mix$alpha, "logmbeta"), "+")
  Cay <- Cdir_plus(ya, mix$alpha, log = TRUE)
  Cd <- exp(outCay - Cay)
  Cs <- D1S(ya, mix$alpha)
  logC <- c(prefactor(1/h + c * D, a0, D))
  r <- Cd * Cs * exp(logC[rep(1, n), ])
  sum(mix$pr * t(r)) / n
}

# The second integral has close-form result if we know true Dirichlet mixture
DSCVf <- function(h, x, mix, c = 1.5, verbose = FALSE) {
  if (any(!is.finite(x))) stop("x must be finite.")
  if (any(x < 0)) stop("x must be non-negative.")
  term1 <- DCV_1st(h, x, c) # First term
  term2 <- DSCV_2nd(h, x, mix, c) # Second term
  re <- term1 - 2 * term2
  if (verbose) print(paste0(round(h, 4), ": ", round(re, 2)))
  re
}
