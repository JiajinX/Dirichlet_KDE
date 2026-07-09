SCV_2nd <- function(h, x, mix) {
  n <- nrow(x)
  ya <- x/h + 1
  Cx <- Cdir(ya, log = TRUE)
  outCay <- outer(Cx, attr(mix$alpha, "logmbeta"), "+")
  Cay <- Cdir_plus(ya, mix$alpha, log = TRUE)
  Cd <- exp(outCay - Cay)
  mean(Cd %*% mix$pr)
}

# The second integral has close-form result if we know true Dirichlet mixture
SCVf <- function(h, x, mix, verbose = FALSE) {
  if (h <= 0) stop("h must be positive.")
  if (any(abs(rowSums(x) - 1) > sqrt(.Machine$double.eps))) stop("x must be on the simplex.")
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
  mixa <- mix$alpha
  m <- nrow(mixa)

  a <- x / h + c
  u <- a - 1
  v <- mixa - 1

  Ca <- Cdir(a, log = TRUE)
  Cmixa <- attr(mixa, "logmbeta")
  if (is.null(Cmixa)) Cmixa <- Cdir(mixa, log = TRUE)

  sum_lgamma <- S1 <- S2 <- Tu <- Tv <- matrix(0, n, m)
  for (k in seq_len(D)) {
    w <- outer(u[, k], v[, k], "+")
    uv <- tcrossprod(u[, k], v[, k])
    iw <- 1 / w

    sum_lgamma <- sum_lgamma + lgamma(w + 1)
    S1 <- S1 + uv / (w * (w - 1))
    S2 <- S2 + uv / w^2
    Tu <- Tu + u[, k] * iw
    Tv <- Tv + sweep(iw, 2, v[, k], "*")
  }

  Cs <- (1 - 1 / D) * S1 + S2 / D - Tu * Tv / D

  suma <- outer(rowSums(a), rowSums(mixa), "+")
  logCd <- outer(Ca, Cmixa, "+") - lgamma(suma - D) + sum_lgamma

  pref <- suma - D - 1
  H <- pref * (pref - 1) * Cs

  mean((exp(logCd) * H) %*% mix$pr)
}

# The second integral has close-form result if we know true Dirichlet mixture
DSCVf <- function(h, x, mix, c = 1.5, verbose = FALSE) {
  if (h <= 0) stop("h must be positive.")
  if (any(abs(rowSums(x) - 1) > sqrt(.Machine$double.eps))) stop("x must be on the simplex.")
  if (any(x <= 0)) stop("x must be positive.")
  term1 <- DCV_1st(h, x, c) # First term
  term2 <- DSCV_2nd(h, x, mix, c) # Second term
  re <- term1 - 2 * term2
  if (verbose) print(paste0(round(h, 4), ": ", round(re, 2)))
  re
}
