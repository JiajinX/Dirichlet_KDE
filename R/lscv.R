CV_1st <- function(h, x) {
  n <- nrow(x)
  ya <- x / h + 1
  Cx <- Cdir(ya, log = TRUE)
  # Lower-triangular pairs
  idx1 <- unlist(lapply(seq_len(n), function(i) i:n), use.names = FALSE)
  idx2 <- rep(seq_len(n), n:1)
  diag <- idx1 == idx2

  aij <- ya[idx1, , drop = FALSE] + ya[idx2, , drop = FALSE] - 1
  log_terms <- Cx[idx1] + Cx[idx2] - Cdir(aij, log = TRUE) + ifelse(diag, 0, log(2))
  exp(log_sum_exp(log_terms) - 2 * log(n))
}

LSCV_2nd <- function(h, x) {
  n <- nrow(x)
  u <- x / h
  ya <- u + 1
  logx <- log(x)
  has_zero <- any(x == 0)
  if (has_zero) {
    zero_x <- x == 0
    logx[zero_x] <- 0 # Temporarily treat 0^0 as one
  }
  log_douter <- tcrossprod(logx, u)
  if (has_zero) {
    # x[i,k] = 0 and u[j,k] > 0 implies zero density
    invalid <- tcrossprod(zero_x * 1, (u > 0) * 1) > 0
    log_douter[invalid] <- -Inf
  }
  log_douter <- sweep(log_douter, 2, Cdir(ya, log = TRUE), "+")
  diag(log_douter) <- -Inf
  exp(log_sum_exp(log_douter) - log(n * (n - 1)))
}

LSCVf <- function(h, x, verbose = FALSE) {
  if (h <= 0) stop("h must be positive.")
  if (any(abs(rowSums(x) - 1) > sqrt(.Machine$double.eps))) stop("x must be on the simplex.")
  if (any(x < 0)) stop("x must be non-negative.")
  term1 <- CV_1st(h, x)
  term2 <- LSCV_2nd(h, x)
  re <- term1 - 2 * term2
  if (verbose) print(paste0(round(h, 4), ": ", round(re, 2)))
  re
}

# LSCV0 <- function(h, x) {
#   n <- nrow(x)
#   D <- ncol(x)
#   gh <- lgamma(1/h + D)
#   g2h <- lgamma(2/h + D)
#   gx <- lgamma(x/h + 1)
#   idx1 <- unlist(lapply(1:n, function(i) i:n))
#   idx2 <- rep(1:n, n:1)
#   diag <- idx1 == idx2
#
#   x1 <- x[idx1, , drop = FALSE]
#   x2 <- x[idx2, , drop = FALSE]
#   xx <- x1 + x2
#   gxx <- lgamma(xx/h + 1)
#   gx1 <- gx[idx1, , drop = FALSE]
#   gx2 <- gx[idx2, , drop = FALSE]
#   outergx <- gx1 + gx2
#   term1 <- matrix(NA, n, n)
#   term1[lower.tri(term1, diag = TRUE)] <- exp(gh - g2h + rowSums(gxx - outergx))
#   term1[upper.tri(term1)] <- term1[lower.tri(term1)]
#
#   x1 <- x[rep(seq_len(n), n), ]
#   x2 <- x[rep(seq_len(n), each = n), ]
#   term2 <- 2 * (1 - 1/n)^(-1) * exp(rowSums(x2/h * log(x1) - lgamma(x2/h + 1)))
#   dim(term2) <- c(n, n)
#   diag(term2) <- 0
#   exp(gh - 2 * log(n)) * sum(term1 - term2)
# }

# Integral of the frist term in DLSCV
# If x = 0, then function requires c > 1.5.
DCV_1st <- function(h, x, c = 1.5) {
  n <- nrow(x)
  D <- ncol(x)
  ya <- x/h + c
  u <- ya - 1
  Cx <- Cdir(ya, log = TRUE)

  # Only calculate lower-triangular pairs
  idx1 <- unlist(lapply(seq_len(n), function(i) i:n), use.names = FALSE)
  idx2 <- rep(seq_len(n), n:1)
  diag <- idx1 == idx2
  ui <- u[idx1, , drop = FALSE]
  uj <- u[idx2, , drop = FALSE]
  usum <- ui + uj
  ri <- ui / usum
  rj <- uj / usum
  Cs1 <- rowSums(ri * uj / (usum - 1))
  Cs2 <- rowSums(ri * rj)
  Cs <- (1 - 1 / D) * Cs1 + Cs2 / D - rowSums(ri) * rowSums(rj) / D
  logCd <- Cx[idx1] + Cx[idx2] - Cdir(usum + 1, log = TRUE)
  log_terms <- log(abs(Cs)) + logCd + ifelse(diag, 0, log(2))
  total <- signed_log_sum_exp(log_terms, sign(Cs))

  alpha0 <- 1 / h + c * D
  pref <- 2 * alpha0 - D - 1
  log_prefactor <- log(pref) + log(pref - 1)

  unname(total["sign"] * exp(total["logabs"] + log_prefactor - 2 * log(n)))
}

DLSCV_2nd <- function(h, x, c = 1.5) {
  n <- nrow(x)
  D <- ncol(x)
  ya <- x / h + c
  u <- ya - 1

  invx <- 1 / x
  invx2 <- invx^2

  # log f_Dir(x_i; ya_j)
  log_douter <- tcrossprod(log(x), u)
  log_douter <- sweep(log_douter, 2, Cdir(ya, log = TRUE), "+")

  # sum(g), sum(g^2), sum(q) for every (i, j)
  sum_g <- tcrossprod(invx, u)
  sum_g2 <- tcrossprod(invx2, u^2)
  sum_q <- tcrossprod(invx2, u)

  bracket <- sum_g2 - sum_g^2 / D - (1 - 1 / D) * sum_q

  diag(log_douter) <- -Inf
  diag(bracket) <- 0

  sum(exp(log_douter) * bracket) / (n * (n - 1))
}

DLSCVf <- function(h, x, c = 1.5, verbose = FALSE) {
  if (h <= 0) stop("h must be positive.")
  if (any(abs(rowSums(x) - 1) > sqrt(.Machine$double.eps))) stop("x must be on the simplex.")
  if (any(x <= 0)) stop("x must be positive.")
  term1 <- DCV_1st(h, x, c) # First term
  term2 <- DLSCV_2nd(h, x, c) # Second term
  re <- term1 + 2 * term2
  if (verbose) print(paste0(round(h, 4), ": ", round(re, 2)))
  re
}

# DLSCVf <- function(h, x, c = 1.5, verbose = FALSE) {
#   if (any(!is.finite(x))) stop("x must be finite.")
#   if (any(x < 0)) stop("x must be non-negative.")
#   n <- nrow(x)
#   D <- ncol(x)
#   ya <- x/h + c
#
#   # First term
#   term1 <- DLSCV_1st(h, x, c)
#   # Second term
#   log_douter <- ddir(x, ya, log = TRUE, outer = TRUE)
#   diag2 <- seq(1, n * n, by = n + 1)
#   # Delete rows of diagonal index.
#   rx <- x[rep(1:n, each = n - 1), , drop = FALSE]
#   rya <- ya[rep(1:n, n)[-diag2], , drop = FALSE]
#   log_rya_1 <- log(rya - 1)
#   log_rx <- log(rx)
#   log_gm <- log_rya_1 - log_rx # log(g)
#
#   idx3 <- unlist(lapply(1:D, function(i) i:D))
#   idx4 <- rep(1:D, D:1)
#   diag3 <- idx3 == idx4 # diag3: diagonal columns == TRUE
#   # g^{\otimes 2} = gg^{\top}: D × D matrix
#   # Upper triangle element of log_D2 = log(gg^{\top})
#   # Each column is log(g_i * g_j)
#   log_D2 <- matrix(NA, nrow = nrow(log_gm), ncol = length(idx3))
#   log_D2[, !diag3] <- log_gm[, idx3[!diag3], drop = FALSE] +
#     log_gm[, idx4[!diag3], drop = FALSE] # Off-diagonal columns
#   log_D2[, diag3] <- log_rya_1 + log(abs(rya - 2)) - 2 * log_rx
#   # Diagonal columns: log((a_k - 1)(a_k - 2)/x_k^2)
#   sign_D2 <- matrix(1, nrow = nrow(log_D2), ncol = ncol(log_D2))
#   sign_D2[, diag3] <- sign(rya - 2) # Records possible negative signs on diagonal terms
#
#   log_terms <- log_D2 + c(t(log_douter))[-diag2] # -diag2 for leave one out
#   # log(D^{\otimes 2} f_h(x_i; x_j/h + c)): n(n-1) × D(D+1)/2 [upper triangle elements]
#   P <- diag3 - 1/D
#   coef <- ifelse(diag3, 1, 2) * -P
#   # Since only the elements of the upper triangle are recorded,
#   # off-diagonal elements * 2.
#   log_terms <- t(t(log_terms) + log(abs(coef)))
#   sign_terms <- t(t(sign_D2) * sign(coef))
#   term2 <- signed_log_sum_exp(c(log_terms), c(sign_terms))
#   term2 <- unname(term2["sign"] * exp(term2["logabs"] - log(n * (n - 1))))
#
#   re <- c(term1 - 2 * term2)
#   if (verbose) print(paste0(round(h, 4), ": ", round(re, 2)))
#   attr(re, "term") <- c(term1, term2)
#   re
# }
