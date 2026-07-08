CV_1st <- function(h, x) {
  n <- nrow(x)
  ya <- x/h + 1
  Cx <- Cdir(ya, log = TRUE)
  outCxx <- outer(Cx, Cx, "+")
  Cxx <- Cdir_plus(ya, ya, log = TRUE)
  exp(log_sum_exp(outCxx - Cxx) - 2 * log(n))
}

LSCV_2nd <- function(h, x) {
  n <- nrow(x)
  ya <- x/h + 1
  douter <- ddir(x, ya, log = TRUE, outer = TRUE)
  diag(douter) <- -Inf
  exp(log_sum_exp(douter) - log(n * (n - 1)))
}

LSCVf <- function(h, x, verbose = FALSE) {
  if (any(!is.finite(x))) stop("x must be finite.")
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
DLSCV_1st <- function(h, x, c = 1.5) {
  n <- nrow(x)
  D <- ncol(x)
  ya <- x/h + c
  Cx <- Cdir(ya, log = TRUE)
  outCxx <- outer(Cx, Cx, "+")
  Cxx <- Cdir_plus(ya, ya, log = TRUE)
  logCd <- outCxx - Cxx
  
  idx1 <- unlist(lapply(1:n, function(i) i:n))
  idx2 <- rep(1:n, n:1)
  diag <- idx1 == idx2
  # Indices of lower triangular elements in symmetric matrix
  log_ya_1 <- log(ya - 1)
  Csup <- log_ya_1[idx1, , drop = FALSE] + log_ya_1[idx2, , drop = FALSE]
  sumya <- ya[idx1, , drop = FALSE] + ya[idx2, , drop = FALSE]
  log_sumya_2 <- log(sumya - 2)
  log_sumya_3 <- log(sumya - 3)
  Cs1 <- rowSums(exp(Csup - (log_sumya_2 + log_sumya_3)))
  Cs2 <- rowSums(exp(Csup - 2 * log_sumya_2))
  tx <- rowSums(exp(log_ya_1[idx1, , drop = FALSE] - log_sumya_2))
  ty <- rowSums(exp(log_ya_1[idx2, , drop = FALSE] - log_sumya_2))
  Cscut <- tx * ty / D
  Cs <- (1 - 1/D) * Cs1 + Cs2/D - Cscut
  log_r1 <- log(abs(Cs)) + logCd[lower.tri(logCd, diag = TRUE)] +
    log(ifelse(diag, 1, 2))
  r1 <- signed_log_sum_exp(log_r1, sign(Cs))
  logC <- prefactor(1/h + c * D, 1/h + c * D, D)
  unname(r1["sign"] * exp(r1["logabs"] + logC - 2 * log(n)))
}

DLSCV_2nd <- function(h, x, c = 1.5) {
  n <- nrow(x)
  D <- ncol(x)
  ya <- x / h + c
  
  log_douter <- ddir(x, ya, log = TRUE, outer = TRUE)
  diag2 <- seq(1, n * n, by = n + 1)
  
  rx <- x[rep(1:n, each = n - 1), , drop = FALSE]
  rya <- ya[rep(1:n, n)[-diag2], , drop = FALSE]
  
  g <- (rya - 1) / rx
  q <- (rya - 1) / rx^2
  
  bracket <- rowSums((g - rowMeans(g))^2) - (1 - 1 / D) * rowSums(q)
  contrib <- exp(c(t(log_douter))[-diag2]) * bracket
  mean(contrib)
}

DLSCVf <- function(h, x, c = 1.5, verbose = FALSE) {
  if (any(!is.finite(x))) stop("x must be finite.")
  if (any(x <= 0)) stop("x must be strictly positive.")
  term1 <- DLSCV_1st(h, x, c) # First term
  term2 <- DLSCV_2nd(h, x, c) # Second term
  re <- term1 + 2 * term2
  if (verbose) print(paste0(round(h, 4), ": ", round(re, 2)))
  re
}

# DLSCVf <- function(h, x, c = 1.5, verbose = FALSE) {
#   if (any(!is.finite(x))) stop("x must be finite.")
#   if (any(x <= 0)) stop("x must be strictly positive.")
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
