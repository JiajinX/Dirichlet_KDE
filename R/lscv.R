LSCVf <- function(h, x, verbose = FALSE) {
  n <- nrow(x)
  Cx <- Cdir(x/h + 1, log = TRUE)
  outCxx <- outer(Cx, Cx, "+")
  Cxx <- Cdir_plus(x / h + 1, x / h + 1, log = TRUE)
  term1 <- sum(exp(outCxx - Cxx)) / n^2
  douter <- ddir(x, x/h + 1, outer = TRUE)
  diag(douter) <- 0
  term2 <- sum(douter) / (n * (n - 1))
  re <- term1 - 2 * term2
  if (verbose) print(paste0(round(h, 4), ": ", round(re, 2)))
  re
}

# The second integral has close-form result if we know true Dirichlet mixture
LSCVcf <- function(h, x, mix, verbose = FALSE) {
  n <- nrow(x)
  ya <- x/h + 1
  Cx <- Cdir(ya, log = TRUE)
  outCxx <- outer(Cx, Cx, "+")
  Cxx <- Cdir_plus(x / h + 1, x / h + 1, log = TRUE)
  # First term
  term1 <- sum(exp(outCxx - Cxx)) / n^2
  # Second term
  outCay <- outer(Cx, attr(mix$alpha, "logmbeta"), "+")
  Cay <- Cdir_plus(ya, mix$alpha, log = TRUE)
  Cd <- exp(outCay - Cay)
  term2 <- sum(mix$pr * t(Cd)) / n
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

LSCVmc <- function(h, x, n = 1e5) {
  fhat <- dirmix(x/h + 1)
  samp <- rdirmix(n, fhat)
  kn <- dirkn(samp, h)
  term1 <- mean(ddirkn(samp, kn))
  douter <- ddir(x, x/h + 1, outer = TRUE)
  diag(douter) <- 0
  n <- nrow(x)
  term2 <- sum(douter) / (n * (n - 1))
  term1 - 2 * term2
}

# Integral of the frist term in DLSCV
DLSCV_1st <- function(h, x, c = 1.5) {
  n <- nrow(x)
  D <- ncol(x)
  ya <- x/h + c
  Cx <- Cdir(ya, log = TRUE)
  outCxx <- outer(Cx, Cx, "+")
  Cxx <- Cdir_plus(ya, ya, log = TRUE)
  Cd <- exp(outCxx - Cxx)

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
  r1 <- Cs * Cd[lower.tri(Cd, diag = TRUE)]
  r1 <- 2 * sum(r1[!diag]) + sum(r1[diag])
  logC <- prefactor(1/h + c * D, 1/h + c * D, D)
  r1 * exp(logC) / n^2
}

DLSCVf <- function(h, x, c = 1.5, verbose = FALSE) {
  n <- nrow(x)
  D <- ncol(x)
  ya <- x/h + c

  # First term
  term1 <- DLSCV_1st(h, x, c)
  # Second term
  douter <- ddir(x, ya, outer = TRUE)
  diag(douter) <- 0
  diag2 <- seq(1, n * n, by = n + 1)
  # Delete rows of diagonal index, because diag(douter) = 0
  rx <- x[rep(1:n, each = n - 1), , drop = FALSE]
  rya <- ya[rep(1:n, n)[-diag2], , drop = FALSE]
  log_rya_1 <- log(rya - 1)
  log_rx <- log(rx)
  gm <- exp(log_rya_1 - log_rx)
  idx3 <- unlist(lapply(1:D, function(i) i:D))
  idx4 <- rep(1:D, D:1)
  diag3 <- idx3 == idx4
  D2 <- t(apply(gm, 1, function(g) g[idx3] * g[idx4]))
  Dg <- exp(log_rya_1 - 2 * log_rx)
  D2[, diag3] <- D2[, diag3] - Dg
  r2 <- c(t(douter))[-diag2] * D2
  r2 <- colSums(r2) / (n * (n - 1))
  P <- diag3 - 1/D
  Pr2 <- -P * r2
  term2 <- 2 * sum(Pr2[!diag3]) + sum(Pr2[diag3])
  re <- c(term1 - 2 * term2)
  if (verbose) print(paste0(round(h, 4), ": ", round(re, 2)))
  re
}

# The second integral has close-form result if we know true Dirichlet mixture
DLSCVcf <- function(h, x, mix, c = 1.5, verbose = FALSE) {
  n <- nrow(x)
  D <- ncol(x)
  ya <- x/h + c
  Cx <- Cdir(ya, log = TRUE)
  # First term
  term1 <- DLSCV_1st(h, x, c)
  # Second term
  a0 <- rowSums(mix$alpha)
  outCay <- outer(Cx, attr(mix$alpha, "logmbeta"), "+")
  Cay <- Cdir_plus(ya, mix$alpha, log = TRUE)
  Cd <- exp(outCay - Cay)
  Cs <- D1S(ya, mix$alpha)
  logC <- prefactor(1/h + c * D, a0, D)
  r <- Cd * Cs * exp(logC[rep(1, n), ])
  term2 <- sum(mix$pr * t(r)) / n
  re <- c(term1 - 2 * term2)
  if (verbose) print(paste0(round(h, 4), ": ", round(re, 2)))
  re
}

# DLSCVmc <- function(h, x, c, mix, n = 1e4) {
#   samp <- rdirmix(n, dirmix(a = rep(1,3))) # Sample from uniform
#   kn <- dirkn(x, h = h, c) # f_hat
#   dfhat <- pddir(samp, obj = kn)
#   dfm <- pddir(samp, obj = mix)
#   term1 <- mean(rowSums(dfhat * dfhat)) / 2
#   term2 <- mean(rowSums(dfhat * dfm)) / 2
#   # cat(term1, " ", term2, "\n")
#   term1 - 2 * term2
# }
