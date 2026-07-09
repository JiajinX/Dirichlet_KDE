DCdh <- function(h, mix, c) {
  size <- attr(attr(mix, "samp"), "size")
  y <- attr(h, "y")
  if (is.list(y)) y <- y[[1]]
  D <- ncol(mix$alpha)

  ya2 <- 2 * y$a
  C2y <- lgamma(2 / h + 2 * c * D - D) - rowSums(lgamma(ya2 - 1))
  Cd <- exp(2 * y$Ca - C2y)
  Cs <- (1 - 1/D) * rowSums((y$a - 1) / (2 * (ya2 - 3))) + (1 - D) / 4
  r <- Cd * Cs
  dim(r) <- rev(size)

  p <- 2 * (1 / h + c * D) - D - 1
  sum(colMeans(r) * mix$pr) * p * (p - 1)
}

DPsi0 <- function(mix) {
  a <- mix$alpha
  m <- nrow(a)
  D <- ncol(a)
  Ca <- attr(a, "logmbeta")
  if (is.null(Ca)) Ca <- Cdir(a, log = TRUE)

  a0 <- rowSums(a)
  idx1 <- unlist(lapply(seq_len(m), function(i) i:m), use.names = FALSE)
  idx2 <- rep(seq_len(m), m:1)
  u1 <- a[idx1, , drop = FALSE] - 1
  u2 <- a[idx2, , drop = FALSE] - 1
  w <- u1 + u2
  r1 <- u1 / w
  r2 <- u2 / w
  Cs1 <- rowSums(r1 * u2 / (w - 1))
  Cs2 <- rowSums(r1 * r2)
  Cs <- (1 - 1 / D) * Cs1 + Cs2 / D - rowSums(r1) * rowSums(r2) / D
  logCd <- Ca[idx1] + Ca[idx2] - Cdir(w + 1, log = TRUE)
  p <- a0[idx1] + a0[idx2] - D - 1
  value <- exp(logCd + log(p) + log(p - 1)) * Cs

  out <- matrix(0, m, m)
  out[lower.tri(out, diag = TRUE)] <- value
  out[upper.tri(out)] <- t(out)[upper.tri(out)]
  out
}

DPsi1 <- function(h, mix, c) {
  size <- attr(attr(mix, "samp"), "size")
  y <- attr(h, "y")
  if (is.list(y)) y <- y[[1]]
  m <- size[1]
  S <- size[2]
  D <- ncol(mix$alpha)
  a <- y$a
  u <- a - 1
  Ca <- y$Ca

  mixa <- mix$alpha
  Cmixa <- attr(mixa, "logmbeta")
  if (is.null(Cmixa)) Cmixa <- Cdir(mixa, log = TRUE)

  group <- rep(seq_len(m), each = S)
  mixa0 <- rowSums(mixa)
  out <- matrix(NA, m, m)

  kernel_a0 <- 1 / h + c * D

  for (j in seq_len(m)) {
    v <- mixa[j, ] - 1
    w <- sweep(u, 2, v, "+")
    r1 <- u / w
    r2 <- sweep(1 / w, 2, v, "*")

    Cs1 <- rowSums(r1 * r2 * w / (w - 1))
    Cs2 <- rowSums(r1 * r2)
    Cs <- (1 - 1 / D) * Cs1 + Cs2 / D - rowSums(r1) * rowSums(r2) / D
    logCd <- Ca + Cmixa[j] - Cdir(w + 1, log = TRUE)

    p <- kernel_a0 + mixa0[j] - D - 1
    value <- exp(logCd + log(p) + log(p - 1)) * Cs

    out[, j] <- rowsum(value, group, reorder = FALSE)[, 1] / S
  }
  out
}

DPsi2 <- function(h, mix, c) {
  size <- attr(attr(mix, "samp"), "size")
  y <- attr(h, "y")
  m <- size[1]
  S <- size[2]
  D <- ncol(mix$alpha)
  s <- seq_len(S)

  out <- matrix(NA, m, m)
  kernel_a0 <- 1 / h + c * D
  p <- 2 * kernel_a0 - D - 1
  log_prefactor <- log(p) + log(p - 1)

  for (i in seq_len(m)) {
    ii <- (i - 1) * S + s

    a1 <- y[[1]]$a[ii, , drop = FALSE]
    u1 <- a1 - 1
    Ca1 <- y[[1]]$Ca[ii]

    for (j in seq_len(m)) {
      jj <- (j - 1) * S + s
      a2 <- y[[2]]$a[jj, , drop = FALSE]
      u2 <- a2 - 1
      w <- u1 + u2
      r1 <- u1 / w
      r2 <- u2 / w

      Cs1 <- rowSums(r1 * u2 / (w - 1))
      Cs2 <- rowSums(r1 * r2)
      Cs <- (1 - 1 / D) * Cs1 + Cs2 / D - rowSums(r1) * rowSums(r2) / D
      logCd <- Ca1 + y[[2]]$Ca[jj] - Cdir(w + 1, log = TRUE)

      out[i, j] <- mean(exp(logCd + log_prefactor) * Cs)
    }
  }
  out
}

DMISEf <- function(h, mix, n, c = 1.5, verbose = FALSE) {
  samp <- attr(mix, "samp")
  if (is.null(samp)) stop("No samples generated.")
  attr(h, "y") <- lapply(samp, function(xi) {
    ya <- xi / h + c
    Cya <- Cdir(ya, log = TRUE)
    list(a = ya, Ca = Cya)
  })

  cdh <- DCdh(h, mix, c)
  psi0 <- DPsi0(mix)
  psi1 <- DPsi1(h, mix, c)
  psi2 <- DPsi2(h, mix, c)

  re <- c(cdh / n + t(mix$pr) %*% ((1 - 1/n) * psi2 - 2 * psi1 + psi0) %*% mix$pr)
  if (verbose) print(paste0(round(h, 4), ": ", round(re, 2)))
  re
}
