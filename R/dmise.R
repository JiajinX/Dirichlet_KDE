prefactor <- function(x, y, D) {
  t <- outer(x, y, "+") - D - 1
  log(t) + log(t - 1)
}

D1Scut <- function(logx, logy, logxy) {
  D <- ncol(logx)
  tx <- rowSums(exp(logx - logxy))
  ty <- rowSums(exp(logy - logxy))
  tx * ty / D
}

D1S <- function(x, y) {
  D <- ncol(x)
  nx <- nrow(x)
  ny <- nrow(y)
  rx <- (x - 1)[rep(seq_len(nx), ny), , drop = FALSE]
  ry <- (y - 1)[rep(seq_len(ny), each = nx), , drop = FALSE]
  log_rx <- log(rx)
  log_ry <- log(ry)
  up <- log_rx + log_ry
  rxy <- rx + ry
  rxy_1 <- rxy - 1
  if (any(rxy_1 <= 0)) stop("x or y is too small.")
  log_rxy <- log(rxy)
  log_rxy_1 <- log(rxy_1)
  down1 <- log_rxy + log_rxy_1
  down2 <- 2 * log_rxy
  s1 <- rowSums(exp(up - down1))
  s2 <- rowSums(exp(up - down2))
  scut <- D1Scut(log_rx, log_ry, log_rxy)
  r <- (1 - 1/D) * s1 + s2/D - scut
  dim(r) <- c(nx, ny)
  r
}

DCdh <- function(h, mix, c) {
  size <- attr(attr(mix, "samp"), "size")
  y <- attr(h, "y")
  if (is.list(y)) y <- y[[1]]
  D <- ncol(mix$alpha)

  ya2 <- 2 * y$a
  C2y <- Cdir(ya2 - 1, log = TRUE)
  Cd <- exp(2 * y$Ca - C2y)
  Cs <- (1 - 1/D) * # nolint: infix_spaces_linter. # nolint: infix_spaces_linter.
    rowSums(exp(log(y$a - 1) - log(2) - log(ya2 - 3))) + (1 - D)/4 # nolint: infix_spaces_linter.
  r <- Cd * Cs
  dim(r) <- rev(size)

  logc <- prefactor(1/h + c * D, 1/h + c * D, D) # nolint: infix_spaces_linter.
  sum(colMeans(r) * mix$pr) * exp(c(logc))
}

DPsi0 <- function(mix) {
  Ca <- attr(mix$alpha, "logmbeta")
  outCaa <- outer(Ca, Ca, "+")
  Caa <- Cdir_plus(mix$alpha, mix$alpha, log = TRUE)
  Cd <- exp(outCaa - Caa)
  Cs <- D1S(mix$alpha, mix$alpha)
  a0 <- rowSums(mix$alpha)
  D <- ncol(mix$alpha)
  logC <- prefactor(a0, a0, D = D)
  exp(log(Cd) + logC) * Cs
}

DPsi1 <- function(h, mix, c) {
  size <- attr(attr(mix, "samp"), "size")
  Ca <- attr(mix$alpha, "logmbeta")
  y <- attr(h, "y")
  if (is.list(y)) y <- y[[1]]
  a0 <- rowSums(mix$alpha)
  D <- ncol(mix$alpha)

  outCay <- outer(y$Ca, Ca, "+")
  Cay <- Cdir_plus(y$a, mix$alpha, log = TRUE)
  Cd <- exp(outCay - Cay)
  Cs <- D1S(y$a, mix$alpha)
  group <- rep(seq_len(size[1]), each = size[2])
  r <- rowsum(Cd * Cs, group) / size[2]
  logC <- prefactor(1/h + c * D, a0, D) # nolint: infix_spaces_linter.
  r * exp(logC[rep(1, size[1]), ])
}

DPsi2 <- function(h, mix, c) {
  size <- attr(attr(mix, "samp"), "size")
  y <- attr(h, "y")
  D <- ncol(mix$alpha)

  v <- seq_len(prod(size))
  i <- rep(v, size[1])
  mtx <- matrix(v, nrow = size[2])
  j <- c(mtx[rep(seq_len(size[2]), size[1]), ])

  up <- y[[1]]$Ca[i] + y[[2]]$Ca[j]
  rya1 <- y[[1]]$a[i, ]
  rya2 <- y[[2]]$a[j, ]
  sumya <- rya1 + rya2
  down <- Cdir(sumya - 1, log = TRUE)
  Cd <- exp(up - down)

  log_rya1_1 <- log(rya1 - 1)
  log_rya2_1 <- log(rya2 - 1)
  Csup <- log_rya1_1 + log_rya2_1
  log_sumya_2 <- log(sumya - 2)
  log_sumya_3 <- log(sumya - 3)
  Csdown1 <- log_sumya_2 + log_sumya_3
  Csdown2 <- 2 * log_sumya_2
  Cs1 <- rowSums(exp(Csup - Csdown1))
  Cs2 <- rowSums(exp(Csup - Csdown2))
  Cscut <- D1Scut(log_rya1_1, log_rya2_1, log_sumya_2)
  Cs <- (1 - 1/D) * Cs1 + Cs2/D - Cscut # nolint: infix_spaces_linter. # nolint: infix_spaces_linter.

  dim(Cd) <- dim(Cs) <- c(prod(size), size[1])
  group <- rep(seq_len(size[1]), each = size[2])
  r <- rowsum(Cd * Cs, group) / size[2]

  logC <- prefactor(1/h + c * D, 1/h + c * D, D) # nolint: infix_spaces_linter. # nolint: infix_spaces_linter.
  r * exp(c(logC))
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

  re <- c(cdh / n + t(mix$pr) %*% ((1 - 1/n) * psi2 - 2 * psi1 + psi0) %*% mix$pr) # nolint: infix_spaces_linter.
  if (verbose) print(paste0(round(h, 4), ": ", round(re, 2)))
  re
}
