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
  logrx <- log(rx)
  logry <- log(ry)
  up <- logrx + logry
  rxy <- rx + ry
  rxy_1 <- rxy - 1
  if (any(rxy_1 < 0)) stop("x or y is too small.")
  logrxy <- log(rxy)
  logrxy_1 <- log(rxy_1)
  down1 <- logrxy + logrxy_1
  down2 <- 2 * logrxy
  r <- (1 - 1/D) * rowSums(exp(up - down1)) + 1/D * rowSums(exp(up - down2))
  cut <- D1Scut(logrx, logry, logrxy)
  r <- r - cut

  dim(r) <- c(nx, ny)
  r
}

DCdh <- function(h, mix, c) {
  samp <- attr(mix, "samp")
  size <- attr(samp, "size")
  D <- ncol(mix$alpha)
  Cy <- attr(h, "Cy")
  if (is.list(Cy)) {
    Cy <- Cy[[1]]
    samp <- samp[[1]]
  }

  sh <- samp / h + c
  C2y <- Cdir(2 * sh - 1, log = TRUE)
  Cd <- exp(2 * Cy - C2y)
  Cs <- (1 - 1/D) *
    rowSums(exp(2 * log(sh - 1) - log(2 * sh - 2) - log(2 * sh - 3))) +
    (1 - D)/4
  r <- Cd * Cs
  dim(r) <- rev(size)

  logc <- prefactor(1/h + c * D, 1/h + c * D, D)
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
  Cd * Cs * exp(logC)
  # exp(log(Cd) + log(Cs) + logC)
}

DPsi1 <- function(h, mix, c) {
  samp <- attr(mix, "samp")
  size <- attr(samp, "size")
  Ca <- attr(mix$alpha, "logmbeta")
  Cy <- attr(h, "Cy")
  if (is.list(Cy)) {
    Cy <- Cy[[1]]
    samp <- samp[[1]]
  }

  outCay <- outer(Cy, Ca, "+")
  sh <- samp / h + c
  Cya <- Cdir_plus(sh, mix$alpha, log = TRUE)
  Cd <- exp(outCay - Cya)
  Cs <- D1S(sh, mix$alpha)
  K <- kronecker(diag(size[1]), matrix(1/size[2], 1, size[2]))
  r <- K %*% (Cd * Cs)

  a0 <- rowSums(mix$alpha)
  D <- ncol(mix$alpha)
  logC <- prefactor(1/h + c * D, a0, D)
  r * exp(c(logC))
}

DPsi2 <- function(h, mix, c) {
  samp <- attr(mix, "samp")
  size <- attr(samp, "size")
  Cy <- attr(h, "Cy")
  D <- ncol(mix$alpha)

  v <- seq_len(prod(size))
  i <- rep(v, size[1])
  mtx <- matrix(v, nrow = size[2])
  j <- c(mtx[rep(seq_len(size[2]), size[1]), ])

  up <- Cy[[1]][i] + Cy[[2]][j]
  sh1 <- samp[[1]][i, ] / h + c
  sh2 <- samp[[2]][j, ] / h + c
  sumsha <- sh1 + sh2
  down <- Cdir(sumsha - 1, log = TRUE)
  Cd <- exp(up - down)

  logsh1_1 <- log(sh1 - 1)
  logsh2_1 <- log(sh2 - 1)
  Csup <- logsh1_1 + logsh2_1
  logsumsh_2 <- log(sumsha - 2)
  logsumsh_3 <- log(sumsha - 3)
  Csdown1 <- logsumsh_2 + logsumsh_3
  Csdown2 <- 2 * logsumsh_2
  Cs <- (1 - 1/D) * rowSums(exp(Csup - Csdown1)) +
    rowSums(exp(Csup - Csdown2)) / D
  cut <- D1Scut(logsh1_1, logsh2_1, logsumsh_2)
  Cs <- Cs - cut

  dim(Cd) <- dim(Cs) <- c(prod(size), size[1])
  K <- kronecker(diag(size[1]), matrix(1/size[2], 1, size[2]))
  r <- K %*% (Cd * Cs)

  logC <- prefactor(1/h + c * D, 1/h + c * D, D)
  r * exp(c(logC))
}

DMISEf <- function(h, mix, n, c = 1) {
  samp <- attr(mix, "samp")
  if (is.null(samp)) stop("No samples generated.")
  attr(h, "Cy") <- lapply(samp, function(xi) Cdir(xi / h + c, log = TRUE))

  cdh <- DCdh(h, mix, c)
  psi0 <- DPsi0(mix)
  psi1 <- DPsi1(h, mix, c)
  psi2 <- DPsi2(h, mix, c)

  c(cdh / n + t(mix$pr) %*% ((1 - 1/n) * psi2 - 2 * psi1 + psi0) %*% mix$pr)
}
