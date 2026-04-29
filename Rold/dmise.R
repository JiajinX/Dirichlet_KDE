prefactor <- function(x, y, D) {
  t <- outer(x, y, "+") - D - 1
  log(t) + log(t - 1)
}

D1S <- function(x, y) {
  nx <- nrow(x)
  ny <- nrow(y)
  rx <- (x - 1)[rep(seq_len(nx), ny), , drop = FALSE]
  ry <- (y - 1)[rep(seq_len(ny), each = nx), , drop = FALSE]
  up <- log(rx) + log(ry)
  rxy <- rx + ry
  rxy_1 <- rxy - 1
  if (any(rxy_1 < 0)) stop("x or y is too small.")
  down <- log(rxy) + log(rxy_1)
  r <- rowSums(exp(up - down))
  dim(r) <- c(nx, ny)
  r
}

DCdh <- function(h, mix, c) {
  size <- attr(attr(mix, "samp"), "size")
  y <- attr(h, "y")
  if (is.list(y)) y <- y[[1]]
  D <- ncol(mix$alpha)

  C2y <- Cdir(2 * y$a - 1, log = TRUE)
  Cd <- exp(2 * y$Ca - C2y)
  Cs <- rowSums(exp(2 * log(y$a - 1) - log(2 * y$a - 2) - log(2 * y$a - 3)))
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
  exp(log(Cd) + log(Cs) + logC)
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
  logC <- prefactor(1/h + c * D, a0, D)
  exp(log(r) + logC[rep(1, size[1]), ])
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

  Csup <- log(rya1 - 1) + log(rya2 - 1)
  Csdown <- log(sumya - 2) + log(sumya - 3)
  Cs <- rowSums(exp(Csup - Csdown))

  dim(Cd) <- dim(Cs) <- c(prod(size), size[1])
  group <- rep(seq_len(size[1]), each = size[2])
  r <- rowsum(Cd * Cs, group) / size[2]

  logC <- prefactor(1/h + c * D, 1/h + c * D, D)
  exp(log(r) + c(logC))
}

DMISEf <- function(h, mix, n, c = 1) {
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

  c(cdh / n + t(mix$pr) %*% ((1 - 1/n) * psi2 - 2 * psi1 + psi0) %*% mix$pr)
}
