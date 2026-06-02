Dir_mix <- function(no, D = 3, c = 1) {
  if (D < 3)  stop("Dimension setting error.")
  p <- switch(no,
    # Same bandwidth
    # 1.
    list(mu = c(0.2, 0.3, 0.5), pr = 1, h = 0.1), # Centered
    # 2.
    list(mu = rbind(c(0.5, 0.3, 0.2), c(0.2, 0.3, 0.5)),
         pr = rep(1, 2) / 2, h = 0.05), # Bimodal
    # 3.
    {
      i <- 1:6
      x1 <- i / 6
      x1 <- x1[1:4]
      r <- 1 - x1
      mu <- cbind(x1, r * x1, r^2)
      colnames(mu) <- NULL
      list(mu = mu, pr = rep(1, 4) / 4, h = 0.08)
    },

    # Different bandwidth
    # 4.
    list(mu = rbind(c(1/3, 1/3 + 0.08, 1/3 - 0.08),
                    c(1/3, 0.1, 0.9 - 1/3),
                    c(1/3 - 0.2, 1/3 + 0.18, 1/3 + 0.02),
                    c(0.44, 0.42, 0.14)),
         pr = c(3, 7, 4, 1) / 15, h = c(0.1, 0.05, 0.035, 0.04)),
    # 5.
    {
      i <- 1:3
      l <- length(i)
      x1 <- i / l * 0.6
      xr <- 1 - x1
      xmu <- cbind(x1, xr * 2/9, xr * 7/9)
      y1 <- i / l * 0.35
      yr <- 1 - y1
      ymu <- cbind(y1, yr * 8/10, yr * 2/10)
      colnames(xmu) <- colnames(ymu) <- NULL
      list(mu = rbind(xmu, ymu), pr = rep(1, 2 * l) / (2 * l),
           h = c(0.1 * rev(x1), 0.2 * rev(y1)))
    },
    # 6.
    {
      i <- 1:6
      l <- length(i)
      x1 <- i / (l + 1.5)
      x1 <- x1[1:5]
      xr <- 1 - x1
      xmu <- cbind(x1, xr * x1, xr^2)
      ymu <- cbind(xmu[, 2] + (1 / 7.5)^2,
                   xmu[, 1] - (1 / 7.5)^2,
                   xmu[, 3])
      colnames(xmu) <- colnames(ymu) <- NULL
      list(mu = rbind(xmu, ymu),
           pr = c(rep(0.45 / l, l), rep(0.55 / l, l)),
           h = c(sqrt(x1) * 0.09, sqrt(rev(x1)) * 0.1))
    },
  )
  if (D > 3) {
    if (is.vector(p$mu)) p$mu <- t(p$mu)
    p$mu <- p$mu * 3 / D
    newcol <- matrix((1 - 3 / D) / (D - 3), nrow = nrow(p$mu), ncol = (D - 3))
    p$mu <- cbind(p$mu, newcol)
  }
  dirmix(a = p$mu / p$h + c, pr = p$pr)
}

# Create multivariate normal distribution from Chacon'S parameters
# D: dimension of the mixture
Norm_mix <- function(no, D = 2) {
  if (D < 2) stop("Dimension setting error.")
  mix <- switch(no,
    # 1.
    list(mu = rbind(c(0, 0)),
         Sigma = ks::invvech(c(1/4, 0, 1)),
         pr = 1),
    # 2.
    list(mu = rbind(c(1, 0), c(-1, 0)) * 3/4,
         Sigma = rbind(diag(4/9, 2), diag(4/9, 2)),
         pr = c(1, 1) / 2),
    # 3.
    list(mu = rbind(c(1, -1), c(-1, 1)) * 2/5,
         Sigma = rbind(ks::invvech(c(4/9, 14/45, 4/9)), diag(4/9, 2)),
         pr = c(1, 1) / 2),
    # 4.
    list(mu = rbind(c(-3/5, 0), c(1/2, 3/5), c(1/2, -3/5)),
         Sigma = rbind(ks::invvech(c(9, 63 / 10, 49 / 4)),
                       ks::invvech(c(9, 0, 49 / 4)),
                       ks::invvech(c(9, 0, 49 / 4))) / 25,
         pr = c(3, 3, 1) / 7),
    # 5.
    list(mu = rbind(c(-1, 0), c(1, 0), c(0, 0)),
         Sigma = matrix(c(rep(ks::invvech(c(4/9, 4/15, 4/9)), 2),
                          ks::invvech(c(1, 3/5, 1)) / 9),
                        byrow = TRUE, nrow = 3 * 2),
         pr = c(11/25, 11/25, 3/25)),
    # 6.
    list(mu = matrix(c(rep(0, 2),
                       c(-1, -1, 1, -1, -1, 1, 1, 1) * 4/5),
                     byrow = TRUE, nrow = 5),
         Sigma = matrix(rep(diag(1 / 10, 2), 5),
                        byrow = TRUE, nrow = 5 * 2),
         pr = rep(1/5, 5))
  )
  if (D > 2) {
    m <- length(mix$pr)
    newcol <- matrix(0, nrow = m, ncol = D - 2)
    mix$mu <- cbind(mix$mu, newcol) # Add 0 for mu of each component
    spSig <- split(mix$Sigma, rep(1:m, each = 2))
    exspSig <- lapply(spSig, function(x) {
      dim(x) <- c(2, 2)
      eigenvalue <- eigen(x)$values
      neweig <- seq(min(eigenvalue), max(eigenvalue), len = D)[-c(1, D)]
      # Add eigen value by sequence
      exmtrix <- matrix(0, nrow = D, ncol = D)
      exmtrix[1:2, 1:2] <- x # Original Sigma
      diag(exmtrix)[3:D] <- neweig
      exmtrix
    })
    mix$Sigma <- do.call(rbind, exspSig)
  }
  class(mix) <- "normix"
  mix
}

plot_simu <- function(type = c("dirichlet", "normal"), prop = TRUE, no = 1:6) {
  type <- match.arg(type)
  fun <- if (type == "dirichlet") Dir_mix else Norm_mix
  oldpar <- par(mfrow = c(2, 3))
  for (i in no) plot(fun(i), prop = prop)
  par(oldpar)
}
