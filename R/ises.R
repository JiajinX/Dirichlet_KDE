# Got (D)MISE by calculating mean value of (D)ISEs

# Generate data to do kernel estimation
gene_data <- function(n, mix, B, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x <- replicate(B, rmix(n, mix), simplify = FALSE)
  if (inherits(mix, "normix")) x <- lapply(x, function(rx) ilrts(rx)$sx)
  attr(mix, "data") <- x
  mix
}

# Generate evaluate points for Monte Carlo method
gene_eval <- function(mix, r = 0, neva = 1e3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  D <- ifelse(inherits(mix, "dirmix"), ncol(mix$alpha), ncol(mix$mu) + 1)
  eva <- gtools::rdirichlet(neva, rep(1, D))
  deva <- if (r == 0) dspx(eva, mix) else pdspx(eva, mix)
  attr(mix, "evapoint") <- list(eva = eva, deva = deva)
  mix
}

# Calculate MISE or DMISE by Monte Carlo method - mean(ISEs/DISEs)
MISEmc <- function(h, mix, r = 0, c = NULL, verbose = FALSE) {
  if (r == 1 && is.null(c)) c <- 1.5
  D <- ifelse(inherits(mix, "dirmix"), ncol(mix$alpha), ncol(mix$mu) + 1)
  evapoint <- attr(mix, "evapoint")
  data <- attr(mix, "data")
  if (is.null(evapoint) || is.null(data))
    stop("No sample generated from the mixture.")

  eva <- evapoint$eva
  deva <- evapoint$deva
  B <- length(data)
  ISEs <- numeric(B)
  for (i in seq_len(B)) {
    if (verbose) cat(i, " ")
    kn <- dirkn(data[[i]], h, c)
    if (r == 0) {
      df <- ddirpred(eva, kn)
      ISEs[i] <- mean((df - deva)^2) / gamma(D)
    } else if (r == 1) {
      df <- pddir(eva, kn)
      ISEs[i] <- mean(rowSums((df - deva)^2)) / gamma(D)
    }
  }
  if (verbose) cat("\n")
  MISE <- mean(ISEs)
  attr(MISE, "ISE") <- ISEs
  MISE
}

# Generate deterministic grid on 2D-simplex
grid_spx <- function(split = 32) {
  n <- split^2
  ll <- split * 3
  area <- 1 / 2 / n
  up <- seq(1, ll, by = 3)
  down <- up[-length(up)] + 1
  upcomb <- expand.grid(up, up)
  upcomb[, 3] <- ll - rowSums(upcomb)
  downcomb <- expand.grid(down, down)
  downcomb[, 3] <- ll - rowSums(downcomb)
  comb <- rbind(upcomb, downcomb)
  comb <- comb[comb[, 3] > 0, ]
  comb <- as.matrix(comb)
  rownames(comb) <- colnames(comb) <- NULL
  x <- comb / ll
  list(x = x, area  = area)
}

# Generate grids
gene_grid <- function(mix, r = 0, split = 32) {
  # Generate grid and corresponding density
  grid <- grid_spx(split)
  dgrid <- if (r == 0) ddirpred(grid$x, mix) else pddir(grid$x, mix)
  attr(mix, "gridinfo") <- list(x = grid$x, dx = dgrid, area = grid$area)
  mix
}

# Calculate MISE or DMISE by determined grids - mean(ISEs/DISEs)
# d = (f1 - f2)^2 of center of sub triangles; area = area of sub triangles.
# ISE = sum(d * area)
MISEdt <- function(h, mix, r = 0, c = NULL, verbose = FALSE) {
  if (r == 1 && is.null(c)) c <- 1.5
  data <- attr(mix, "data")
  grid <- attr(mix, "gridinfo")
  if (is.null(grid) || is.null(data)) stop("No data or grid provided.")

  gd <- grid$x
  dmix <- grid$dx
  B <- length(data)
  ISEs <- numeric(B)
  for (i in seq_len(B)) {
    if (verbose) cat(i, " ")
    kn <- dirkn(data[[i]], h, c)
    df <- if (r == 0) ddirpred(gd, kn) else pddir(gd, kn)
    ISEs[i] <- sum((df - dmix)^2) * grid$area
  }
  if (verbose) cat("\n")
  MISE <- mean(ISEs)
  attr(MISE, "ISE") <- ISEs
  MISE
}

# Calculate ISEs and MISE of a set of bandwidth h
# fun: MISEmc, MISEdt
MISEh <- function(vh, mix, fun, r = 0, c = NULL,
                  exhaustive = FALSE, verbose = FALSE) {
  if (r > 1) stop("Only density and its first derivative are being considered.")
  if (r == 1 && is.null(c)) c <- 1.5
  lh <- length(vh)
  info <- rep(NA, lh)
  ise <- matrix(nrow = length(attr(mix, "data")), ncol = lh)
  colnames(ise) <- names(info) <- round(vh, 4)
  rold <- Inf
  opt <- FALSE
  for (i in seq_along(vh)) {
    if (verbose) cat(names(info)[i], " ")
    mu <- fun(vh[i], mix, r, c)
    info[i] <- mu
    ise[, i] <- attr(mu, "ISE")
    if (info[i] > rold && opt == FALSE) {
      opt <- TRUE
      imin <- i - 1
      if (!exhaustive) break
    }
    rold <- info[i]
  }
  if (verbose) cat("\n")
  if (!opt) imin <- i # not find opt until end
  if (imin == 1) opt <- FALSE # first h is opt
  info <- info[1:i]
  attr(info, "h") <- vh[1:i]
  ise <- ise[, 1:i, drop = FALSE]
  tp <- ifelse(r == 0, "Density", "First derivative")
  re <- list(type = tp, info = info, opt = opt, h = vh[imin],
             MISE = as.numeric(info[imin]))
  attr(re, "ISE") <- ise
  re
}

# Plot lines of bandwidth h vs MISE
plot_mise <- function(mise, col = 1, ylim = NULL, type = "o", add = FALSE) {
  h <- attr(mise, "h")
  if (!add) {
    plot(h, mise, xlim = rev(range(h)), ylim = ylim, col = col, type = type)
  } else lines(h, mise, col = col, type = type)
}

# Test significant difference between MISE value and a set of ISE values
MISEtest <- function(ISE, MISE) {
  sapply(seq_len(length(MISE)), function(i) {
    r <- t.test(ISE[[i]], mu = MISE[i])
    r$p.value
  })
}
