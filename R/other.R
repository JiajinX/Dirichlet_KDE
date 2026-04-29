# Check if points out of simplex boundary after updating and do back tracking
stepup <- function(old, update, bound = 0) {
  if (is.vector(old)) old <- t(old)
  if (is.vector(update)) update <- t(update)
  new <- old + update
  j <- which(rowSums(new < bound) > 0) # Find components with new < boundary
  if (!identical(j, integer(0))) { # Pull new back to boundary
    sold <- old[j, , drop = FALSE]
    supdate <- update[j, , drop = FALSE]
    lf <- (sold - bound) / -supdate
    # Case of update alpha:
    # For a_k = 1, update_k will be >= 0:
    # when a_k = 1 & update_k = 0, lf = NAN
    # when a_k = 1 & update_k > 0, lf = 0
    lf[lf < 0 | is.na(lf)] <- Inf
    # lf < 0: a_k > 1 & update_k > 0
    # lf = NA: update_k = 0
    rmin <- apply(lf, 1, function(r) which.min(r))
    minloc <- cbind(seq_along(j), j, rmin)
    new[j, ] <- sold + lf[minloc[, -2, drop = FALSE]] * supdate
    # Back tracking, pull back to boundary
    new[minloc[, -1, drop = FALSE]] <- bound
  }
  new
}

LCVf <- function(h, x) {
  douter <- ddir(x, x/h + 1, outer = TRUE)
  diag(douter) <- 0
  d <- rowSums(douter)
  sum(log(d[d > 0]))
}

# Find h by max pseudo-log-likelihood (LCV)
optLCVh <- function(x, interval = c(1, 1e-3), tol = 1e-3) {
  opt <- optimize(LCVf, x = x, maximum = TRUE, interval = interval, tol = tol)
  list(type = "Density", h = opt$maximum, LCV = opt$objective)
}

# Generate data from Dirichlet/Normal mixture
rmix <- function(n, mix) {
  fun <- if (inherits(mix, "normix")) rnormix else rdirmix
  fun(n, mix)
}

# Density on the simplex for Dirichlet/Normal mixture/kernel
dspx <- function(x, obj, log = FALSE) {
  fun <- if (inherits(obj, "normix")) dSnormix else ddirpred
  fun(x, obj, log)
}

# Partial derivative of Dirichlet/Normal mixture/kernel
pdspx <- function(x, obj, proj = TRUE) {
  fun <- if (inherits(obj, "normix")) pdnormix else pddir
  fun(x, obj, proj)
}

plotdg <- function(obj, prop = TRUE, scale = NULL, axld = 100, axlg = 30) {
  plot(obj, prop = prop, axl = axld)
  vfplot(obj, add = TRUE, scale = scale, axl = axlg)
}
