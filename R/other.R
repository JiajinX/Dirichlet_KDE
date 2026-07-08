# Stably calculate log(sum(exp(x))).
log_sum_exp <- function(x, na.rm = FALSE) {
  if (na.rm) x <- x[!is.na(x)]
  if (length(x) == 0) return(-Inf)
  mx <- max(x)
  if (is.na(mx) || !is.finite(mx)) return(mx)
  mx + log(sum(exp(x - mx)))
}

# Stably calculate sum(sgn * exp(x)) as sign and log absolute value.
signed_log_sum_exp <- function(x, sgn) {
  if (length(x) != length(sgn)) stop("x and sgn must have the same length.")
  keep <- sgn != 0 & x > -Inf
  x <- x[keep]
  sgn <- sgn[keep]
  if (length(x) == 0) return(c(sign = 0, logabs = -Inf))
  if (any(is.na(x)) || any(is.na(sgn))) return(c(sign = NA, logabs = NA))
  if (any(!is.finite(sgn))) stop("sgn must be finite.")

  lpos <- log_sum_exp(x[sgn > 0])
  lneg <- log_sum_exp(x[sgn < 0])
  if (lpos == lneg) return(c(sign = 0, logabs = -Inf))
  if (lpos > lneg) return(c(sign = 1, logabs = lpos + log1p(-exp(lneg - lpos))))
  c(sign = -1, logabs = lneg + log1p(-exp(lpos - lneg)))
}

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
  if (!is.matrix(x)) x <- as.matrix(x)
  opt <- optimize(LCVf, x = x, maximum = TRUE, interval = interval, tol = tol)
  list(type = "Density", h = opt$maximum, LCV = opt$objective)
}

# Generate data from Dirichlet/Normal mixture
rmix <- function(n, mix, seed = NULL) {
  set.seed(seed)
  fun <- if (inherits(mix, "normix")) rnormix else rdirmix
  fun(n, mix)
}

# Density on the simplex for Dirichlet/Normal mixture/kernel
dspx <- function(x, obj, log = FALSE) {
  if (is.data.frame(x)) x <- as.matrix(x)
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

# Set color transparency
transcols <- function(cols, alpha) {
  paste0(cols, as.hexmode(round(alpha * 255)))
}

closure <- function(x) {
  if (is.vector(x)) x <- t(x)
  if (any(x < 0)) stop("Each value should >= 0.")
  x / rowSums(x)
}

co_aggregate <- function(x, comb, closure = TRUE) {
  if (closure) x <- closure(x)
  cols <- sort(unique(do.call(c, comb)))
  if (any(seq_len(ncol(x)) != cols)) stop("Problem with the comb settings.")
  sapply(comb, function(cols) rowSums(x[, cols, drop = FALSE]))
}

dir_aggregate <- function(obj, comb) {
  if (class(obj) == "dirmix") {
    dirmix(co_aggregate(obj$alpha, comb, FALSE), obj$pr)
  } else dirkn(co_aggregate(obj$x, comb), obj$h, obj$c)
}
