# Find the optimal h by grid search to minimize the MISE or DMISE(r=1)
MISEfh <- function(vh, mix, n, r = 0, c = NULL,
                   exhaustive = FALSE, verbose = FALSE) {
  if (r > 1) stop("Only density and its first derivative are being considered.")
  if (r == 1 && is.null(c)) c <- 1.5
  if (is.null(attr(mix, "samp"))) mix <- gene_compsamp(mix)
  info <- rep(NA, length(vh))
  names(info) <- round(vh, 4)
  rold <- Inf
  opt <- FALSE
  for (i in seq_along(vh)) {
    if (verbose) cat(names(info)[i], " ")
    info[i] <- if (r == 0) MISEf(vh[i], mix, n) else DMISEf(vh[i], mix, n, c)
    if (info[i] > rold && opt == FALSE) {
      opt <- TRUE
      imin <- i - 1
      if (!exhaustive) break
    }
    rold <- info[i]
  }
  if (verbose) cat("\n")
  if (!opt) imin <- i
  if (imin == 1) opt <- FALSE
  info <- info[1:i]
  attr(info, "h") <- vh[1:i]
  tp <- ifelse(r == 0, "Density", "First derivative")
  list(type = tp, info = info, opt = opt, h = vh[imin],
       MISE = as.numeric(info[imin]))
}

# Find the optimal h by grid search to minimize LSCV or DLSCV(r=1)
LSCVfh <- function(vh, x, r = 0, c = NULL,
                   exhaustive = FALSE, verbose = FALSE) {
  if (r > 1) stop("Only density and its first derivative are being considered.")
  if (r == 1 && is.null(c)) c <- 1.5
  info <- rep(NA, length(vh))
  names(info) <- round(vh, 4)
  rold <- Inf
  opt <- FALSE
  for (i in seq_along(vh)) {
    if (verbose) cat(names(info)[i], " ")
    info[i] <- if (r == 0) LSCVf(vh[i], x) else DLSCVf(vh[i], x, c)
    if (info[i] > rold && opt == FALSE) {
      opt <- TRUE
      imin <- i - 1
      if (!exhaustive) break
    }
    rold <- info[i]
  }
  if (verbose) cat("\n")
  if (!opt) imin <- i
  info <- info[1:i]
  attr(info, "h") <- vh[1:i]
  tp <- ifelse(r == 0, "Density", "First derivative")
  list(type = tp, info = info, opt = opt, h = vh[imin],
       LSCV = as.numeric(info[imin]))
}

# Find the optimal h by optimize function to minimize the MISE or DMISE(r=1)
optMISEh <- function(mix, n, r = 0, c = NULL,
                     interval = c(1, 1e-3), tol = 1e-3, verbose = FALSE) {
  if (r > 1) stop("Only density and its first derivative are being considered.")
  if (is.null(attr(mix, "samp"))) mix <- gene_compsamp(mix)
  if (r == 0) {
    tp <- "Density"
    opt <- optimize(MISEf, mix = mix, n = n,
                    interval = interval, tol = tol, verbose = verbose)
  } else if (r == 1) {
    tp <- "First derivative"
    if (is.null(c)) c <- 1.5
    opt <- optimize(DMISEf, mix = mix, n = n, c = c,
                    interval = interval, tol = tol, verbose = verbose)
  }
  list(type = tp, h = opt$minimum, MISE = opt$objective)
}

# Find the optimal h by optimize function to minimize the MISE or DMISE(r=1)
optLSCVh <- function(x, r = 0, c = NULL,
                     interval = c(1, 1e-3), tol = 1e-3, verbose = FALSE) {
  if (r > 1) stop("Only density and its first derivative are being considered.")
  if (!is.matrix(x)) x <- as.matrix(x)
  if (r == 0) {
    tp <- "Density"
    opt <- optimize(LSCVf, x = x,
                    interval = interval, tol = tol, verbose = verbose)
  } else if (r == 1) {
    tp <- "First derivative"
    if (is.null(c)) c <- 1.5
    opt <- optimize(DLSCVf, x = x, c = c,
                    interval = interval, tol = tol, verbose = verbose)
  }
  list(type = tp, h = opt$minimum, LSCV = opt$objective)
}
