# Log-likelihood of Dirichlet mixture
logLik.dirmix <- function(mix, x) {
  sum(ddirmix(x, mix, log = TRUE))
}

# Fit one Dirichlet distribution for the data
# Modify from compositions::fitDirichlet with restrict of alpha >= 1
fitDir <- function(x, bound = 1, tol = 1e-4, maxit = 50) {
  if (is.vector(x)) x <- t(x)
  n <- nrow(x)
  D <- ncol(x)
  ultx <- compositions::ult(x)
  if (is.null(dim(ultx))) ultx <- t(ultx)
  elog <- colMeans(ultx)
  alpha <- rep(2, D)
  ll <- ll.old <-
    -n * (sum(lgamma(alpha)) - lgamma(sum(alpha)) - sum(elog * (alpha - 1)))
  for (i in 1:maxit) {
    E <- digamma(alpha) - digamma(sum(alpha))
    V <- diag(trigamma(alpha)) - trigamma(sum(alpha))
    delta = sqrt(sum((elog - E)^2))
    update <- solve(V, (elog - E))
    update[update < 0 & alpha == bound] <- 0
    alpha <- stepup(alpha, update, bound = bound)[1, ] # bound = 1.5????
    ll <-
      -n * (sum(lgamma(alpha)) - lgamma(sum(alpha)) - sum(elog * (alpha - 1)))
    if (ll <= ll.old + tol) break
  }
  list(mix = dirmix(alpha), ll = ll)
}

initEMDir <- function(data, m, bound = 1, minsize = 2, try = 50) {
  for (t in seq_len(try)) { # Avoid a cluster has only one point
    km <- kmeans(data, m)
    sz <- tabulate(km$cluster)
    if (all(sz >= minsize)) break
  }
  if (!all(sz >= minsize)) {
    stop("Failed to find a kmeans partition.")
  }
  ia <- t(sapply(seq_len(m), function(i) {
    fitDir(data[km$cluster == i, ,  drop = FALSE], bound = bound)$mix$alpha
  }))
  pr <- tabulate(km$cluster) / length(km$cluster)
  dirmix(a = ia, pr = pr) # What if a cluster has one point?
}

# Delete support points with pr close to 0
dropsp <- function(mix) {
  j <- mix$pr < 1e-10
  if (any(j)) { # Delete support points with pr close to 0
    mix$pr <- mix$pr[!j]
    mix$pr <- mix$pr / sum(mix$pr)
    logmbeta <- attr(mix$alpha, "logmbeta")[!j]
    mix$alpha <- mix$alpha[!j, ]
    attr(mix$alpha, "logmbeta") <- logmbeta
    p <- p[, !j, drop = FALSE]
  }
  mix
}

# Objective function needs to be maximized by iteration alpha (EM -> iter_alpha)
# Also with the first and second derivative of objective function
# x: data; a: alpha; p: value get from EM
calfun_alpha <- function(x, a, p, deriv = 0) {
  n <- nrow(x)
  D <- ncol(x)
  m <- nrow(a)
  logx <- attr(x, "log")
  if (is.null(logx)) logx <- pmax(log(x), log(1e-100))
  if (deriv != 0) dim(p) <- NULL
  xp <-
    switch(deriv + 1,
           rowSums((tcrossprod(a, logx) - rowSums(lgamma(a))) * t(p)),
           # Original objective function
           (logx[rep(1:n, m), ] - digamma(a)[rep(1:m, each = n), ]) * p,
           # p_ij * [log(x_ik) - ψ(α_jk)] - first derivative
           -trigamma(a)[rep(1:m, each = n), ] * p
           # -p_ij * ψ'(α_jk) - second derivative
    )
  if (deriv != 0) {
    dim(xp) <- c(n, m * D)
    xp <- colSums(xp)
    dim(xp) <- c(m, D)
  }
  a0 <- rowSums(a)
  dim(p) <- c(n, m)
  at <- switch(deriv + 1, # Add a0 term when a0 is not fixed
               rowSums(t(p) * lgamma(a0)),
               rowSums(t(p) * digamma(a0)),
               rowSums(t(p) * trigamma(a0)))
  xp <- if (deriv == 2) {
    lapply(1:m, function(i) diag(xp[i, ]) + at[i])
  } else xp + at
  xp
}

# Updating alpha by using Newton-Raphson method (EM)
# x: data; a: alpha; p: value get from EM
iter_alpha <- function(x, a, p, bound = 1,
                       maxit = 100, maxbt = 50, tol = 1e-4) {
  m <- nrow(a)
  obj.old <- obj <- calfun_alpha(x, a, p)
  fa.old <- fa <- sum(obj)
  # cat("\nalpha-function:"); print(fa, digits = 16)
  for (i in seq_len(maxit)) {
    g <- calfun_alpha(x, a, p, deriv = 1)
    H <- calfun_alpha(x, a, p, deriv = 2)

    a.old <- a
    a <- t(sapply(1:m, function(comp) {
      P <- H[[comp]] %*% (a[comp, ] - bound) - g[comp, ]
      Q <- -H[[comp]]
      lsei::pnnqp(Q, P)$x + bound # 1 -> bound (set by user)
    }))
    obj <- calfun_alpha(x, a, p)

    # Back Tracking part for heterogeneous finite mixture
    ibt <- which(obj < obj.old) # The True components need back tracking
    if (!identical(ibt, integer(0))) {
      dea <- 0.333 *
        rowSums(g[ibt, , drop = FALSE] *
                  (a[ibt, , drop = FALSE] - a.old[ibt, , drop = FALSE]))
      alpha <- 1
      a.up <- a
      for (k in seq_len(maxbt)) { # back tracking
        obj[ibt] <- calfun_alpha(x, a[ibt, , drop = FALSE],
                                 p[, ibt, drop = FALSE])
        idx <- which(obj[ibt] < obj.old[ibt] + dea * alpha)
        # TRUE for components still need back tracking
        if (identical(idx, integer(0))) break
        ibt <- ibt[idx]
        dea <- dea[idx]
        alpha <- alpha * 0.5
        a[ibt, ] <- (1 - alpha) * a.old[ibt, , drop = FALSE] +
          alpha * a.up[ibt, , drop = FALSE]
      }
    }
    fa <- sum(obj)
    # cat("\nalpha-function:"); print(fa, digits = 16)
    if (fa <= fa.old + tol) break
    obj.old <- obj
    fa.old <- fa
  } # The alpha-function seems to be monotonically increasing for each iteration
  attr(a, "logmbeta") <- lgamma(rowSums(a)) - rowSums(lgamma(a))
  a
}

# EM of heterogeneous Dirichlet mixture
EMDir <- function(data, m, bound = 1, maxit = 100, tol = 1e-4) {
  n <- nrow(data)
  mix <- initEMDir(data, m, bound) # Initial part
  ll <- logLik(mix, data)
  for (i in seq_len(maxit)) {
    logd <- ddir(data, mix$alpha, log = TRUE, outer = TRUE) +
      rep(log(mix$pr), rep(n, length(mix$pr)))
    # log f(x_i;α_j) + log π_j = log [π_j * f(x_i;α_j)]
    ma <- lsei::matMaxs(logd)
    pi.d <- exp(logd - ma)
    pi.s <- rowSums(pi.d)
    p <- pi.d / pi.s
    # -ma or not will not change p, which is the same as
    mix$pr <- colMeans(p) # π_j'
    mix <- dropsp(mix)
    mix$alpha <- iter_alpha(data, mix$alpha, p, bound) # Update alpha

    ll.old <- ll
    ll <- logLik(mix, data)
    if (ll <= ll.old + tol) break
  }
  at <- mix$alpha - 1
  attr(at, "logmbeta") <- NULL
  attr(mix, "mode") <- at / rowSums(at)
  list(mix = mix, ll = ll)
}

# Merge the similar components and drop tiny components (collapse_points)
# prec: threshold of merge and drop
unique.dirmix <- function(mix, prec = 0) {
  if (nrow(mix$alpha) == 1) return(mix)
  prec <- rep(prec, len = 2)
  # prec[1] for difference in each pair of support points, prec[2] for pr
  j <- mix$pr <= prec[2]
  a <- mix$alpha[!j, , drop = FALSE]
  pr <- mix$pr[!j] # Drop components with small pr
  if (prec[1] == 0) return(dirmix(a, pr))
  repeat {
    if (nrow(a) == 1) break
    d <- abs(a[-1, , drop = FALSE] - a[-nrow(a), , drop = FALSE])
    # Calculate difference of each pair of support points
    i <- lsei::matMaxs(d) <= prec[1] # Will use odd indexes of TRUE
    if (sum(i) == 0) break
    noti <- c(!i, TRUE) # not i
    even <- seq(2, by = 2, len = length(i))
    i[even] <- i[even] & noti[even - 1] & noti[even + 1]
    # Use evens if only isolated
    ind <- which(i) # Indexes of points need to be merged
    pr2 <- pr[ind] + pr[ind + 1]
    a[ind, ] <- (a[ind, , drop = FALSE] * pr[ind] +
                   a[ind + 1, , drop = FALSE] * pr[ind + 1]) / pr2
    # Collapse neighboring pairs to the first
    pr[ind] <- pr2
    a <- a[-(ind + 1), , drop = FALSE]
    pr <- pr[-(ind + 1)] # Remove the second of a pair
  }
  dirmix(a, pr)
}

collapse_points <- function(mix, data, tol = 1e-3) {
  ll <- logLik(mix, data)
  j <- mix$pr == 0 # Drop components with mass = 0
  if (any(j)) {
    logmbeta <- attr(mix$alpha, "logmbeta")[!j]
    mix$alpha <- mix$alpha[!j, , drop = FALSE]
    mix$pr <- mix$pr[!j]
    attr(mix$alpha, "logmbeta") <- logmbeta
    attr(mix, "mode") <- attr(mix, "mode")[!j, , drop = FALSE]
    ll <- logLik(mix, data) # Update log-likelihood after drop components
  }
  m <- length(mix$pr)
  if (m == 1 || tol <= 0) return(list(mix = mix, ll = ll))
  o <- do.call(order, split(mix$alpha, col(mix$alpha))) # Order by each column
  logmbeta <- attr(mix$alpha, "logmbeta")[o]
  mix$alpha <- mix$alpha[o, , drop = FALSE]
  mix$pr <- mix$pr[o]
  attr(mix$alpha, "logmbeta") <- logmbeta
  prec <- 1e-10 # Initial precision value
  repeat {
    if (nrow(mix$alpha) <= 1) break
    d <- abs(mix$alpha[-1, , drop = FALSE] -
               mix$alpha[-nrow(mix$alpha), , drop = FALSE])
    prec <- 10 * max(min(lsei::matMaxs(d)), prec)
    mixt <- unique(mix, prec = c(prec, 0))
    # Merge all pairs of components with difference less than
    # the current smallest difference * 10
    llt <- logLik(mixt, data)
    if (llt >= ll - tol) mix <- mixt
    else break
  }
  ll <- logLik(mix, data)
  list(mix = mix, ll = ll)
}

optfm <- function(data, bound = 1, criterion = c("BIC", "AIC", "AICc"), mtol = 3,
                  maxrd = 10, verbose = FALSE) {
  if (is.data.frame(data)) data <- as.matrix(data)
  criterion <- match.arg(criterion)
  n <- nrow(data) # sample size
  D <- ncol(data) # dimension
  maxit <- n
  ic <- double(maxit)
  icmin <- Inf
  info <- matrix(nrow = maxit, ncol = 3,
                 dimnames = list(seq_len(maxit), c("m", "ll", criterion)))
  for (i in seq_len(maxit)) {
    validmix <- FALSE
    for (j in seq_len(maxrd)) {
      r <- EMDir(data, i, bound)
      if (all(r$mix$alpha < 1e4) && is.finite(r$ll)) {
        validmix <- TRUE
        break
      }
    }
    if (!validmix) {
      warning("Can not find a valid finite mixture result.")
      return(NULL)
    }
    r <- collapse_points(r$mix, data)
    # Merge similar mixture components and delete components with tiny pr
    p <- D * nrow(r$mix$alpha) - 1
    ic[i] <- -2 * r$ll +
      switch(criterion,
             AIC = {2 * p},
             BIC = {log(n) * p},
             AICc = {2 * p + 2 * p * (p + 1) / max(n - p - 1, 0)})
    info[i, ] <- c(i, r$ll, ic[i])
    if (maxit == 1 || ic[i] < icmin) {
      icmin <- ic[i]
      imin <- i
      mix <- r$mix
      ll <- r$ll
    }
    if (verbose) {
      cat("Iter.#", i,
          ":\tm = ", i,
          "\tll = ", r$ll, "\t", criterion, " = ", ic[i], "\n", sep = "")
    }
    wdmin <- min(info[imin:min(n, imin + mtol), criterion])
    if (is.na(wdmin)) wdmin <- Inf
    if (ic[i] == Inf || wdmin == icmin) break
  }
  list(info = info[1:i, ], m = imin, ll = ll, criterion = icmin, mix = mix)
}
