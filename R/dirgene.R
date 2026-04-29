# Sort components of Dirichlet mixture
sort.dirmix <- function(mix, decreasing = FALSE) {
  if (length(mix$pr) <= 1) return(mix)
  logmbeta <- attr(mix$alpha, "logmbeta")
  o <- order(mix$alpha[, 1], decreasing = decreasing)
  mix$alpha <- mix$alpha[o, ]
  mix$pr <- mix$pr[o]
  attr(mix, "mode") <- attr(mix, "mode")[o, ]
  attr(mix$alpha, "logmbeta") <- logmbeta[o]
  mix
}

# Constant term of Dirichlet distribution
Cdir <- function(a, log = FALSE) {
  if (is.vector(a)) a <- t(a)
  r <- lgamma(rowSums(a)) - rowSums(lgamma(a))
  if (!log) r <- exp(r)
  r
}

# Create a Dirichlet mixture given alpha and proportion of all components
dirmix <- function(a = rep(2, 3), pr = 1, sort = TRUE) {
  if (is.vector(a)) a <- t(a)
  m <- nrow(a)
  D <- ncol(a)

  pr <- rep(pr, len = m)
  pr <- pr / sum(pr)
  at <- a - 1 # To avoid some numeric issue in R
  if (any(at < -1e-8)) stop("Error: alpha < 1")
  # Alpha >= 1 for concave and unimodal
  at[at < 0] <- 0
  mode <- at / rowSums(at)
  a <- at + 1

  if (is.null(attr(a, "logmbeta"))) attr(a, "logmbeta") <- Cdir(a, log = TRUE)
  # If no log-multibeta or alpha has been modified
  mix <- list("alpha" = a, "pr" = pr)
  attr(mix, "mode") <- mode
  class(mix) <- "dirmix"
  if (sort) mix <- sort(mix)
  mix
}

# Create a Dirichlet kernel given x, h and c
dirkn <- function(x, h, c = NULL) {
  if (is.vector(x)) x <- t(x)
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(c)) c <- 1
  o <- order(x[, 1])
  x <- x[o, , drop = FALSE]
  kn <- list(alpha = x/h + c, x = x, h = h, c = c)
  class(kn) <- "dirkn"
  kn
}

# Generate random data from Dirichlet mixture distribution
rdirmix <- function(n, mix) {
  m <- nrow(mix$alpha)
  D <- ncol(mix$alpha)
  idx <- sample.int(m, n, prob = mix$pr, replace = TRUE)
  data <- matrix(nrow = n, ncol = D)
  start <- 1 # index
  for (i in seq_len(m)) {
    cnt <- sum(idx == i)
    if (cnt == 0) next
    end <- start + cnt - 1
    data[start:end, ] <- gtools::rdirichlet(cnt, mix$alpha[i, ])
    start <- end + 1
  }
  data[sample(n), ]
}

rdircomp <- function(n, mix) {
  alpha <- mix$alpha
  lapply(seq_len(nrow(alpha)), function(i) {
    gtools::rdirichlet(n, alpha[i, ])
  })
}

# outer = FALSE: calculate vector log D(x_i; α_i) with length n
# outer = TRUE: calculate matrix log D(x_i; α_j) with n * m
ddir <- function(x, a, log = FALSE, outer = FALSE) {
  if (is.vector(x)) x <- t(x)
  if (is.vector(a)) a <- t(a)
  if (any(is.nan(c(x, a)))) stop("The input has NaN.")
  n <- nrow(x)
  m <- nrow(a)

  logx <- attr(x, "log")
  logmbeta <- attr(a, "logmbeta")
  if (is.null(logx)) logx <- log(x)
  if (is.null(logmbeta)) logmbeta <- Cdir(a, log = TRUE) # Constant part

  if (!outer) { # Got a vector by adopting cycling rule on x, a
    a <- a[rep_len(1:m, n), , drop = FALSE]
    logmbeta <- rep_len(logmbeta, n)
  } else { # Got a matrix like outer function
    a <- a[rep(1:m, each = n), , drop = FALSE]
    logmbeta <- rep(logmbeta, each = n)
    logx <- logx[rep(1:n, m), , drop = FALSE]
  }
  xa <- (a - 1) * logx
  xa[is.na(xa)] <- 0 # When x = 0 and a = 1, x^(a-1) = 1 => (a-1) * logx = 0
  d <- logmbeta + rowSums(xa)
  if (outer) dim(d) <- c(n, m)
  if (!log) d <- exp(d)
  d
}

# Calculate c(log sum(π_j * D(x_1; α_j)), ..., log sum(π_j * D(x_n; α_j)))
# Length of vector = the number of data
ddirmix <- function(x, mix, log = FALSE) {
  if (is.vector(x)) x <- t(x)
  j <- mix$pr == 0
  if (any(j)) { # Drop mass with 0
    mix$alpha <- mix$alpha[!j, , drop = FALSE]
    mix$pr <- mix$pr[!j]
    attr(mix$alpha, "logmbeta") <- attr(mix$alpha, "logmbeta")[!j]
  }
  unipi <- unique(mix$pr)
  logpi <- if (length(unipi) == 1) {
    log(unipi)
  } else rep(log(mix$pr), each = nrow(x))
  logd <- ddir(x, mix$alpha, log = TRUE, outer = TRUE) + logpi
  # log D(x_i;α_j) + log π_j = log [π_j * D(x_i;α_j)]
  ma <- lsei::matMaxs(logd)
  ma[ma == -Inf] <- 0 # Otherwise -Inf-(-Inf) = NaN ???
  pid <- exp(logd - ma)
  pis <- rowSums(pid)
  r <- ma + log(pis)
  # c(log sum(π_j * D(x_1; α_j)), ..., log sum(π_j * D(x_n; α_j)))
  if (!log) r <- exp(r)
  r
}

ddirkn <- function(x, kn, log = FALSE) {
  if (is.vector(x)) x <- t(x)
  logd <- ddir(x, kn$alpha, log = TRUE, outer = TRUE)
  ma <- lsei::matMaxs(logd)
  ma[ma == -Inf] <- 0
  r <- log(rowMeans(exp(logd - ma))) + ma
  if (!log) r <- exp(r)
  r
}

ddirpred <- function(x, obj, log = FALSE) {
  d <- if (class(obj) == "dirkn") {
    ddirkn(x, obj, log = log)
  } else ddirmix(x, obj, log = log)
  d
}
