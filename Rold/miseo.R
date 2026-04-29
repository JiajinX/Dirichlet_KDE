ddirkn <- function(s, X, h, c = 1, log = FALSE, change = FALSE) {
  if (is.vector(s)) s <- t(s)
  if (is.vector(X)) X <- t(X)

  if (change) {
    logd <- ddir(X, s/h + c, log = TRUE, outer = TRUE)
    ma <- lsei::matMaxs(logd, dim = 2)
    ma[ma == -Inf] <- 0
    r <- log(colMeans(exp(sweep(logd, 2, ma)))) + ma
  } else {
    logd <- ddir(s, X/h + c, log = TRUE, outer = TRUE)
    ma <- lsei::matMaxs(logd)
    ma[ma == -Inf] <- 0
    r <- log(rowMeans(exp(logd - ma))) + ma
  }
  if (!log) r <- exp(r)
  r
}

sampo <- function(mix, mn = 1e5) {
  m <- length(mix$pr)
  D <- ncol(mix$alpha)
  samp <- do.call(rbind, rdircomp(mn, mix))
  attr(samp, "size") <- c(m, mn)
  attr(mix, "samp") <- samp
  attr(mix, "unisamp") <- gtools::rdirichlet(mn, rep(1, D))
  mix
}

Cdho <- function(h, mix) {
  Ca <- attr(mix$alpha, "logmbeta")
  unisamp <- attr(mix, "unisamp")
  C2uni <- attr(h, "C2uni")
  D <- ncol(mix$alpha)

  C2unia <- Cdir_plus(2 * unisamp / h + 1,  mix$alpha, log = TRUE)
  sum(mix$pr * exp(Ca - lgamma(D) + log(colMeans(exp(C2uni - C2unia)))))
}

# Symmetric
Psi2o <- function(h, mix) {
  D <- ncol(mix$alpha)
  m <- length(mix$pr)
  Ca <- attr(mix$alpha, "logmbeta")
  unisamp <- attr(mix, "unisamp")
  n <- nrow(unisamp)
  C2uni <- attr(h, "C2uni")

  outCaa <- outer(Ca, Ca, "+")
  up <- outer(c(outCaa), C2uni, "+")
  Cunia <- Cdir_plus(unisamp / h + 1,  mix$alpha, log = TRUE)
  down <- rep(c(Cunia), m) + c(Cunia[rep(seq_len(n), m), ])
  dim(down) <- c(n, m^2)
  r <- rowMeans(exp(up - t(down) - lgamma(D)))
  dim(r) <- c(m, m)
  r
}

MISEfo <- function(h, mix, n) {
  samp <- attr(mix, "samp")
  unisamp <- attr(mix, "unisamp")
  if (is.null(samp) || is.null(unisamp)) stop("No samples generated.")
  attr(h, "Cy") <- Cdir(samp / h + 1, log = TRUE)
  attr(h, "C2uni") <- 2 * Cdir(unisamp / h + 1, log = TRUE)

  cdh <- Cdho(h, mix)
  psi0 <- Psi0(mix)
  psi1 <- Psi1(h, mix)
  psi2 <- Psi2o(h, mix)

  c(cdh / n + t(mix$pr) %*% ((1 - 1/n) * psi2 - 2 * psi1 + psi0) %*% mix$pr)
}

MISEmco <- function(h, mix) {
  evapoint <- attr(mix, "evapoint")
  data <- attr(mix, "data")
  if (is.null(evapoint) || is.null(data))
    stop("No sample generated from the mixture.")

  eva <- evapoint$eva
  dmix <- evapoint$d
  B <- length(data)
  ISEs <- numeric(B)
  for (i in seq_len(B)) {
    df <- ddirkn(eva, data[[i]], h, change = TRUE)
    ISEs[i] <- mean((df - dmix)^2 / dmix)
  }
  MISE <- mean(ISEs)
  attr(MISE, "ISE") <- ISEs
  MISE
}

MISEmco <- function(h, mix) {
  D <- ncol(mix$alpha)
  evapoint <- attr(mix, "evapoint")
  data <- attr(mix, "data")
  if (is.null(evapoint) || is.null(data))
    stop("No sample generated from the mixture.")

  eva <- evapoint$eva
  dmix <- evapoint$d
  B <- length(data)
  ISEs <- numeric(B)
  for (i in seq_len(B)) {
    df <- ddirkn(eva, data[[i]], h, change = TRUE)
    ISEs[i] <- mean((df - dmix)^2) / gamma(D)
  }
  MISE <- mean(ISEs)
  attr(MISE, "ISE") <- ISEs
  MISE
}

MISEdto <- function(h, mix) {
  data <- attr(mix, "data")
  grid <- attr(mix, "gridinfo")
  if (is.null(grid) || is.null(data)) stop("No data or grid provided.")

  gd <- grid$x
  dmix <- grid$dx
  B <- length(data)
  ISEs <- numeric(B)
  for (i in seq_len(B)) {
    df <- ddirkn(gd, data[[i]], h, change = TRUE)
    ISEs[i] <- sum((df - dmix)^2) * grid$area
  }
  MISE <- mean(ISEs)
  attr(MISE, "ISE") <- ISEs
  MISE
}
