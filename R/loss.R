# 1. Create proposal density q as a smoothed version of the true density
# 2. Generate large samples of q
IE_samp <- function(mix, n = 1e4, fct = 1.2) {
  isdir <- class(mix) == "dirmix"
  q <- mix
  if (isdir) {
    mode <- attr(mix, "mode")
    h <- fct / (rowSums(mix$alpha) - ncol(mix$alpha)) # Larger bandwidth
    q$alpha <- mode / h + 1
  } else { # normix
    m <- length(mix$pr)
    D <- ncol(mix$Sigma)
    spSig <- split(mix$Sigma, rep(1:m, each = D)) # Split sigma
    heavySig <- lapply(spSig, function(x) {
      dim(x) <- c(D, D)
      diag(x) <- diag(x) * fct
      x
    })
    q$Sigma <- do.call(rbind, heavySig)
  }
  sp <- rmix(n, q)
  if (!isdir) sp <- ilrts(sp)
  dmix <- dspx(sp, mix)
  idx <- dmix > 1e-10
  sp <- if (isdir) sp[idx, ] else select_ilrx(sp, idx)
  # Delete points with too tiny density
  dq <- dspx(sp, q)
  dmix <- dmix[idx]
  list(x = sp, dmix = dmix, dhv = dq)
  # Density of mix and heavier tail mixture
}

IE <- function(samp, est) {
  x <- samp$x
  if (inherits(samp$x, "ilrx")) x <- x$sx
  ed <- ddirpred(x, est)
  td <- samp$dmix
  hvd <- samp$dhv
  iae <- mean(abs(td - ed) / hvd)
  ise <- mean((td - ed)^2 / hvd)
  list(iae = iae, ise = ise)
}

# Sample from uniform distribution
DIE_samp <- function(mix, n = 1e4) {
  D <- if (inherits(mix, "dirmix")) ncol(mix$alpha) else ncol(mix$mu) + 1
  sp <- gtools::rdirichlet(n, rep(1, D))
  pdmix <- pdspx(sp, mix)
  idx <- rowSums(abs(pdmix) < 1e-10) == 0
  sp <- sp[idx, ]
  pdmix <- pdmix[idx, ]
  list(x = sp, pdmix = pdmix, vol = exp(-lfactorial(D - 1)))
}

DIE <- function(samp, est) {
  epd <- pddir(samp$x, est)
  tpd <- samp$pdmix
  vol <- samp$vol
  iae <- mean(rowSums(abs(tpd - epd))) * vol
  ise <- mean(rowSums((tpd - epd)^2)) * vol
  list(iae = iae, ise = ise)
}

# Measure the estimation
msdst <- function(type = c("dirichlet", "normal"), no, D = 3, n = 100,
                  nsamp = 1e4, rep = 1) {
  type <- match.arg(type)
  cn <- c("seed", "PI", "CV") # Estimation methods
  miae <- mise <- matrix(nrow = rep, ncol = length(cn),
                         dimnames = list(NULL, cn))
  mix0 <- switch(type,
                 dirichlet = Dir_mix(no, D),
                 normal = Norm_mix(no, D - 1))
  samp <- IE_samp(mix0, nsamp)
  # Density estimation
  for (i in seq_len(rep)) {
    print(i)
    seed <- i
    set.seed(seed)
    data <- rmix(n, mix0)
    if (type == "normal") data <- ilrts(data)$sx

    fm <- optfm(data)
    if (any(fm$mix$alpha > 1e2)) {
      print("Invalid finite mixture.")
      next
    }
    pih <- optMISEh(fm$mix, n)
    cvh <- optLSCVh(data)

    knpi <- dirkn(x = data, h = pih$h)
    kncv <- dirkn(x = data, h = cvh$h)
    IEpi <- IE(samp, knpi)
    IEcv <- IE(samp, kncv)
    miae[i, ] <- c(seed, IEpi$iae, IEcv$iae)
    mise[i, ] <- c(seed, IEpi$ise, IEcv$ise)
  }
  miae <- cbind(no, miae)
  mise <- cbind(no, mise)
  list(miae = miae, mise = mise)
}

msfst <- function(type = c("dirichlet", "normal"), no, D = 3, n = 100, c = NULL,
                  nsamp = 1e4, rep = 1) {
  type <- match.arg(type)
  if (is.null(c)) c <- 1.5
  cn <- c("seed", "PI", "CV") # Estimation methods
  miae <- mise <- matrix(nrow = rep, ncol = length(cn),
                         dimnames = list(NULL, cn))
  mix0 <- switch(type,
                 dirichlet = Dir_mix(no, D),
                 normal = Norm_mix(no, D - 1))
  samp <- DIE_samp(mix0, nsamp)
  # Density estimation
  for (i in seq_len(rep)) {
    print(i)
    seed <- i
    set.seed(seed)
    data <- rmix(n, mix0)
    if (type == "normal") data <- ilrts(data)$sx

    fm <- optfm(data)
    if (any(fm$mix$alpha > 1e2 | fm$mix$alpha <= 1.5)) {
      print("Invalid finite mixture.")
      next
    }
    pih <- optMISEh(fm$mix, n, r = 1, c = c)
    cvh <- optLSCVh(data, r = 1, c = c)

    knpi <- dirkn(x = data, h = pih$h, c = c)
    kncv <- dirkn(x = data, h = cvh$h, c = c)
    IEpi <- DIE(samp, knpi)
    IEcv <- DIE(samp, kncv)
    miae[i, ] <- c(seed, IEpi$iae, IEcv$iae)
    mise[i, ] <- c(seed, IEpi$ise, IEcv$ise)
  }
  miae <- cbind(no, miae)
  mise <- cbind(no, mise)
  list(miae = miae, mise = mise)
}
