msdst <- function(type = c("dirichlet", "normal"), no, D = 3, n = 100, r = 0,
                  nsamp = 1e4, rep = 1) {
  type <- match.arg(type)
  if (r == 0) {
    c <- 1
    IEsamp <- IE_samp
    IEmeas <- IE
  } else {
    c <- 1.5
    IEsamp <- DIE_samp
    IEmeas <- DIE
  }
  cn <- c("seed", "PI", "CV", "PIh", "CVh") # Estimation methods
  miae <- mise <- matrix(nrow = rep, ncol = length(cn), dimnames = list(NULL, cn))
  mix0 <- switch(type,
                 dirichlet = Dir_mix(no, D),
                 normal = Norm_mix(no, D - 1))
  samp <- IEsamp(mix0, nsamp)
  # Density estimation
  for (i in seq_len(rep)) {
    seed <- id * 2e2 + i
    cat(i, ": ", seed, "\n")
    set.seed(seed)
    data <- rmix(n, mix0)
    if (type == "normal") data <- ilrts(data)$sx

    fm <- optfm(data)
    if (any(fm$mix$alpha > 1e4 | fm$mix$alpha <= c)) {
      miae[i, "seed"] <- mise[i, "seed"] <- seed
      print("Invalid finite mixture.")
      next
    }
    pih <- optMISEh(fm$mix, n, r)
    cvh <- optLSCVh(data, r)

    knpi <- dirkn(data, pih$h, c)
    kncv <- dirkn(data, cvh$h, c)
    IEpi <- IEmeas(samp, knpi)
    IEcv <- IEmeas(samp, kncv)
    miae[i, ] <- c(seed, IEpi$iae, IEcv$iae, pih$h, cvh$h)
    mise[i, ] <- c(seed, IEpi$ise, IEcv$ise, pih$h, cvh$h)
  }
  miae <- cbind(no, miae)
  mise <- cbind(no, mise)
  list(miae = miae, mise = mise)
}
