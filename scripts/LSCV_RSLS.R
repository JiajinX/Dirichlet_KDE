invisible(lapply(dir("R", full.names = TRUE), source))

cvcurve <- function(no, n = 100, r = 0, vh, nsamp = 1e4, rep = 1) {
  c <- 2
  mix0 <- Dir_mix(no)
  mlscv <- mscv <- mrfcv <- matrix(nrow = rep, ncol = length(vh) + 1)

  for (i in seq_len(rep)) {
    seed <- i * 2e2
    data <- rmix(n, mix0, seed)

    # FM
    fm <- optfm(data)
    fmix <- fm$mix
    invalid <- is.null(fm) || any(fmix$alpha <= c)
    if (invalid) {
      miae[i, "seed"] <- mise[i, "seed"] <- seed
      print("Invalid finite mixture.")
      next
    }
    # CV
    t1 <- sapply(vh, function(h) CV1(h, data, r, c))
    t2_lscv <- sapply(vh, function(h) LSCV2(h, data, r, c))
    t2_scv <- sapply(vh, function(h) SCV2(h, data, fmix, r, c))
    t2_rfcv <- sapply(vh, function(h) SCV2(h, data, mix0, r, c))
    lscv <- t1 - 2 * t2_lscv
    scv <- t1 - 2 * t2_scv
    rfcv <- t1 - 2 * t2_rfcv

    mlscv[i, ] <- c(seed, lscv)
    mscv[i, ] <- c(seed, scv)
    mrfcv[i, ] <- c(seed, rfcv)
  }
  re <- list(lscv = mlscv, rsls = mscv, orsls = mrfcv)
  attr(re, "h") <- vh
  re
}

no <- 2
n <- 500
r <- 0 # 0,1: density / first derivative
vh <- seq(0.02, 0.005, len = 31)
rep <- 100
cvcurve(no, n, r, vh, rep = rep)
# Results are saved in: results/LSCV_RSLS/
