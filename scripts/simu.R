invisible(lapply(dir("R", full.names = TRUE), source))

msdst <- function(type = c("dirichlet", "normal"), no, D = 3, n = 100, r = 0,
                  nsamp = 1e4, rep = 1) {
  type <- match.arg(type)
  if (r == 0) {
    c <- 1
    IEsamp <- IE_samp
    IEmeas <- IE
    cn <- c("seed", "RF", "PI", "LSCV", "SCV", "LK",
            "RFh", "PIh", "LSCVh", "SCVh", "LKh")
  } else {
    c <- 1.5
    IEsamp <- DIE_samp
    IEmeas <- DIE
    cn <- c("seed", "RF", "PI", "LSCV", "SCV", "RFh", "PIh", "LSCVh", "SCVh")
  }
  miae <- mise <- matrix(nrow = rep, ncol = length(cn), dimnames = list(NULL, cn))
  mix0 <- switch(type,
                 dirichlet = Dir_mix(no, D),
                 normal = Norm_mix(no, D - 1))
  vh <- c(seq(0.06, 0.01, len = 50), seq(0.01, 0.001, len = 30)[-1])
  rf <- readRDS("results/raw/refh")
  rfh <- rf[rf$type == type & rf$no == no & rf$r == r & rf$D == D & rf$n == n, "h"]
  for (i in seq_len(rep)) {
    seed <- i * 2e2
    samp <- IEsamp(mix0, nsamp) # IE_samp / DIE_samp
    data <- rmix(n, mix0, seed)
    if (type == "normal") data <- ilrts(data)$sx
    # FM
    fm <- optfm(data)
    fmix <- fm$mix
    invalid <- is.null(fm) || ifelse(r == 0, any(fmix$alpha < c), any(fmix$alpha <= c))
    if (invalid) {
      miae[i, "seed"] <- mise[i, "seed"] <- seed
      print("Invalid finite mixture.")
      next
    }
    # PI
    pih <- optMISEh(fmix, n, r)
    # CV
    t1 <- sapply(vh, function(h) CV1(h, data, r, c))
    t2_lscv <- sapply(vh, function(h) LSCV2(h, data, r, c))
    t2_scv <- sapply(vh, function(h) SCV2(h, data, fmix, r, c))
    lscvh <- vh[which.min(t1 - 2 * t2_lscv)]
    scvh <- vh[which.min(t1 - 2 * t2_scv)]
    if (r == 0) {
      lkh <- optLCVh(data)
      knlk <- dirkn(data, lkh$h, c)
      IElk <- IEmeas(samp, knlk)
    }
    knrf <- dirkn(data, rfh, c)
    knpi <- dirkn(data, pih$h, c)
    knlscv <- dirkn(data, lscvh, c)
    knscv <- dirkn(data, scvh, c)
    IErf <- IEmeas(samp, knrf)
    IEpi <- IEmeas(samp, knpi) # IE/DIE
    IElscv <- IEmeas(samp, knlscv)
    IEscv <- IEmeas(samp, knscv)
    if (r == 0) {
      miae[i, ] <- c(seed, IErf$iae, IEpi$iae, IElscv$iae, IEscv$iae,
                     IElk$iae, rfh, pih$h, lscvh, scvh, lkh$h)
      mise[i, ] <- c(seed, IErf$ise, IEpi$ise, IElscv$ise, IEscv$ise,
                     IElk$ise, rfh, pih$h, lscvh, scvh, lkh$h)
    } else {
      miae[i, ] <- c(seed, IErf$iae, IEpi$iae, IElscv$iae, IEscv$iae,
                     rfh, pih$h, lscvh, scvh)
      mise[i, ] <- c(seed, IErf$ise, IEpi$ise, IElscv$ise, IEscv$ise,
                     rfh, pih$h, lscvh, scvh)
    }
  }
  miae <- cbind(no, miae)
  mise <- cbind(no, mise)
  list(miae = miae, mise = mise)
}

type <- "dirichlet" # dirichlet / normal
no <- 1 # 1~6: No. of test mixture
D <- 3 # 3,5: dimension
n <- 500 # 100,500: sample size
r <- 0 # 0,1: density / first derivative
rep <- 1e3 # Number of repetitions
simu <- msdst(type, no, D, n, r, rep = 1e3)
# The results are saved in: results/Simulation/
