ISErl <- function(obj, test) {
  a <- obj$alpha
  n <- nrow(a)
  Ca <- Cdir(a, log = TRUE)
  outCaa <- outer(Ca, Ca, "+")
  Caa <- Cdir_plus(a, a, log = TRUE)
  term1 <- if (class(obj) == "dirmix") {
    pr <- obj$pr
    sum(exp(outCaa - Caa) * outer(pr, pr))
  } else exp(log_sum_exp(outCaa - Caa) - 2 * log(n))
  term2 <- mean(dspx(test, obj))
  term1 - 2 * term2
}

rlCV <- function(data, k = 10, seed = NULL) {
  if (is.data.frame(data)) data <- as.matrix(data)
  data <- closure(data)
  if (!is.null(seed)) set.seed(seed)
  fid <- sample(rep(seq_len(k), length = nrow(data)))
  mthd <- c("PI", "CV", "RS", "LK", "FM")
  ms <- matrix(nrow = k, ncol = length(mthd), dimnames = list(NULL, mthd))
  vh <- c(seq(0.06, 0.01, len = 50), seq(0.01, 0.001, len = 30)[-1])
  for (i in seq_len(k)) {
    print(i)
    train <- data[fid != i, ]
    test <- data[fid == i, ]

    fm <- optfm(train)
    fmix <- fm$mix
    invalid <- is.null(fm) || any(fmix$alpha < 1)
    if (invalid) {
      ms[i, "seed"] <- seed
      print("Invalid finite mixture.")
      next
    }
    pih <- optMISEh(fmix, nrow(train)) # PI
    # CV & RS
    t1 <- sapply(vh, function(h) CV1(h, train))
    t2_lscv <- sapply(vh, function(h) LSCV2(h, train))
    t2_scv <- sapply(vh, function(h) SCV2(h, train, fmix))
    lscvh <- vh[which.min(t1 - 2 * t2_lscv)]
    scvh <- vh[which.min(t1 - 2 * t2_scv)]
    lkh <- optLCVh(train) # LK

    bw <- c(pih$h, lscvh, scvh, lkh$h)
    knre <- sapply(bw, function(h) ISErl(dirkn(train, h), test))
    fmre <- ISErl(fmix, test)
    ms[i, ] <- c(knre, fmre)
  }
  if (!is.null(seed)) ms <- cbind(seed, ms)
  ms
}
