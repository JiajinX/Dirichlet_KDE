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
  mthd <- c("PI", "CV", "LK", "FM")
  ms <- matrix(nrow = k, ncol = length(mthd))
  colnames(ms) <- mthd
  for (i in seq_len(k)) {
    print(i)
    train <- data[fid != i, ]
    test <- data[fid == i, ]

    fm <- optfm(train)
    pih <- optMISEh(fm$mix, nrow(train))
    cvh <- optLSCVh(train)
    lkh <- optLCVh(train)

    knpi <- dirkn(train, pih$h)
    kncv <- dirkn(train, cvh$h)
    knlk <- dirkn(train, lkh$h)

    ms[i, ] <- c(ISErl(knpi, test), ISErl(kncv, test),
                 ISErl(knlk, test), ISErl(fm$mix, test))
  }
  if (!is.null(seed)) ms <- cbind(seed, ms)
  ms
}
