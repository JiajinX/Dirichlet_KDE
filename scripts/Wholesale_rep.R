library(caret)

get_acc <- function(est, prior, test) {
  cotest <- closure(test[, 3:8])
  dc <- sapply(est, function(e) dspx(cotest, e))
  likeli <- dc * prior
  post <- closure(likeli)
  est_Channel <- apply(post, 1, which.max)
  mean(est_Channel == test$Channel)
}

classify_wholesale <- function(wholesale, interval = c(1, 0), rep = 1, p = 0.8) {
  cn <- c("seed", "FM", "MISE", "LSCV", "LCV",
          "FMm1", "FMm2", "MISEh1", "MISEh2",
          "LSCVh1", "LSCVh2", "LCVh1", "LCVh2")
  clsre <- matrix(nrow = rep, ncol = length(cn))
  colnames(clsre) <- cn
  for (i in seq_len(rep)) {
    print(i)
    seed <- i
    set.seed(seed)
    train_idx <- createDataPartition(wholesale$Channel, p = p, list = FALSE)
    train <- wholesale[train_idx, ]
    cotrain <- closure(train[, 3:8])
    cotrain <- split(cotrain, train$Channel)

    fm <- lapply(cotrain, function(cls) optfm(cls))
    fmix <- lapply(fm, `[[`, "mix")

    mise <- lapply(1:2, function(i) optMISEh(fmix[[i]], nrow(cotrain[[i]]), interval = interval))
    misekn <- lapply(1:2, function(i) dirkn(cotrain[[i]], mise[[i]]$h))

    lscv <- lapply(cotrain, function(x) optLSCVh(x, interval = interval))
    lscvkn <- lapply(1:2, function(i) dirkn(cotrain[[i]], lscv[[i]]$h))

    lcv <- lapply(cotrain, function(x) optLCVh(x, interval = interval))
    lcvkn <- lapply(1:2, function(i) dirkn(cotrain[[i]], lcv[[i]]$h))

    prior <- c(closure(tabulate(train$Channel)))
    test <- wholesale[-train_idx, ]
    fm_acc <- get_acc(fmix, prior, test)
    mise_acc <- get_acc(misekn, prior, test)
    lscv_acc <- get_acc(lscvkn, prior, test)
    lcv_acc <- get_acc(lcvkn, prior, test)

    clsre[i, ] <- c(seed, fm_acc, mise_acc, lscv_acc, lcv_acc,
                    sapply(fm, `[[`, "m"), sapply(misekn, `[[`, "h"),
                    sapply(lscvkn, `[[`, "h"), sapply(lcvkn, `[[`, "h"))
  }
  as.data.frame(clsre)
}


wholesale <- read.csv("data/app/Wholesale.csv")
clsre <- classify_wholesale(wholesale, c(0.05, 0.005), 3)
