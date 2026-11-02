mean_shift <- function(obj, x0, a, pth = FALSE, stepmax = 0.01,
                       maxit = 1e3, tol = 1e-6, add = FALSE, col = "grey") {
  if (is.vector(x0)) x0 <- t(x0)
  if (is.data.frame(x0)) x0 <- as.matrix(x0)
  if (!add && pth) simplexpanel()
  xt <- x0   # Current locations
  # Rows that still need updating in the next round
  update <- rep(TRUE, nrow(xt))
  for (i in seq_len(maxit)) {
    up_idx <- which(update)
    if (length(up_idx) == 0) break
    xup <- xt[up_idx, , drop = FALSE]
    ngd <- pddir(xup, obj) / ddirpred(xup, obj)
    ngd_norm <- sqrt(rowSums(ngd^2))
    conv <- ngd_norm <= tol
    if (all(conv)) break
    if (any(conv)) {
      update[up_idx[conv]] <- FALSE
      up_idx <- up_idx[!conv]
      xup <- xup[!conv, , drop = FALSE]
      ngd <- ngd[!conv, , drop = FALSE]
    }
    ratio <- matrix(Inf, nrow = nrow(xup), ncol = ncol(xup))
    neg <- ngd < 0
    ratio[neg] <- xup[neg] / (-ngd[neg])
    amax <- apply(ratio, 1, min)
    ngd_norm <- sqrt(rowSums(ngd^2))
    xt[up_idx, ] <- xup + ngd * pmin(a, 0.3 * amax, stepmax / ngd_norm)
    if (pth) arrows_simplex(xup, xt[up_idx, ],
                            col = col, add = TRUE, length = 0, tol = 0)
  }
  xt
}

# Find local maximum by modal EM
MEM <- function(st, obj, tem = 5e2, tol = 1e-13,
                pth = FALSE, add = FALSE, col = "grey50") {
  if (is.vector(st)) st <- t(st)
  a <- obj$alpha
  a_1 <- a - 1
  pr <- if (class(obj) == "dirmix") obj$pr else 1 / nrow(obj$alpha)
  logpr <- log(pr)
  if (pth && !add) simplexpanel()
  # if (pth) points_simplex(st, TRUE, col)
  x <- st
  act <- rep(TRUE, nrow(x))
  eps <- 1e-15
  for (i in seq_len(tem)) {
    xact <- x[act, , drop = FALSE]
    fk <- ddir(xact, a, log = TRUE, outer = TRUE)
    pid <- t(t(fk) + logpr)
    ma <- lsei::matMaxs(pid)
    logpis <- ma + log(rowSums(exp(pid - ma)))
    pk <- exp(pid - logpis)
    bk <- pk %*% a_1
    xnew <- bk / rowSums(bk)
    if (pth) arrows_simplex(xact, xnew, TRUE, col = col, length = 0)

    Qnew <- rowSums(bk * log(pmax(xnew, eps)))
    Qact <- rowSums(bk * log(pmax(xact, eps)))
    convact <- abs(Qnew - Qact) < tol
    x[act, ] <- xnew
    act[which(act)[convact]] <- FALSE
    if (!any(act)) break
  }
  x
}

getcluster <- function(ed, digits = 4) {
  red <- round(ed, digits)
  cc <- unique(red)
  cno <- match(do.call(paste, as.data.frame(red)),
               do.call(paste, as.data.frame(cc)))
  cluster <- list(centers = cc, cluster = cno)
  class(cluster) <- "meanshift"
  cluster
}

dvalley_map <- function(clu, obj, size = NULL, k = 3) {
  cc <- clu$centers
  cno <- clu$cluster
  cn <- table(cno)
  if (is.null(size)) size <- length(cno) * 0.05
  sci <- which(cn <= size)
  dst <- as.matrix(dist(cc))
  r <- seq(0.01, 0.99, len = 99)
  den <- dst
  den[, ] <- NA
  for (i in sci) {
    cnear <- setdiff(order(dst[i, ]), sci)[1:k]
    for (j in cnear) {
      pth <- outer(1 - r, cc[i, ]) + outer(r, cc[j, ])
      d <- dspx(pth, obj)
      den[i, j] <- min(d)
    }
  }
  map <- apply(den[sci, , drop = FALSE], 1, which.max)
  map <- cbind(as.integer(names(map)), as.integer(map))
  colnames(map) <- c("From", "To")
  map
}

dist_pair <- function(x, y) {
  if (is.vector(x)) x <- t(x)
  if (is.vector(y)) y <- t(y)
  x2 <- rowSums(x^2)
  y2 <- rowSums(y^2)
  dist <- outer(x2, y2, "+") - 2 * x %*% t(y)
  sqrt(dist)
}

dist_merge <- function(x, tol) {
  dist <- as.matrix(dist(x))
  idx <- which(dist < tol & upper.tri(dist), arr.ind = TRUE)
  hc <- hclust(stats::dist(x), method = "complete")
  cno <- cutree(hc, h = tol)
  cc <- aggregate(x, by = list(group = cno), FUN = mean)
  cluster <- list(centers = closure(cc[, -1]), cluster = cno)
  class(cluster) <- "meanshift"
  cluster
}

distmin_map <- function(clu, x, size = NULL) {
  cc <- clu$centers
  cno <- clu$cluster
  cn <- table(cno)
  if (is.null(size)) size <- length(cno) * 0.05
  sci <- which(cn <= size)
  bci <- which(cn > size)
  dist <- matrix(nrow = length(sci), ncol = length(bci),
                 dimnames = list(sci, bci))
  for (i in seq_along(sci)) {
    scx <- x[cno == sci[i], , drop = FALSE]
    for (j in seq_along(bci)) {
      bcx <- x[cno == bci[j], , drop = FALSE]
      d <- dist_pair(scx, bcx)
      dist[i, j] <- mean(apply(d, 1, min))
    }
  }
  map <- apply(dist, 1, which.min)
  map <- cbind(as.integer(names(map)), as.integer(map))
  colnames(map) <- c("From", "To")
  map
}

mgcluster <- function(clu, map) {
  cc <- clu$centers
  cno <- clu$cluster
  change <- map[, "From"]
  cc <- cc[-change, , drop = FALSE]
  cno[cno %in% change] <- map[match(cno[cno %in% change], map[, "From"]), "To"]
  lbs <- unique(cno)
  lbsmap <- setNames(order(lbs), lbs)
  cno <- as.integer(lbsmap[as.character(cno)])
  cluster <- list(centers = cc, cluster = cno)
  class(cluster) <- "meanshift"
  cluster
}
