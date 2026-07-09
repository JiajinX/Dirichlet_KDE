# Function from x to the modes of mixture
# findmode <- function(st, obj, step_init = NULL, maxit = 1e3, maxrd = 4,
#                      tol = 1e-13, plot = 0, add = FALSE, col = "grey30") {
#   if (is.vector(st)) st <- matrix(st, nrow = 1)
#   if (is.data.frame(st)) st <- as.matrix(st)
#   if (is.null(step_init)) stop("No initial step size provided.")
#   if (!add && plot != 0) simplexpanel()
#   if (plot > 1) points_simplex(st, add = TRUE, col = col, pch = 20)
#   ed <- st   # Current locations
#   # Rows that still need updating in the next round
#   pending <- rep(TRUE, nrow(ed))
#   step_len <- step_init
#   for (r in seq_len(maxrd)) {
#     if (!any(pending)) break
#     active <- which(pending)
#     next_pending <- rep(FALSE, nrow(ed))
#     for (i in seq_len(maxit)) {
#       if (length(active) == 0) break
#       sti <- ed[active, , drop = FALSE]
#       nDf <- pddir(sti, obj) / ddirpred(sti, obj)
#       edi <- stepup(sti, step_len * nDf)
#       # Rows whose updated point contains 0
#       hit0 <- rowSums(edi == 0) > 0
#       # rows hitting 0 leave the current round,
#       # but continue in the next round with half step size
#       if (any(hit0)) next_pending[active[hit0]] <- TRUE
#       keep <- !hit0 # Rows that do not hit 0 can be updated
#       if (any(keep)) {
#         old_keep <- sti[keep, , drop = FALSE]
#         new_keep <- edi[keep, , drop = FALSE]
#         ed[active[keep], ] <- new_keep # Accept update
#         if (plot == 3) {
#           arrows_simplex(old_keep, new_keep,
#                          col = col, add = TRUE, length = 0, tol = 0)
#         }
#         # Check whether still moving
#         moved <- lsei::matMaxs(abs(new_keep - old_keep)) > tol
#         # Continue in this round if still moving, otherwise stop for this row
#         next_pending[active[keep]] <- moved
#         active <- active[keep][moved]
#       } else active <- integer(0)
#     }
#     pending <- next_pending
#     step_len <- step_len / 2
#   }
#   if (plot > 0) {
#     points_simplex(ed, add = TRUE,
#                    pch = 8, cex = 1.2, lwd = 2, col = 2)
#   }
#   ed
# }

# plot.meanshift <- function(clu, st, col = NULL,
#                            pth = FALSE, obj = NULL, steplen = NULL) {
#   m <- max(clu$cluster)
#   if (is.null(col)) col <- paletteer::paletteer_d("ggthemes::gdoc", m)
#   simplexpanel()
#   for (i in seq_len(m)) {
#     if (pth) {
#       findmode(st[clu$cluster == i, ], obj, steplen,
#                plot = 3, add = TRUE, col = col[i])
#     } else {
#       points_simplex(st[clu$cluster == i, ], col = col[i], add = TRUE)
#       points_simplex(clu$centers, add = TRUE,
#                      pch = 8, cex = 1.2, lwd = 2, col = "grey30")
#     }
#   }
# }

# Find local maximum by modal EM
MEM <- function(st, obj, tem = 5e2, tol = 1e-13,
                pth = FALSE, add = FALSE, col = "grey50") {
  if (is.vector(st)) st <- t(st)
  a <- obj$alpha
  a_1 <- a - 1
  pr <- if (class(obj) == "dirmix") obj$pr else 1 / nrow(obj$alpha)
  logpr <- log(pr)
  if (pth && !add) simplexpanel()
  if (pth) points_simplex(st, TRUE, col)
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

plot.meanshift <- function(clu, st, pth = FALSE, obj = NULL, contour = FALSE,
                           col = NULL, lab = NULL) {
  m <- max(clu$cluster)
  if (is.null(col)) col <- paletteer::paletteer_d("ggthemes::gdoc", m)
  simplexpanel(lab, labcex = 1.2)
  if (contour) plot(obj, data = st, add = TRUE, col = "grey60",
                    lv = seq(30, 90, len = 4) / 100, axl = 200)
  for (i in seq_len(m)) {
    sti <- st[clu$cluster == i, ]
    if (pth) {
      MEM(sti, obj, pth = TRUE, add = TRUE, col = col[i])
    } else {
      points_simplex(sti, col = col[i], add = TRUE)
    }
    points_simplex(clu$centers, add = TRUE,
                   pch = 3, cex = 1.2, col = 1)
  }
}
