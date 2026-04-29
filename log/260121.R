# Mean shift clustering

invisible(lapply(dir("R", full.names = TRUE), source))

# Function from x0 to the modes of mixture
findmode <- function(st, obj, steplen = NULL, maxit = 1e3, tol = 1e-13,
                     plot = 0, add = FALSE, col = "grey30") {
  if (is.vector(st)) st <- t(st)
  if (is.null(steplen)) stop("No step size provided.")
  if (!add && plot != 0) simplexpanel()
  if (plot > 1) points_simplex(st, add = TRUE, col = col, pch = 20)
  ed <- st
  ed[] <- NA
  idx <- seq_len(nrow(st))
  for (i in seq_len(maxit)) {
    sti <- st[idx, , drop = FALSE]
    nDf <- pddir(sti, obj) / ddirpred(sti, obj)
    ed[idx, ] <- edi <- stepup(sti, steplen * nDf)
    if (plot == 3) {
      arrows_simplex(sti, edi, col = col, add = TRUE, length = 0, tol = 0)
    }
    idx <- lsei::matMaxs(ed - st) > tol
    if (sum(idx) == 0) break
    st <- ed
    print(i)
  }
  if (plot > 0) points_simplex(ed, add = TRUE,
                               pch = 8, cex = 1.2, lwd = 2, col = 2)
  ed
}

getcluster <- function(ed, digits = 5) {
  red <- round(ed, digits)
  cc <- unique(red)
  cno <- as.integer(interaction(as.data.frame(red), drop = TRUE))
  cluster <- list(centers = cc, cluster = cno)
  class(cluster) <- "meanshift"
  cluster
}

plot.meanshift <- function(clu, st, col = NULL,
                           pth = FALSE, obj = NULL, steplen = NULL) {
  m <- max(clu$cluster)
  if (is.null(col)) col <- viridis::viridis(m)
  simplexpanel()
  for (i in seq_len(m)) {
    if (pth) {
      findmode(st[clu$cluster == i, ], obj, steplen,
               plot = 3, add = TRUE, col = col[i])
    } else {
      points_simplex(st[clu$cluster == i, ], col = col[i], add = TRUE)
      points_simplex(clu$centers, add = TRUE,
                     pch = 8, cex = 1.2, lwd = 2, col = 2)
    }
  }
}

mix0 <- Dir_mix(4)
plot(mix0, prop = F)
vfplot(mix0, add = T)
n <- 250
x <- rdirmix(n, mix0)
est <- optfm(x)
plot(est$mix, prop = F)
vfplot(est$mix, add = T)

estmix <- gene_compsamp(est$mix)
dh <- optDMISEh(vh = seq(0.03, 0.01, len = 10), estmix, n, c = 1.5, verbose = TRUE)
knest <- dirkn(x = x, h = dh$h)
plot(knest, prop = F)
vfplot(knest, add = T)

ed <- findmode(x, knest, steplen = knest$h/5, plot = 0)
clu <- getcluster(ed)
plot(clu, x)
plot(clu, x, pth = T, obj = knest, steplen = knest$h/5)
