# Function from x to the modes of mixture
findmode <- function(st, obj, step_init = NULL, maxit = 1e3, maxrd = 4,
                     tol = 1e-13, plot = 0, add = FALSE, col = "grey30") {
  if (is.vector(st)) st <- matrix(st, nrow = 1)
  if (is.data.frame(st)) st <- as.matrix(st)
  if (is.null(step_init)) stop("No initial step size provided.")
  if (!add && plot != 0) simplexpanel()
  if (plot > 1) points_simplex(st, add = TRUE, col = col, pch = 20)
  ed <- st   # Current locations
  # Rows that still need updating in the next round
  pending <- rep(TRUE, nrow(ed))
  step_len <- step_init
  for (r in seq_len(maxrd)) {
    if (!any(pending)) break
    active <- which(pending)
    next_pending <- rep(FALSE, nrow(ed))
    for (i in seq_len(maxit)) {
      if (length(active) == 0) break
      sti <- ed[active, , drop = FALSE]
      nDf <- pddir(sti, obj) / ddirpred(sti, obj)
      edi <- stepup(sti, step_len * nDf)
      # Rows whose updated point contains 0
      hit0 <- rowSums(edi == 0) > 0
      # rows hitting 0 leave the current round,
      # but continue in the next round with half step size
      if (any(hit0)) next_pending[active[hit0]] <- TRUE
      keep <- !hit0 # Rows that do not hit 0 can be updated
      if (any(keep)) {
        old_keep <- sti[keep, , drop = FALSE]
        new_keep <- edi[keep, , drop = FALSE]
        ed[active[keep], ] <- new_keep # Accept update
        if (plot == 3) {
          arrows_simplex(old_keep, new_keep,
                         col = col, add = TRUE, length = 0, tol = 0)
        }
        # Check whether still moving
        moved <- lsei::matMaxs(abs(new_keep - old_keep)) > tol
        # Continue in this round if still moving, otherwise stop for this row
        next_pending[active[keep]] <- moved
        active <- active[keep][moved]
      } else active <- integer(0)
    }
    pending <- next_pending
    step_len <- step_len / 2
  }
  if (plot > 0) {
    points_simplex(ed, add = TRUE,
                   pch = 8, cex = 1.2, lwd = 2, col = 2)
  }
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
  if (is.null(col)) col <- paletteer::paletteer_d("ggthemes::gdoc", m)
  simplexpanel()
  for (i in seq_len(m)) {
    if (pth) {
      findmode(st[clu$cluster == i, ], obj, steplen,
               plot = 3, add = TRUE, col = col[i])
    } else {
      points_simplex(st[clu$cluster == i, ], col = col[i], add = TRUE)
      points_simplex(clu$centers, add = TRUE,
                     pch = 8, cex = 1.2, lwd = 2, col = "grey30")
    }
  }
}
