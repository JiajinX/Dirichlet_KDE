invisible(lapply(dir("R", full.names = TRUE), source))

pth <- "results/LSCV_RSLS/"
lscv <- readRDS(paste0(pth, "lscv"))
rsls <- readRDS(paste0(pth, "rsls"))
orsls <- readRDS(paste0(pth, "orsls"))

xmin <- 0.005
vh <- attr(lscv, "h")
cvlst <- list(lscv = lscv, rsls = rsls, orsls = orsls)
mu <- lapply(cvlst, colMeans)
mu <- do.call(rbind, mu)

# pdf("figures/6.2-LSCV_RSLS.pdf", width = 6, height = 5)
# Curve
plot.new()
idx <- vh >= xmin
plot.window(xlim = c(xmin, 0.02), ylim = range(mu[, idx]))
title(xlab = "Bandwidth", ylab = "Mean loss", cex.lab = 1.2)
box()
axis(1, cex.axis = 0.8)
axis(2, cex.axis = 0.9)
cols <- c("black", "blue", "red")
for (i in 3:1) {
  col <- cols[i]
  lines(vh, mu[i, ], col = col, lty = ifelse(i == 1, 2, 1))
  abline(v = vh[which.min(mu[i, ])], col = col, lty = ifelse(i == 3, 1, 2))
}
legend("topright", legend = c("LSCV", "RSLS", "ORSLS"),
       lty = c(2, 1, 1), col = cols, bg = "white")

# Bar plot
alpha <- 0.3
cols_alpha <- adjustcolor(cols, alpha.f = alpha)
tabmat <- matrix(nrow = 3, ncol = length(vh),
                 dimnames = list(names(cvlst), as.character(vh)))
for (i in seq_along(cvlst)) {
  h <- vh[apply(cvlst[[i]], 1, which.min)]
  tabmat[i, ] <- as.numeric(table(factor(h, levels = vh)))
}
tabmat <- closure(tabmat)
tabmat <- tabmat[, rev(seq_len(ncol(tabmat)))]
show_id <- seq(1, length(vh), by = 2)
for (i in seq_len(nrow(tabmat))) {
  bp <- barplot(tabmat[i, ], add = i > 1, axes = i == 1, axisnames = FALSE,
                col = cols_alpha[i], border = NA,
                ylim = c(0, max(tabmat) * 1.1), cex.axis = 0.9,
                ylab = if (i == 1) "Probability" else "",  xlab = "", cex.lab = 1.2)
  if (i == 1) {
    axis(side = 1, at = bp[show_id], labels = format(rev(vh[show_id]), nsmall = 3), las = 2, cex.axis = 0.8)
  }
}
legend("topleft", legend = c("LSCV", "RSLS", "ORSLS"),
       fill = cols_alpha, bg = "white")
mtext("Bandwidth", side = 1, line = 3, cex = 1.2)
box()
# dev.off()
