invisible(lapply(dir("R", full.names = TRUE), source))
pth <- "results/MISE_MC/"

# pdf("figures/6.1-MISE_curve.pdf", width = 6, height = 5)
cols = c("red", "black")
for (r in 0:1) {
  xmin <- ifelse(r == 0, 0.005, 0.01)
  pthname <- paste0(pth, "dir_r", r, "_3D_500")
  emi <- readRDS(pthname) # Exact evaluation
  dmi <- colMeans(readRDS(paste0(pthname, "_mc"))) # Direct approximation
  vh <- attr(emi, "h")
  plot.new()
  idx <- vh >= xmin
  plot.window(xlim = c(xmin, 0.02), ylim = range(dmi[idx]))
  title(xlab = "Bandwidth", ylab = "MISE", cex.lab = 1.2)
  box(); axis(1, cex.axis = 0.8); axis(2, cex.axis = 0.9)
  id <- ifelse(r == 0, 999, 125) # Seed number
  emid <- emi[id, ]
  lines(vh, emid, col = cols[2])
  abline(v = vh[which.min(emid)], col = cols[2])
  lines(vh, dmi, lty = 2, col = cols[1])
  abline(v = vh[which.min(dmi)], lty = 2, col = cols[1])
  legend("topright", legend = c("Direct", "Exact"), bg = "white",
         col = cols, lty = c(2, 1))
}
# dev.off()
