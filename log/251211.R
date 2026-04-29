invisible(lapply(dir("R", full.names = TRUE), source))

# MISE vs h =====================================
c <- 1.5
n <- 5
B <- 50
neva <- 5e4
mix <- simumix(11, c = 1.5)
mix <- gene_compsamp(mix)
mix <- gene_data(n, mix, B)
mix <- gene_eval(mix, r = 1, neva)
mix <- gene_grid(mix, r = 1, split = 100)

vh <- seq(0.2, 0.1, len = 15)
Mth <- DMISEfh(vh, mix, n, c = c)
plot_mise(Mth)
Mmc <- MISEh(vh, mix, MISEmc, r = 1, verbose = T)
Mdt <- MISEh(vh, mix, MISEdt, r = 1, verbose = T)

# pdf("figure/251211/Dir8.pdf")
plot_mise(Mth, ylim = range(c(Mth, Mmc$MISE)))
plot_mise(Mmc$MISE, col = 2, add = T)
plot_mise(Mdt$MISE, col = "blue", add = T)
pmc <- as.numeric(MISEtest(Mmc$ISE, Mth) > 0.05)
pdt <- as.numeric(MISEtest(Mdt$ISE, Mth) > 0.05)
title(main = paste0(paste0(pmc, collapse = " "), "\n",
                    paste0(pmc, collapse = " ")))
legend("top", legend = c("theory", "mc", "grid"),
       col = c("black", "red", "blue"), lty = 1, pch = 1, bty = "n")
dev.off()


# MISE vs h vs n =====================================
Dircheck <- function(n, vh, mix, r, c = NULL, B) {
  if (r == 1 && is.null(c)) c <- 1.5
  MISE <- DMISEfh(vh, mix, n, c)
  MISE_min <- min(MISE)
  h_min <- vh[which.min(MISE)]
  cat("\n")
  mix <- gene_data(n, mix, B)
  isemc <- MISEmc(h_min, mix, r, c, verbose = TRUE)
  pv <- t.test(attr(isemc, "ISE"), mu = MISE_min)$p.value
  attr(MISE_min, "search") <- MISE
  list(n = n, h = h_min, MISE_theo = MISE_min, MISE_mc = isemc, pv = pv)
}

r <- 1
B <- 50
neva <- 5e4
# mix <- dirmix(a = c(7, 13, 5))
mix <- Dir_mix(4)
mix <- gene_compsamp(mix)
mix <- gene_eval(mix, r, neva = neva)
vh <- seq(0.04, 0.01, len = 30)
vn <- c(50, 100, 200, 350, 500)

# z <- DMISEfh(vh, mix, 500, c = 1.5)
# plot_mise(z)

re <- list()
for (i in seq_along(vn)) { # For different h
  print(vn[i])
  re[[i]] <- Dircheck(vn[i], vh, mix, B = B, r = r) # r = 0
}
rdf <- do.call(rbind, lapply(re, function(x) {
  data.frame(n = x$n, h = x$h,
             MISE_theo = x$MISE_theo, MISE_mc = x$MISE_mc, pv = x$pv)
}))

# Plot of optimal h =====================================
scale <- 8e-4
# mix <- dirmix(a = c(7, 13, 5))
mix <- Dir_mix(4)
plot(mix, prop = F)
vfplot(mix, add = T, scale = scale)

# nhlist <- readRDS("data/251127/Dir_df")
nhlist <- readRDS("data/251127/Dirmix4_df")
n <- 500
h <- nhlist$h[nhlist$n == n]
c <- 1.5
data <- rdirmix(n, mix)
kn <- dirkn(data, h, c = c)
plot(kn, prop = F)
vfplot(kn, add = T, scale = scale)
dirmode(kn, st = attr(mix, "mode"), plot = 3, add = T)

# Different of gradient vector
vfplot(mix, scale = scale)
vfplot(kn, col = 2, scale = scale, add = T)
dirmode(mix, plot = 1, col = 1, add = T)
dirmode(kn, st = attr(mix, "mode"), plot = 1, col = "blue", add = T)
vfdiffplot(mix, kn, scale = scale, add = F)

# 3D plot
scale <- 5e-4
vfplot3D(mix, scale = scale)
vfplot3D(kn, col = 2, scale = scale, add = T)

