invisible(lapply(dir("R", full.names = TRUE), source))

Dircheck <- function(n, vh, mix, r, c = NULL, B = 50) {
  if (r == 1 && is.null(c)) c <- 1.5
  MISE <- sapply(vh, function(i) {
    cat(round(i, 4), " ")
    if (r == 1) DMISEf(i, mix, n, c) else MISEf(i, mix, n)
  })
  cat("\n")
  names(MISE) <- round(vh, 4)
  MISE_min <- min(MISE)
  h_min <- vh[which.min(MISE)]

  mix <- gene_data(n, mix, B)
  isemc <- MISEmc(h_min, mix, r = r, c = c, verbose = TRUE)
  pv <- t.test(attr(isemc, "ISE"), mu = MISE_min)$p.value
  attr(MISE_min, "search") <- MISE
  list(n = n, h = h_min, MISE_theo = MISE_min, MISE_mc = isemc, pv = pv)
}

mix0 <- dirmix(a = c(7, 13, 5))
# mix0 <- Dir_mix(4)
# plot(mix0, prop = FALSE)
# vfmixplot(mix0, add = TRUE)
vh <- seq(0.02, 0.003, len = 30)
vn <- c(50, 100, 200, 350, 500)
neva <- 5e4
B <- 50
r <- 0
re <- list()
mix <- gene_compsamp(mix0)
mix <- gene_eval(mix, r, neva = neva)
for (i in seq_along(vn)) {
  print(vn[i])
  re[[i]] <- Dircheck(vn[i], vh, mix, B = B, r = r) # r = 0
}
rdf <- do.call(rbind, lapply(re, function(x) {
  data.frame(n = x$n, h = x$h,
             MISE_theo = x$MISE_theo, MISE_mc = x$MISE_mc, pv = x$pv)
}))
saveRDS(re, "data/251126/Dir")
saveRDS(rdf, "data/251126/Dir_df")

n <- 50
vh <- seq(0.05, 0.03, len = 10)
MISE <- sapply(vh, function(i) {
  cat(round(i, 4), " ")
  DMISEf(i, mix, n, c = 1.5)
})
plot(vh, MISE, xlim = rev(range(vh)), type = "o")
# data <- rdirmix(n, mix0)
# knplot(data, h = h_min, c = c, prop = FALSE)
# vfknplot(data, h = h_min, c = c, add = TRUE)

n <- 50
vh <- seq(0.025, 0.01, len = 10)
MISE <- sapply(vh, function(i) {
  cat(round(i, 4), " ")
  MISEf(i, mix, n)
})
plot(vh, MISE, xlim = rev(range(vh)), type = "o")


# =================================================
mix0 <- simumix(1)
n <- 5
mix <- gene_compsamp(mix0)
mix <- gene_eval(mix, 0, neva = 1e5)
mix <- gene_grid(mix, split = 50)
mix <- gene_data(n, mix, B = 50)
vh <- seq(0.06, 0.04, len = 20)
mise <- sapply(vh, function(i) {
  cat(round(i, 4), " ")
  MISEf(i, mix, n)
})
isemc <- MISEh(vh, mix, fise = MISEmc)
misemc <- isemc$MISE
isedt <- MISEh(vh, mix, fise = MISEdt)
misedt <- isedt$MISE

plot.new()
plot.window(xlim = rev(range(vh)), ylim = range(mise, misemc, misedt))
box(); axis(1); axis(2)
lines(vh, mise, col = "blue", type = "o")
lines(vh, misemc, col = "red", type = "o")
lines(vh, misedt, col = "darkgreen", type = "o")
MISEtest(isemc$ISE, mise)
MISEtest(isedt$ISE, mise)
