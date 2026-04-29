invisible(lapply(dir("R", full.names = TRUE), source))

# para <- list(mu = rbind(c(0.1, 0.35, 0.55),
#                         c(0.78, 0.13, 0.09),
#                         c(0.34, 0.14, 0.52)),
#              pr = c(0.24, 0.47, 0.29), h = c(0.08, 0.03, 0.06))
# mix0 <- dirmix(para$mu / para$h + 1, pr = para$pr)
#
# h <- 0.1
# mix <- samp_comp(mix0, mn = 100)
# DMISEf(h, mix, n = 5)
#
# vh <- seq(0.05, 0.01, length.out = 20)
# (dmise <- sapply(vh, DMISEf, mix = mix, n = 5))
# plot(vh, dmise, xlim = rev(range(vh)), type = "o")
#
# mix1 <- simumix(2, c = 1)
# mix2 <- simumix(2, c = 2)
# plot(simumix(2, c = 5))

mix0 <- simumix(2)
# mix0 <- Dir_mix(1)
data <- rdirmix(500, mix0)
h <- 0.1
plot(mix0, prop = FALSE)
vfmixplot(mix0, add = TRUE)
c <- 1
knplot(X = data, h = h, c = c, prop = FALSE)
vfknplot(X = data, h = h, c = c, add = TRUE)

vfknplot(X = data, h = h, c = 1)
vfknplot(X = data, h = h, c = 1.5, add = T, col = "blue")

# X/h + c
mode <- c(0.6, 0.1, 0.3)
h <- 0.1
mix0 <- dirmix(mode/h + 1)
plot(mix0, prop = FALSE, plv = seq(0.9, 0.995, len = 8))
vfmixplot(mix0, add = TRUE)

data <- rdirmix(1000, mix0)
h <- 0.05
c <- 1
knplot(X = data, h = h, c = c, prop = FALSE, plv = seq(0.9, 0.995, len = 8))
vfknplot(X = data, h = h, c = c, add = TRUE)

vfmixplot(mix0)
vfknplot(X = data, h = h, c = 1, add = TRUE, col = "blue")

# MISE with single h ==================================================
c <- 1.5
mix0 <- simumix(1, c = c)
n <- 5

h <- 0.05
mix1 <- samp_comp(mix0, mn = 100)
dmise1 <- DMISEf(h, mix1, n, c)

mix2 <- samp_Dmix(mix0, n, neva = 1e3, B = 1e3)
dmise2 <- DMISEmc(h, c, mix2)

mix3 <- grid_Dmix(mix0, n, split = 50)
dmise3 <- dmise3 <- DMISEdt(h, c, mix3)

plot(density(attr(dmise3, "ISE")))
abline(v = dmise1, col = "red")
abline(v = dmise3, col = "blue")
t.test(attr(dmise3, "ISE"), mu = dmise1)

# MISE with multiple h ==================================================

c <- 1.5
mix0 <- simumix(8, c = c)
# mix0 <- Dir_mix(10, c = c)
n <- 5
vh <- seq(0.2, 0.08, length.out = 10)

mix1 <- samp_comp(mix0)
dmise1 <- sapply(vh, function(i) {
  cat(round(i, 4), " ")
  DMISEf(i, mix1, n, c)
})
plot(vh, dmise1, xlim = rev(range(vh)), type = "o")

mix2 <- samp_Dmix(mix0, n, neva = 1e5, B = 100)
dise2 <- DMISEh(vh, mix2, c, DMISEmc)
dmise2 <- dise2$MISE

mix3 <- grid_Dmix(mix0, n, split = 50)
dise3 <- DMISEh(vh, mix3, c, DMISEdt)
dmise3 <- dise3$MISE

plot.new()
plot.window(xlim = rev(range(vh)), ylim = range(dmise1, dmise2, dmise3))
box(); axis(1); axis(2)
lines(vh, dmise1, type = "o")
lines(vh, dmise2, col = "red", type = "o")
lines(vh, dmise3, col = "blue", type = "o")
(test21 <- as.numeric(MISEtest(dise2$ISE, dmise1) > 0.05))
(test31 <- as.numeric(MISEtest(dise3$ISE, dmise1) > 0.05))
title(main =
        paste0(paste0(test21, collapse = " "),
               "\n",
               paste0(test31, collapse = " ")))
