invisible(lapply(dir("R", full.names = TRUE), source))

# Single h

mix0 <- simumix(2)
# mix0 <- Dir_mix(4)
n <- 5
h <- 0.5

mix1 <- samp_comp(mix0)
(mise1 <- MISEf(h, mix1, n))

mix2 <- samp_mix(mix0, n, B = 1e2, neva = 1e4)
(mise2 <- MISEmc(h, mix2))

mix3 <- grid_mix(mix0, n, split = 50)
(mise3 <- MISEdt(h, mix3))

p <- t.test(attr(mise2, "ISE"), mu = mise1)
p$p.value
p <- t.test(attr(mise3, "ISE"), mu = mise1)
p$p.value

# ==============================================================================
# Vector of h

no <- c(9, 10)
pdf("figure/251113/MISE2.pdf")
for (i in no) {
  mix0 <- simumix(i)
  n <- 10
  (vh <- seq(1, 0.05, len = 10))

  mix1 <- samp_comp(mix0)
  mise1 <- sapply(vh, function(h) {
    cat(h, " ")
    MISEf(h, mix1, n)
  })
  cat("\n")

  mix2 <- samp_mix(mix0, n, neva = 1e6, B = 1e2)
  ise2 <- MISEh(vh, mix2, MISEmc)
  mise2 <- ise2$MISE

  # mix3 <- grid_mix(mix0, n, split = 50)
  # ise3 <- MISEh(vh, mix3, MISEdt)
  # mise3 <- ise3$MISE

  plot.new()
  plot.window(xlim = rev(range(vh)), ylim = range(mise1, mise2, mise3))
  box(); axis(1); axis(2)
  lines(vh, mise1, col = "blue", type = "o")
  lines(vh, mise2, col = "red", type = "o")
  #lines(vh, mise3, col = "darkgreen", type = "o")
  test21 <- as.numeric(MISEtest(ise2$ISE, mise1) > 0.05)
  #test31 <- as.numeric(MISEtest(ise3$ISE, mise1) > 0.05)
  title(main =
          paste0(paste0(test21, collapse = " ")))
        # ,
        #          "\n",
        #          paste0(test31, collapse = " ")))
}
dev.off()

# ==============================================================================
# LSCV
n <- 1e3
mix0 <- simumix(1)
data <- rdirmix(n, mix0)
h <- 0.1
t1 <- Sys.time()
LSCV(h, data)
t2 <- Sys.time()
LSCV1(h, data)
t3 <- Sys.time()
diff(c(t1, t2, t3))

(vh <- seq(1, 0.05, len = 10))

lscv1 <- sapply(vh, function(h) {
  cat(h, " ")
  LSCV(h, data)
  # LSCV1(h, data)
})
lscv2 <- sapply(vh, function(h) {
  cat(h, " ")
  LSCVmc(h, data)
})

plot.new()
plot.window(xlim = rev(range(vh)), ylim = range(lscv1, lscv2))
box(); axis(1); axis(2)
lines(vh, lscv1, col = "blue")
lines(vh, lscv2, col = "red")

# ==============================================================================
# MISE-omit

no <- c(9, 10)
pdf("figure/251113/MISEo2.pdf")
for (i in no) {
  mix0 <- simumix(i)
  n <- 10
  (vh <- seq(1, 0.01, len = 10))

  mixo <- sampo(mix0)
  miseo1 <- sapply(vh, function(h) {
    cat(h, " ")
    MISEfo(h, mixo, n)
  })
  cat("\n")

  mix2 <- samp_mix(mix0, n, neva = 1e6, B = 1e2)
  iseo2 <- MISEh(vh, mix2, MISEmco)
  miseo2 <- iseo2$MISE

  # mix3 <- grid_mix(mix0, n, split = 50)
  # iseo3 <- MISEh(vh, mix3, MISEdto)
  # miseo3 <- iseo3$MISE

  plot.new()
  plot.window(xlim = rev(range(vh)), ylim = range(miseo1, miseo2))
  box(); axis(1); axis(2)
  lines(vh, miseo1, col = "blue")
  lines(vh, miseo2, col = "red")
  # ines(vh, miseo3, col = "darkgreen")
  testo21 <- as.numeric(MISEtest(iseo2$ISE, miseo1) > 0.05)
  # testo31 <- as.numeric(MISEtest(iseo3$ISE, miseo1) > 0.05)
  title(main =
          paste0(paste0(testo21, collapse = " ")))
                 # ,
                 # "\n",
                 # paste0(testo31, collapse = " ")))
}
dev.off()

# ==============================================================================
# Compare
mix0 <- simumix(2)
n <- 1e6
(vh <- seq(1, 0.01, len = 10))

mixa <- samp_comp(mix0)
misea <- sapply(vh, function(h) {
  cat(h, " ")
  MISEf(h, mixa, n)
})

mixo <- sampo(mix0)
miseo <- sapply(vh, function(h) {
  cat(h, " ")
  MISEfo(h, mixo, n)
})

plot.new()
plot.window(xlim = rev(range(vh)), ylim = range(misea, miseo))
box(); axis(1); axis(2)
lines(vh, misea, col = "blue")
lines(vh, miseo, col = "red")
# Ouimet one is better when component is close to the boundary.
