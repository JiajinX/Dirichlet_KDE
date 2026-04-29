B <- 50

# Evaluate f^2
ps0 <- numeric(B)
for (i in 1:B) {
  sp <- rdirmix(1e3, mix0)
  ps0[i] <- mean(ddirmix(sp, mix0))
}
(mps0 <- mean(ps0))

# Evaluate hat(f)f
sp <- rdirmix(1e3, mix0)
ps1 <- numeric(B)
for (i in 1:B) {
  data <- rdirmix(500, mix0)
  f <- dirmix(data / h + 1)
  ps1[i] <- mean(ddirmix(sp, f))
}
(mps1 <- mean(ps1))

# Evaluate hat(f)^2
ps2 <- numeric(B)
for (i in 1:B) {
  data <- rdirmix(500, mix0)
  f <- dirmix(data / h + 1)
  sp <- rdirmix(1e3, f)
  ps2[i] <- mean(ddirmix(sp, f))
}
(mps2 <- mean(ps2))

mps2 + mps0 - 2*mps1

# ----------------------------------------

y1 <- attr(mix, "samp")
Cy1 <- Cdir(y1/h + 1, log = TRUE)

Cy1y1 <- Cdir((y1 + y1)/h + 1, log = TRUE)
mean(exp(Cy1 + Cy1 - Cy1y1))

y2 <- samp_comp(mix0)
y2 <- attr(y2, "samp")
Cy2 <- Cdir(y2/h + 1, log = TRUE)
Cy1y2 <- Cdir((y1 + y2)/h + 1, log = TRUE)
(psi2 <- mean(exp(Cy1 + Cy2 - Cy1y2)))


Cpls <- Cdir_plus(y1/h + 1, y2/h + 1, log = TRUE)
dim(Cpls)
