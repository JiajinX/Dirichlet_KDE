# If we know the true mixtures in simulation,
# what is the optimal bandwidth and lowest MISE, MIAE we can get.
invisible(lapply(dir("R", full.names = TRUE), source))

# Dirichlet mixture
## MISE
no <- 6
n <- 500
r <- 0
D <- 3
mix0 <- Dir_mix(no, D)
(pih <- optMISEh(mix0, n, r, verbose = T)) # Best bandwidth in therom

pth <- paste0("data/simu/260404/dir_r", r, "_", D, "D_", n, "_mise")
data <- as.data.frame(readRDS(pth))
plot(density(data[data$no == no, ]$PIh, na.rm = TRUE))
abline(v = pih$h, col = 2)

## LSCV: the second integral has close-form result if we know true Dirichlet mixture
h <- 0.01
x <- rdirmix(n, mix0)
# Density
LSCVf(h, x)
LSCVcf(h, x, mix0)
# First derivative
DLSCVf(h, x)
DLSCVcf(h, x, mix0)


# Normal mixture
no <- 3
D <- 2
mix0 <- Norm_mix(no, D)
n <- 100
r <- 0
h <- 0.1
B <- 100
nsamp <- 1e4
# B training samples and kernels
trainx <- replicate(B, ilrts(rmix(n, mix0))$sx, simplify = FALSE) # Fixed
kns <- lapply(trainx, function(x) dirkn(x, h))
# Generate one common integration sample
samp <- IE_samp(mix0, nsamp) # Fixed
ISE <- numeric(B)
for (i in seq_len(B)) {
  fh <- dspx(samp$x$sx, kns[[i]])
  ISE[i] <- mean((fh - samp$dmix)^2 / samp$dhv)
}
MISEh <- mean(ISE)



# General version
r <- 1
no <- 3
n <- 100
mix0 <- Dir_mix(no)
B <- 1e2
rep <- 10
opth <- matrix(nrow = rep, ncol = 3)
colnames(opth) <- c("Opt_Theory", "Grid_Theory", "Grid_MC")
optmise <- opth
hrange <- 0.003
lh <- 10
for (i in 1:rep) {
  mix0 <- gene_compsamp(mix0)
  mix0 <- gene_data(n, mix0, B) # Generate B training samples and kernels
  mix0 <- gene_eval(mix0, r = r) # Generate evaluated points
  cat(i, " ")
  optreo <- optMISEh(mix0, n, r = r)
  vh <- seq(optreo$h + hrange, optreo$h - hrange, len = lh)
  gdtreo <- MISEfh(vh, mix0, n, r = r)
  gdmc <- MISEh(vh, mix0, MISEmc, r = r)
  opth[i, ] <- c(optreo$h, gdtreo$h, gdmc$h)
  optmise[i, ] <- c(optreo$MISE, gdtreo$MISE, gdmc$MISE)
  if (!gdtreo$opt) opth[i, "Grid_Theory"] <- optmise[i, "Grid_Theory"] <- NA
  if (!gdmc$opt) opth[i, "Grid_MC"] <- optmise[i, "Grid_MC"] <- NA
}
opth <- cbind(no, opth)
optmise <- cbind(no, optmise)
opt <- list(h = opth, mise = optmise)

# =============================================================================

no <- 6
D <- 3
mix0 <- Norm_mix(no, D - 1)
# mix0 <- Dir_mix(no, D)
rep <- 10
opt <- matrix(nrow = rep, ncol = 2)
colnames(opt) <- c("h", "MISE")
r <- 1
nsamp <- 1e3
mix0 <- gene_eval(mix0, r, nsamp)
lh <- 20
vh <- seq(0.009, 0.003, len = lh)
n <- 500
B <- 100
for (i in seq_len(rep)) {
  cat(i, "\n")
  mix0 <- gene_data(n, mix0, B)
  optre <- MISEh(vh, mix0, MISEmc, r = r)
  if (optre$opt) opt[i, ] <- c(optre$h, optre$MISE)
  print(opt)
}
(re <- cbind(no, opt))
beepr::beep(3)
colMeans(re, na.rm = T)



