# Simulation
invisible(lapply(dir("R", full.names = TRUE), source))

scale <- 8e-4
nsamp <- 1e4
# Test mixture
type = "normal"
no <- 5
isdir <- type == "dirichlet"
mix0 <- if (isdir) Dir_mix(no) else Norm_mix(no)
plotdg(mix0, scale = scale)
# Generate data
set.seed(806)
n <- 100
data <- rmix(n, mix0)
if (!isdir) data <- ilrts(data)$sx
# Density estimate
fm <- optfm(data)$mix
if (any(fm$alpha <= 1.5)) stop("Alpha should > 1.5")
plotdg(fm, scale = scale)
pih <- optMISEh(fm, n) # MISE
cvh <- optLSCVh(data) # LSCV
knpi <- dirkn(x = data, h = pih$h)
plot(knpi, prop = FALSE)
kncv <- dirkn(x = data, h = cvh$h)
plot(kncv, prop = FALSE)
samp <- IE_samp(mix0, nsamp)
(IEfm <- IE(samp, fm))
(IEpi <- IE(samp, knpi))
(IEcv <- IE(samp, kncv))
# Partial derivative estimate
c <- 1.5
pidh <- optMISEh(fm, n, r = 1) # MISE
cvdh <- optLSCVh(data, r = 1) # LSCV
kndpi <- dirkn(x = data, h = pidh$h, c = c)
plotdg(kndpi, scale = scale)
kndcv <- dirkn(x = data, h = cvdh$h, c = c)
plotdg(kndcv, scale = scale)
samp <- DIE_samp(mix0, nsamp)
(IEdfm <- DIE(samp, fm))
(IEdpi <- DIE(samp, kndpi))
(IEdcv <- DIE(samp, kndcv))


# Simulation part
n <- 500
rep <- 5
re <- msdst(type = "dir", no = no, n = n, rep = rep)
dre <- msfst(type = "dir", no = no, n = n, rep = rep)
# no = 1, n = 500, i = 4
# no = 4, n = 500, i = 5
# how about nsp = 1e3??
