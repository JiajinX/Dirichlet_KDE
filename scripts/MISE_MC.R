invisible(lapply(dir("R", full.names = TRUE), source))

vh <- seq(0.02, 0.005, len = 31)
no <- 2
mix0 <- Dir_mix(no) # Test mixture
n <- 500 # Sample size
r <- 0 # 0,1: density / first derivative

# MISE from the direct Monte Carlo approximation
Q <- 1e4
B <- 1e3
mix0 <- gene_eval(mix0, r, Q)
mix0 <- gene_data(n, mix0, B)
dmise <- rep(NA, length(vh))
names(mise) <- vh
for (i in seq_along(vh)) {
  h <- vh[i]
  dmise[i] <- MISEmc(h, mix0, r)
}

# MISE from the exact MISE representation
S <- 1e6
mix0 <- gene_compsamp(mix0, mn = S)
emise <- MISEfh(vh, mix0, n, r, exhaustive = TRUE)$info

# Results are saved in: results/MISE_MC/
