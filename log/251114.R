mix <- simumix(1)
size <- 1e6

eva <- gtools::rdirichlet(size, rep(1, 3))
eva2 <- gtools::rdirichlet(size, mix$alpha)
grid <- grid_spx(50)
deva <- ddirmix(eva, mix)
deva2 <- ddirmix(eva2, mix)
dgrid <- ddirmix(grid$x, mix)

# plot(density(dgrid), col = 2)
# lines(density(deva), col = "blue")
# lines(density(deva2), col = "darkgreen")

h <- 1.5
n <- 1
data <- rdirmix(n, mix)
dfeva <- ddirkn(eva, data, h)
dfeva2 <- ddirkn(eva2, data, h)
dfgd <- ddirkn(grid$x, data, h)

# plot(density((dfeva - deva)^1), col = "blue")
# lines(density((dfgd - dgrid)^1), col = 2)
# lines(density((dfeva2 - deva2)^1), col = "darkgreen")

# m1 <- MISEf(h, samp_comp(mix), n)
m2 <- mean((dfeva - deva)^2) / gamma(3)
m3 <- mean((dfgd - dgrid)^2) / gamma(3)
m4 <- mean((dfeva2 - deva2)^2 / deva2)
mean(deva2)
# paste0("Theory: ", round(m1, 5))
paste0("Grid: ", round(m3, 5))
paste0("Unif: ", round(m2, 5))
paste0("MC: ", round(m4, 5))
