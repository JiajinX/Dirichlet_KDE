# In the Gaussian mixture case, the optimal bandwidth for the first derivative can be
# smaller than that for the density because the true derivative is highly concentrated
# near the boundary of the simplex. The logistic-normal density may have steep boundary
# behavior after transformation back to the simplex. A larger bandwidth over smooths
# this steep derivative and produces substantial boundary bias. Therefore, the
# derivative-optimal bandwidth may be shifted toward a smaller value.

mix0 <- Norm_mix(1)
# mix0 <- Dir_mix(6)
r <- 1
nsamp <- 1e4
mix0 <- gene_eval(mix0, r, nsamp)
eva <- attr(mix0, "evapoint")$eva
deva <- attr(mix0, "evapoint")$deva
smeva <- apply(eva, 1, min)
# gn <- sqrt(rowSums(deva^2)) # Length of vector
# summary(smeva)
# summary(gn)
# plot(smeva, gn, xlab = "min component of x", ylab = "norm of true derivative")

gn2 <- rowSums(deva^2)
sum(gn2[smeva < 0.1]) / sum(gn2)

# Gaussian
# 1. 0.78  2. 0.87  3. 0.74
# 4. 0.76  5. 0.90  6. 0.57
# Dirichlet
# 1. 0.39  2. 0.19  3. 0.62
# 4. 0.64  5. 0.62  6. 0.62

# thd <- seq(1e-3, 3e-1, by = 5e-3)
# sumgn2 <- sum(gn2)
# prop <- sapply(thd, function(i) sum(gn2[smeva < i])) / sum(gn2)
# plot(thd, prop, type = "l")
# names(prop) <- thd
# round(prop, 4)

check_boundary <- function(h, mix, c = 1.5) {
  eva <- attr(mix, "evapoint")$eva
  deva <- attr(mix, "evapoint")$deva
  data <- attr(mix, "data")
  smeva <- apply(eva, 1, min)
  B <- length(data)
  prop <- rep(NA, B)
  for (i in seq_len(B)) {
    # print(i)
    kn <- dirkn(data[[i]], h, c)
    df <- pddir(eva, kn)
    err2 <- rowSums((df - deva)^2)
    prop[i] <- sum(err2[smeva < 0.1]) / sum(err2)
  }
  mean(prop)
}

mix0 <- Norm_mix(1)
mix0 <- Dir_mix(4)
r <- 1
nsamp <- 1e3
mix0 <- gene_eval(mix0, r, nsamp)
mix0 <- gene_data(n = 100, mix0, B = 100)
check_boundary(0.0127, mix0)
check_boundary(0.0125, mix0)
