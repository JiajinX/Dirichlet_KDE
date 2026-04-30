invisible(lapply(dir("R", full.names = TRUE), source))
library(dplyr)
library(ggplot2)
library(patchwork)

# Close form of left-hand side
DLSCVcf <- function(h, x, mix, c = 1.5) {
  n <- nrow(x)
  D <- ncol(x)
  ya <- x/h + c
  Cx <- Cdir(ya, log = TRUE)
  outCxx <- outer(Cx, Cx, "+")
  Cxx <- Cdir_plus(ya, ya, log = TRUE)
  Cd <- exp(outCxx - Cxx)

  # First term
  idx1 <- unlist(lapply(1:n, function(i) i:n))
  idx2 <- rep(1:n, n:1)
  diag <- idx1 == idx2
  # Indices of lower triangular elements in symmetric matrix
  log_ya_1 <- log(ya - 1)
  Csup <- log_ya_1[idx1, , drop = FALSE] + log_ya_1[idx2, , drop = FALSE]
  sumya <- ya[idx1, , drop = FALSE] + ya[idx2, , drop = FALSE]
  log_sumya_2 <- log(sumya - 2)
  log_sumya_3 <- log(sumya - 3)
  Cs1 <- rowSums(exp(Csup - (log_sumya_2 + log_sumya_3)))
  Cs2 <- rowSums(exp(Csup - 2 * log_sumya_2))
  tx <- rowSums(exp(log_ya_1[idx1, , drop = FALSE] - log_sumya_2))
  ty <- rowSums(exp(log_ya_1[idx2, , drop = FALSE] - log_sumya_2))
  Cscut <- tx * ty / D
  Cs <- (1 - 1/D) * Cs1 + Cs2/D - Cscut
  r1 <- Cs * Cd[lower.tri(Cd, diag = TRUE)]
  r1 <- 2 * sum(r1[!diag]) + sum(r1[diag])
  logC <- prefactor(1/h + c * D, 1/h + c * D, D)
  term1 <- r1 * exp(logC) / n^2

  # Second term
  a0 <- rowSums(mix$alpha)
  outCay <- outer(Cx, attr(mix$alpha, "logmbeta"), "+")
  Cay <- Cdir_plus(ya, mix$alpha, log = TRUE)
  Cd <- exp(outCay - Cay)
  Cs <- D1S(ya, mix$alpha)
  logC <- prefactor(1/h + c * D, a0, D)
  r <- Cd * Cs * exp(logC[rep(1, n), ])
  term2 <- sum(mix$pr * t(r)) / n
  c(term1 - 2 * term2)
}

# Plot LSCV lines for different bandwidth
hLSCV <- function(data, rg = NULL) {
  if (is.null(rg)) rg <- range(data$vh)
  minmc <- data$vh[which.min(data$mc)]
  mincls <- data$vh[which.min(data$cls)]
  print(paste0("Minh ", "[close: ", round(mincls, 3),
               " | ", "mc: ", round(minmc, 3), "]"))
  sub <- subset(data, vh <= rg[2] & vh >= rg[1])
  P1 <- ggplot(sub, aes(x = vh)) +
    geom_line(aes(y = mc, color = "mc")) +
    geom_line(aes(y = cls, color = "cls")) +
    geom_vline(aes(xintercept = minmc, color = "mc"), linetype = "dashed") +
    geom_vline(aes(xintercept = mincls, color = "cls"), linetype = "dashed") +
    scale_x_reverse(limits = range(sub$vh)) +
    scale_color_manual(name = "Method",
                       values = c("mc" = "red", "cls" = "black"),
                       labels = c("mc" = "MC", "cls" = "Close")) +
    labs(x = "h", y = "LSCV")
  P2 <- ggplot(sub, aes(x = vh)) +
    geom_line(aes(y = cls - mc)) +
    scale_x_reverse() +
    labs(x = NULL, y = "Close - MC") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  print(P2 / P1 + plot_layout(heights = c(4, 6)))
  sub
}

n = 500
mix <- Dir_mix(1)
seed <- 206 # [1: 205, 206] [3: 5609]
set.seed(seed) # 206
data <- rdirmix(n, mix)
CVh <- optLSCVh(data, 1)
(opth <- CVh$h) # Optimal bandwidth h from LSCV

# Test for different h value
vh <- c(seq(0.5, 0.1, length = 10),
        seq(0.09, 0.01, length = 25),
        seq(0.009, 0.001, length = 25))
hle <- length(vh)
cls <- mc <- rep(NA, rep = hle)
for (i in 1:hle) {
  print(i)
  ih <- vh[i]
  mc[i] <- DLSCVf(ih, data)
  cls[i] <- DLSCVcf(ih, data, mix)
}
vhre <- as.data.frame(cbind(seed, vh, cls, mc))
hLSCV(vhre, rg = c(0.001, 0.5))



