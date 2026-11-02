invisible(lapply(dir("R", full.names = TRUE), source))

# ======================================
# Comparison of Selectors
wholesale <- read.csv("data/Wholesale.csv")
channel <- NA # NA for full; 1 for Horeca; 2 for Retail
if (!is.na(channel)) wholesale <- wholesale[wholesale$Channel == channel, ]
data <- closure(wholesale[, 3:8])
rep <- 100
res <- list()
for (i in seq_len(rep)) res[[i]] <- rlCV(data)
res <- do.call(rbind, res)
# The results are saved in: results/Wholesale/


# ======================================
# Figure 5
full <- closure(wholesale[, 3:8])
horeca <- full[wholesale$Channel == 1, ]
retail <- full[wholesale$Channel == 2, ]
fm <- optfm(full)$mix
h_full <- optMISEh(fm, nrow(full))

# (1,4; 3,5; 2,6)
x3 <- c(2, 6)
x2 <- c(1, 4)
x1 <- setdiff(seq_len(ncol(full)), c(x2, x3))
comb <- list(x1, x2, x3)
mgfu <- co_aggregate(full, comb)
mghc <- co_aggregate(horeca, comb)
mgrt <- co_aggregate(retail, comb)
# pdf("figures/8.1-wholesale.pdf", width = 8, height = 5)
cols <- c("#A0CBE8", "#F28E2B", "grey30")
axl <- 300
lv <- c(30, 60, 85) / 100
lab <- list("Grocery, detergents & paper", "Fresh,\nfrozen", "Milk,\ndelicatessen")
simplexpanel(lab, labcex = 1.2)
points_simplex(mghc, pch = 20, col = cols[1], add = TRUE)
points_simplex(mgrt, pch = 20, col = cols[2], add = TRUE)
plot(dirkn(mgfu, h_full$h), lv = lv, axl = axl,
     col = cols[3], add = TRUE, labcex = 1.2)
legend("topright", c("Horeca", "Retail"), pch = 20, col = cols,
       cex = 1.2, inset = c(0.08, 0))
# dev.off()


