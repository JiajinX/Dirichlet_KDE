invisible(lapply(dir("R", full.names = TRUE), source))
library(caret)

get_acc <- function(est, prior, test) {
  cotest <- closure(test[, 3:8])
  dc <- sapply(est, function(e) dspx(cotest, e))
  likeli <- dc * prior
  post <- closure(likeli)
  est_Channel <- apply(post, 1, which.max)
  mean(est_Channel == test$Channel)
}

wholesale <- read.csv("data/app/Wholesale.csv")
set.seed(1601)
train_idx <- createDataPartition(wholesale$Channel, p = 0.8, list = FALSE)
train <- wholesale[train_idx, ]
cotrain <- closure(train[3:8])
cotrain <- split(cotrain, train$Channel)

fm <- lapply(cotrain, function(cls) optfm(cls))
fmix <- lapply(fm, `[[`, "mix")
# fmix <- list(EMDir(cotrain[[1]], 4)$mix, EMDir(cotrain[[2]], 4)$mix)

interval <- c(0.05, 0.005)
mise <- lapply(1:2, function(i) optMISEh(fmix[[i]], nrow(cotrain[[i]]), interval = interval))
misekn <- lapply(1:2, function(i) dirkn(cotrain[[i]], mise[[i]]$h))

lscv <- lapply(cotrain, function(x) optLSCVh(x, interval = interval))
lscvkn <- lapply(1:2, function(i) dirkn(cotrain[[i]], lscv[[i]]$h))

lcv <- lapply(cotrain, function(x) optLCVh(x, interval = interval))
lcvkn <- lapply(1:2, function(i) dirkn(cotrain[[i]], lcv[[i]]$h))

prior <- c(closure(tabulate(train$Channel)))
test <- wholesale[-train_idx, ]
fm_acc <- get_acc(fmix, prior, test)
mise_acc <- get_acc(misekn, prior, test)
lscv_acc <- get_acc(lscvkn, prior, test)
lcv_acc <- get_acc(lcvkn, prior, test)
c(fm_acc, mise_acc, lscv_acc, lcv_acc)

est <- lscvkn
# Visualization
# \Delta^2
comb <- list(c(1,2), c(3,5), c(4,6)) # , x1 = c(0.3, 1), x2 = c(0.2, 1)
# comb <- list(1, c(3,5), c(2,4,6)) # , x1 = c(0.01, 0.5), x2 = c(0.5, 0.7)
lv <- c(2, 4, 6)
agg_cotrain <- lapply(cotrain, function(x) co_aggregate(x, comb))
points_simplex(agg_cotrain[[1]], col = transcols("#9FB6CD", 0.4))
points_simplex(agg_cotrain[[2]], col = transcols("#9ACD32", 0.4), add = TRUE)
agg_est <- lapply(est, function(e) dir_aggregate(e, comb))
plot(agg_est[[1]], prop = FALSE, lv = lv, scale = prior[1],
     col = transcols("#36648B", 0.9), add = TRUE)
plot(agg_est[[2]], prop = FALSE, lv = lv, scale = prior[2],
     col = transcols("#2E8B57", 0.9), add = TRUE)
decbondary(agg_est[[1]], agg_est[[2]], prior = prior, add = TRUE,
           col = 2, lty = 2, lwd = 2)

# \Delta^3
comb <- list(1, 2, c(3,5), c(4,6))
agg_cotrain <- lapply(cotrain, function(x) co_aggregate(x, comb))
points_tetrahedron(agg_cotrain[[1]], col = "#FFB900")
points_tetrahedron(agg_cotrain[[2]], col = "#5773CC", add = TRUE)
agg_est <- lapply(est, function(e) dir_aggregate(e, comb))
lv <- c(10, 40, 70)
contour3D(agg_est[[1]], lv = lv, col = "#FFB900", add = TRUE)
contour3D(agg_est[[2]], lv = lv, col = "#5773CC", add = TRUE)
decbondary3D(agg_est[[1]], agg_est[[2]], prior = prior,
             color = "grey90", axl = 60, alpha = 1, add = TRUE)

old_par <- par(no.readonly = TRUE)
par(mfrow = c(2, 4))
xall <- 1:6
xcomb <- combn(xall, 2)
for (i in seq_len(ncol(xcomb))) {
  xi <- xcomb[, i]
  comb <- list(setdiff(xall, xi), xi[1], xi[2])
  agg_cotrain <- lapply(cotrain, function(x) co_aggregate(x, comb))
  agg_est <- lapply(est, function(e) dir_aggregate(e, comb))
  comb[[1]] <- paste(comb[[1]], collapse = ",")
  simplexpanel(lab = comb)
  points_simplex(agg_cotrain[[1]], col = transcols("#9FB6CD", 0.4), add = TRUE)
  points_simplex(agg_cotrain[[2]], col = transcols("#9ACD32", 0.4), add = TRUE)
  plot(agg_est[[1]], prop = FALSE, lv = lv, scale = prior[1],
       col = transcols("#36648B", 0.9), add = TRUE)
  plot(agg_est[[2]], prop = FALSE, lv = lv, scale = prior[2],
       col = transcols("#2E8B57", 0.9), add = TRUE)
  decbondary(agg_est[[1]], agg_est[[2]], prior = prior, add = TRUE,
             col = 2, lty = 2, lwd = 2)
}
par(old_par)

# \Delta^1
x <- 1
comb <- list(x, setdiff(1:6, x))
agg_cotrain <- lapply(cotrain, function(x) co_aggregate(x, comb))
agg_est <- lapply(est, function(e) dir_aggregate(e, comb))
betaplot(agg_est[[1]], data = agg_cotrain[[1]][, 1], ymax = 5)
betaplot(agg_est[[2]], data = agg_cotrain[[2]][, 1], add = TRUE, col = 2)

# Histogram of two channel
cosale <- closure(wholesale[3:8])
ggplot(cosale, aes(x = Delicassen, fill = factor(wholesale$Channel))) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30, alpha = 0.5, position = "identity") +
  labs(y = "Density", fill = "Channel")




