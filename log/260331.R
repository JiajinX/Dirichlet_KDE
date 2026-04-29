invisible(lapply(dir("R", full.names = TRUE), source))

pth <- "data/simu/260404"
summary <- T
summary <- ifelse(summary, "/summary/", "/")

r <- 0
D <- 3
n <- 500
type <- "norm"
ms <- "mise"

doc <- paste0(pth, summary, type, "_r", r, "_", D, "D_", n, "_", ms)
(data <- readRDS(doc))
data <- data[rowSums(is.na(data)) == 0, ]

# Bandwidth Diff ===========================
diff <- data[, 5] - data[, 6]
drop <- which(abs(diff) > 0.5)
diff1 <- diff[-drop]
plot(density(diff1))

library(dplyr)
data1 <-  data[-drop, ]
data1 <- as.data.frame(data1)
data1 %>% group_by(no) %>% summarise(mean(PI, na.rm = TRUE),
                                     mean(CV, na.rm = TRUE))

# Loss Diff =================================
diff <- data[, 3] - data[, 4]
data[which.max(diff),]

