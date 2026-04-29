# Organize results from simulation

merge_simu <- function(type = c("dir", "norm"), D, n, r,
                       metrics = c("miae", "mise")) {
  type <- match.arg(type)
  pattern <- paste0(type, "[0-9]+_r", r, "_", D, "D_", n)
  files <- list.files(path = "simu", pattern = pattern, full.names = FALSE)
  fileno <- as.numeric(gsub("_.*", "", basename(files)))
  od <- order(fileno)
  print(files[od])
  loss <- setdiff(1:max(fileno), fileno)
  cat("\nLost in file", min(fileno), "~", max(fileno), ": ")
  print(loss)
  filelist <- lapply(paste0("simu/", files), readRDS)
  mglist <- lapply(metrics, function(m) { # For different metrics
    mats <- lapply(filelist, `[[`, m) # Extract metrics matrix
    cmats <- do.call(rbind, mats) # Combine
    cat("\n", m, ":", dim(cmats))
    as.data.frame(cmats[order(cmats[, "seed"]), ])
  })
  cat("\n")
  mglist <- setNames(mglist, metrics)
  mglist
}

summary_simu <- function(mglist) {
  lapply(mglist, function(mg) {
    mg <- mg[, -c(2, grep("h", colnames(mg))), drop = FALSE]
    res <- aggregate(. ~ no, data = mg, FUN = function(x) {
      m <- mean(x, na.rm = TRUE)
      se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
      c(mean = m, se = se)
    })
    res <- do.call(data.frame, res)
    res[, c("no",
            grep("\\.mean$", names(res), value = TRUE),
            grep("\\.se$", names(res), value = TRUE))]
  })
}

save_simu <- function(filelist, type = c("dir", "norm"),
                      tag, D, n, r, summary = FALSE) {
  type <- match.arg(type)
  sm <- ifelse(summary, "/summary", "")
  path <- paste0("simu/mg_", tag, sm) # Path
  filename <- paste0(type, "_r", r, "_", D, "D_", n) # File name
  pathfile <- paste0(path, "/", filename)
  for (i in names(filelist)) saveRDS(filelist[[i]], file = paste0(pathfile, "_", i))
}

proc_simu <- function(type = c("dir", "norm"), D, n, r,
                      tag, metrics = c("miae", "mise")) {
  type <- match.arg(type)
  res <- merge_simu(type, D, n, r)
  save_simu(res, type, tag, D, n, r)
  res <- summary_simu(res)
  save_simu(res, type, tag, D, n, r, summary = TRUE)
}

# t.test for the simulation results between different methods
testre <- function(re, mthd = NULL) {
  D <- length(mthd)
  n <- length(unique(re[, "no"]))
  pmtx <- matrix(nrow = n, ncol = D)
  colnames(pmtx) <- mthd
  spre <- split.data.frame(re[, mthd], re[, "no"])
  for (i in seq_len(n)) {
    x <- spre[[i]]
    mu <- colMeans(x, na.rm = TRUE)
    idx <- which.min(mu)
    for (j in seq_len(D)) {
      tt <- t.test(x[, idx], x[, j], alternative = "less", paired = TRUE)
      pmtx[i, j] <- tt$p.value
    }
  }
  equal <- pmtx > 0.05
  attr(equal, "p-value") <- pmtx
  equal
}

# Transform simulation results to latex format
letaxre <- function(data, scale = 1, rd = 2, mthd = c("PI", "CV", "LK", "FM")) {
  rd <- paste0("%.", rd, "f")
  data[, -1] <- lapply(data[, -1] * scale, function(x) sprintf(rd, x))
  colname <- as.vector(rbind(paste0(mthd, ".mean"), paste0(mthd, ".se")))
  data <- data[, colname, drop = FALSE]
  mucol <- grepl("mean", names(data))
  secol <- grepl("se", names(data))
  data[mucol] <- lapply(data[mucol], function(x) paste0(" & $", x))
  data[secol] <- lapply(data[secol], function(x) paste0("_{", x, "}$"))
  apply(data, 1, function(x) paste0(x, collapse = ""))
}
