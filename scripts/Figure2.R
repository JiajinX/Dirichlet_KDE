invisible(lapply(dir("R", full.names = TRUE), source))
pth <- "figure/"

axld <- 300
axlg <- 16
scale <- 4e-2

pdf(paste0(pth, "simumix.pdf"), height = 5, width = 5)
for (i in 1:6) plotdg(Dir_mix(i), scale = scale, axld = axld, axlg = axlg)
for (i in 1:6) plotdg(Norm_mix(i), scale = scale, axld = axld, axlg = axlg)
dev.off()
