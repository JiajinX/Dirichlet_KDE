invisible(lapply(dir("R", full.names = TRUE), source))

pdf("figures/7-simumix.pdf", width = 5, height = 5)
no <- 1:6
for (i in no) plot(Dir_mix(i))
for (i in no) plot(Norm_mix(i))
dev.off()
