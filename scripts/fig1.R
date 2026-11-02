invisible(lapply(dir("R", full.names = TRUE), source))
pth <- "figures/"

a <- rep(5, 3)
mix <- dirmix(a)
# mix <- Dir_mix(1)

# Projected first derivative
# pdf(paste0(pth, "5.1-projected_gradient.pdf"), height = 5, width = 5)
plot(mix, axl = 300)
vfplot(mix, add = TRUE, axl = 16, scale = 4e-2)
# dev.off()

# First derivative
zoom <- 0.6
V <- matrix(c(0.685, 0.543, -0.486, 0,
              -0.726, 0.558, -0.4, 0,
              0.053, 0.628, 0.778, 0,
              0, 0, 0, 1), nrow = 4)
rgl::open3d(windowRect = c(0, 0, 1500, 1500))
rgl::triangles3d(diag(3), alpha = 0.2)
rgl::view3d(userMatrix = V, zoom = zoom)
vfplot3D(mix, axl = 18, scale = 1.1e-3, add = TRUE, col = "#454545")
# rgl::rgl.snapshot(paste0(pth, "5.1-gradient.png"), fmt = "png")
