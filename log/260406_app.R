# Explore real data
invisible(lapply(dir("R", full.names = TRUE), source))
library(dplyr)
data(Hydrochem)
hydro <- Hydrochem %>% select(Location, River, Ca, K, Na, Mg)
data <- hydro %>% select(Ca, Na, Mg)
data <- as.matrix(data / rowSums(data))

points_simplex(data[hydro$River == "Anoia", ], col = "blue")
points_simplex(data[hydro$River == "Cardener", ], col = "red", add = T)
points_simplex(data[hydro$River == "UpperLLobregat", ], col = "darkgreen", add = T)
points_simplex(data[hydro$River == "LowerLLobregat", ], col = "grey30", add = T)

points_simplex(data[hydro$Location == "LowerLLobregat", ], col = "blue")
points_simplex(data[hydro$Location == "UpperLLobregat", ], col = "darkgreen", add = T)

# =================
library(loon)
data("olive")
olive <- olive %>% select(-palmitic)
data <- olive %>% select(palmitoleic, linolenic, arachidic)
data <- as.matrix(data / rowSums(data))

points_simplex(data[olive$Region == "South", ], col = "blue")
points_simplex(data[olive$Region == "Sardinia", ], col = "red", add = T)

points_simplex(data[olive$Area == "South-Apulia", ], col = "blue")
points_simplex(data[olive$Area == "Inland-Sardinia", ], col = "red", add = T)
# points_simplex(data[olive$Area == "Calabria", ], col = "darkgreen", add = T)
points_simplex(data[olive$Area == "North-Apulia", ], col = "orange", add = T)

# ====================
load("data/app/BDDSegX.RData")
BDD <- BDDSegX[, 7:11]
BDD[, "S_ADE"] <- BDD[, "S_A"] + BDD[, "S_D"] + BDD[, "S_E"]
coBDD <- as.matrix(BDD[, c("S_ADE", "S_B", "S_C")])

# ====================
wc <- read.csv("data/app/Wholesale customers data.csv")
cowc <- as.matrix(wc[, -c(1,2)])
cowc <- cowc / rowSums(cowc)

x1 <- "Fresh"
x2 <- "Grocery"
other <-  1 - cowc[, x1] - cowc[, x2]
subcowc <- cbind(cowc[, x1], cowc[, x2], other)
points_simplex(subcowc[wc$Channel == 1, ], col = "red")
points_simplex(subcowc[wc$Channel == 2, ], col = "blue", add = TRUE)

fm <- optfm(cowc, bound = 1.51)
pih <- optMISEh(fm$mix, nrow(cowc), r = 1)
knpi <- dirkn(cowc, pih$h, c = 1.5)

# ====================
library("DirichletReg")

# Rocks
rk <- Rocks
cork <- Rocks[, 1:5]
cork <- cork / rowSums(cork)
par_repoint(cork, group = rk$type)
fm <- optfm(cork)
pih <- optMISEh(fm$mix, nrow(cork), 1, verbose = TRUE)
knpi <- dirkn(cork, pih$h, c = 1.5)
par_reden(knpi)

GlacialTills

# ====================
library("robCompositions")
n <- nrow(ageCatWorld)
D <- 3
noise <- exp(rnorm(n * D, sd = 0.04))
range(noise)
ageprop <- ageCatWorld[, 1:3] * noise # Add some perturbations
ageprop <- ageprop / rowSums(ageprop)
points_simplex(ageprop)
fm <- optfm(ageprop)
plot(fm$mix)
pih <- optMISEh(fm$mix, nrow(ageprop), r = 1,
                interval = c(5e-3, 5e-4), tol = 1e-4, verbose = TRUE)
knpi <- dirkn(ageprop, pih$h, c = 1.5)
plot(fm$mix, prop = F, plv = seq(0.8, 1, len = 30))
plot(knpi, prop = F, plv = seq(0.8, 1, len = 10))
# Step is too large: row 8, 9, 16
# Error: 87, 140, 184 -> go to the boundary
ed <- findmode(ageprop, knpi, step_init = 1e-3, plot = 2, add = T)
clu <- getcluster(ed)
plot.meanshift(clu, ageprop, pth = T, obj = knpi, steplen = 1e-3)

pidh <- optMISEh(fm$mix, nrow(ageprop),
                 interval = c(5e-3, 5e-4), tol = 1e-4, verbose = T)
plot(dirkn(ageprop, pidh$h), prop = F, plv = seq(0.8, 1, len = 10))
