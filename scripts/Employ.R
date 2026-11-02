library(WDI)
library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
invisible(lapply(dir("R", full.names = TRUE), source))

# ======================================
# Data clean
emp_ind <- c(
  agriculture = "SL.AGR.EMPL.ZS",
  industry    = "SL.IND.EMPL.ZS",
  services    = "SL.SRV.EMPL.ZS"
)
emp_raw <- WDI(indicator = emp_ind, start = 2025, end = 2025, extra = TRUE)
emp <- emp_raw %>%
  filter(region != "Aggregates") %>%
  select(country, iso2c, iso3c, income, agriculture, industry, services) %>%
  filter(if_all(c(agriculture, industry, services), ~ !is.na(.x))) %>%
  mutate(total = agriculture + industry + services,
         agriculture = agriculture / total,
         industry = industry / total,
         services = services / total) %>%
  select(-total)
saveRDS(emp, "data/Employment")


# ======================================
# Figure 6
emp <- readRDS("data/Employment")
empmtx <- emp %>% select(agriculture, industry, services) %>% as.matrix()
points_simplex(empmtx)
set.seed(4)
fm <- optfm(empmtx, bound = 1.51)
vh <- seq(0.008, 0.004, len = 10)
pih <- MISEfh(vh, fm$mix, n = nrow(empmtx), r = 1, verbose = TRUE)
knpi <- dirkn(empmtx, pih$h, c = 1.5)
plot(knpi)

# Path to the modes
a <- 0.01
ed <- mean_shift(knpi, empmtx, a)
clu <- dist_merge(ed, 0.1)
table(clu$cluster)
mgmap <- distmin_map(clu, empmtx)
clu <- mgcluster(clu, mgmap)
table(clu$cluster)

# Basins of attraction
grid <- simplex_trigrid(7)
ed <- MEM(grid$centers, knpi)
gclu <- getcluster(ed)
grid_mgmap <- distmin_map(gclu, grid$centers)
gclu <- mgcluster(gclu, grid_mgmap)

# Path figure
# pdf("figures/8.2-clustering.pdf", width = 7, height = 6)
cols <- c("#E15759", "#8CD17D", "#A0CBE8", "#F28E2B")
simplexpanel(lab = c("Ag", "In", "Se"))
for (i in seq_len(max(clu$cluster))) {
  st <- empmtx[clu$cluster == i, ]
  points_simplex(st, add = TRUE, cols[i])
  mean_shift(knpi, st, a, pth = TRUE, add = TRUE, col = cols[i])
}
points_simplex(clu$centers, add = TRUE, pch = 3, cex = 1.2, col = 1)

# Attraction figure
cols <- c(cols[c(3, 2, 4, 1)], "grey90")
plot_cluster_region(grid$centers, gclu$cluster, grid$scale, grid$ort, knpi,
                    alpha = c(0.1, 1), gamma = 1, qmax = 0.99,
                    cols = cols, lab = c("Ag", "In", "Se"))
subemp <- c(
  "JP", "MY", "IR", "TR", "SG", "BR", "DZ", "PL", "CZ", "EG", "GR", "CO", "LU",
  "CN", "NG", "PE", "PK", "VN", "CG", "PH", "YE",
  "IN", "BD", "MD", "LR", "AF", "KP",
  "ET", "AO", "TD", "LA", "GQ", "CF", "BI")
subidx <- emp$iso2c %in% subemp
text_simplex(empmtx[subidx, ], emp$iso2c[subidx], col = 1, add = TRUE, cex = 1)
subidx <- emp$iso2c %in% c("US", "TH", "KE", "TZ")
text_simplex(empmtx[subidx, ], emp$iso2c[subidx], col = "blue", add = TRUE, cex = 1)
# dev.off()



# ======================================
# Figure 7
# pdf("figures/8.2-worldmap.pdf", width = 7, height = 6)
clumap <- data.frame(iso3c = emp$iso3c, clust = as.factor(clu$cluster))
world <- ne_countries(scale = "medium", type = "map_units", returnclass = "sf")
plot_clu <- world %>% left_join(clumap, by = c("iso_a3_eh" = "iso3c"))
cols <- c("#E15759", "#8CD17D", "#A0CBE8", "#F28E2B")
# pdf("Map.pdf", width = 7, height = 6)
ggplot(plot_clu) +
  geom_sf(aes(fill = clust), color = "grey70", linewidth = 0.15) +
  scale_fill_manual(values = cols, na.value = "grey90") + # paletteer_d("ggthemes::Tableau_20")
  theme_void() + theme(legend.position = "none")
# dev.off()
