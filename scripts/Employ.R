# Clustering application for the WDI data with employment in different part

library(WDI)
library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
invisible(lapply(dir("R", full.names = TRUE), source))

emp_ind <- c(
  agriculture = "SL.AGR.EMPL.ZS",
  industry    = "SL.IND.EMPL.ZS",
  services    = "SL.SRV.EMPL.ZS"
)

emp <- WDI(indicator = emp_ind, start = 2025, end = 2025, extra = TRUE) %>%
  filter(region != "Aggregates") %>%
  select(country, iso3c, income, agriculture, industry, services) %>%
  filter(if_all(c(agriculture, industry, services), ~ !is.na(.x))) %>%
  mutate(total = agriculture + industry + services,
         agriculture = agriculture / total,
         industry = industry / total,
         services = services / total) %>%
  select(-total)

empmtx <- emp %>% select(agriculture, industry, services) %>% as.matrix()
points_simplex(empmtx)
set.seed(4)
fm <- optfm(empmtx, bound = 1.51)
# fm <- optfm(empmtx)
plot(fm$mix)
vh <- seq(0.008, 0.004, len = 10)
pih <- MISEfh(vh, fm$mix, n = nrow(empmtx), r = 1, verbose = TRUE)
knpi <- dirkn(empmtx, pih$h, c = 1.5)
plot(knpi)
ed <- MEM(empmtx, knpi, pth = TRUE)
clu <- getcluster(ed)
table(clu$cluster)
# mgmap <- dvalley_map(clu, knpi)
mgmap <- distmin_map(clu, empmtx)
clu <- mgcluster(clu, mgmap)
table(clu$cluster)
# pdf("Employment.pdf", width = 7, height = 6)
cols <- c("#F28E2B", "#8CD17D", "#A0CBE8", "#E15759")
plot.meanshift(clu, empmtx, pth = TRUE, obj = knpi, contour = FALSE,
               col = cols, lab = c("Ag", "Se", "In"))

# Mapping color by employment
clumap <- data.frame(iso3c = emp$iso3c, clust = as.factor(clu$cluster))
world <- ne_countries(scale = "medium", type = "map_units", returnclass = "sf")
plot_clu <- world %>% left_join(clumap, by = c("iso_a3_eh" = "iso3c"))
ggplot(plot_clu) +
  geom_sf(aes(fill = clust), color = "grey70", linewidth = 0.15) +
  scale_fill_manual(values = cols, na.value = "grey90") + # paletteer_d("ggthemes::Tableau_20")
  theme_void() + theme(legend.position = "none")
# dev.off()

# Mapping color by income
emp$income[emp$income == "Not classified"] <- NA
incmap <- emp %>%
  mutate(income = factor(income,
                         levels = c("Low income", "Lower middle income",
                                    "Upper middle income", "High income")))

plot_inc <- world %>% left_join(incmap, by = c("iso_a3" = "iso3c"))
ggplot(plot_inc) +
  geom_sf(aes(fill = income), color = "grey70", linewidth = 0.15) +
  scale_fill_manual(values = cols[c(4,1,2,3)], na.value = "grey90") +
  theme_void() + labs(fill = "Income group") + theme(legend.position = "none")
