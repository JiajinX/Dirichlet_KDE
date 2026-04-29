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

emp <- WDI(indicator = emp_ind, start = 2023, end = 2023, extra = TRUE) %>%
  filter(region != "Aggregates") %>%
  select(country, iso3c, income, agriculture, industry, services) %>%
  filter(if_all(c(agriculture, industry, services), ~ !is.na(.x))) %>%
  mutate(total = agriculture + industry + services,
         agriculture = agriculture / total,
         industry = industry / total,
         services = services / total) %>% select(-total)

empmtx <- emp %>% select(agriculture, industry, services) %>% as.matrix()
points_simplex(empmtx)
fm <- optfm(empmtx, bound = 1.51)
plot(fm$mix)
pih <- optMISEh(fm$mix, nrow(empmtx), r = 1, interval = c(0.008, 0.004), tol = 1e-4, verbose = TRUE)
knpi <- dirkn(empmtx, pih$h, c = 1.5)
plot(knpi)
ed <- findmode(empmtx, knpi, step_init = 1e-3, plot = 2, add = T)
clu <- getcluster(ed)
cols <- c("#E15759", "#F28E2B", "grey90",  "#8CD17D", "#A0CBE8")
plot.meanshift(clu, empmtx, pth = F, obj = knpi, steplen = 1e-3, col = cols)
table(clu$cluster)

# Mapping color by employment
clumap <- data.frame(iso3c = emp$iso3c, clust = as.factor(clu$cluster))
world <- ne_countries(scale = "medium", type = "map_units", returnclass = "sf")
plot_clu <- world %>% left_join(clumap, by = c("iso_a3" = "iso3c"))
ggplot(plot_clu) +
  geom_sf(aes(fill = clust), color = "grey70", linewidth = 0.15) +
  scale_fill_manual(values = cols, na.value = "grey90") + # paletteer_d("ggthemes::Tableau_20")
  theme_void() + theme(legend.position = "none")

# Mapping color by income
incmap <- emp %>%
  mutate(income = factor(income,
                         levels = c("Not classified", "Low income", "Lower middle income",
                                    "Upper middle income", "High income")))
plot_inc <- world %>% left_join(incmap, by = c("iso_a3" = "iso3c"))
ggplot(plot_inc) +
  geom_sf(aes(fill = income), color = "grey70", linewidth = 0.15) +
  scale_fill_manual(values = cols[c(3, 1, 2, 4, 5)], na.value = "grey90") +
  theme_void() + labs(fill = "Income group") + theme(legend.position = "none")





