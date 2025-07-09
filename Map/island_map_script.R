# Load required libraries
library(ggplot2)
library(sf)
library(viridis)
library(rnaturalearth)
library(rnaturalearthdata)
library(rmapshaper)
library(ggpubr)
library(grid)
library(dplyr)

# Load and simplify world shapefile
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(geounit != "Antarctica") %>%
  ms_simplify()

#Island plotting data from Julian
load("Map\\island_map_plot.RData")

dat <- read.csv("Data\\data_global_SAR.csv")
# dat <- dat[-which(dat$geo_entity == "Greenland"),]
# dat <- filter(dat, suit_geo == 1)
dat <- filter(dat, entity_class != "continent")
# dat <- filter(dat, native_count > 0)
# dat <- filter(dat, area_prop >=0.80)

# Load and transform island point data
islands <- dat %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(st_crs(GADM_islands)) %>%
  rename(`Native richness` = "native_count",
         `Endemic richness` = "endemic_count",
         `Island type` = "category")

# Sort islands by native_count so larger values are plotted last
islands <- islands %>%
  arrange("Native richness")  # ascending order

# Plot
map <- ggplot() +
  geom_sf(data = world, lwd = 0, color = NA) +
  geom_sf(data = GADM_islands, lwd = 0, color = NA) +
  geom_sf(data = islands,
          aes(size = `Native richness`,
              color = `Endemic richness`,
          shape = `Island type`), stroke = 1) +
  coord_sf(crs = "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") +
  scale_x_continuous(breaks = c(0, -90, -180, 90, 180)) +
  scale_y_continuous(breaks = c(0, -23.4, 23.4, 89, -71.3)) +
  scale_color_viridis(
    option = "inferno",
    trans = "log10",
    breaks = c(1, 10, 100, 1000, max(na.omit(islands$`Endemic richness`))),
    limits = c(1, max(na.omit(islands$`Endemic richness`)))
  ) +
  scale_size_continuous(breaks = c(1, 10, 100, 1000, 10000)) +
  theme_minimal() +
  ggtitle("c)") +
  theme(plot.title = element_text(size=22, hjust = 0,
                                          vjust = 2),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),
        legend.title=element_text(size=13))

jpeg("Map.jpeg", width = 30, height = 13, units = "cm",
     res = 300)
map
dev.off()

#save(map, GADM_islands, islands, world, file = "island_map_plot.RData")
