##map cities
source("01_data_import.R")
source("02_sample_coverage_richness.R")

library(dplyr)
library(sf)
library(ggplot2)
library(patchwork)
library(rnaturalearth)

#load coordinate data
city_coord <- readRDS("data/urban_center_coordinates.rds")

# Keep only cities included in the analysis
cities_map <- city_coord %>%filter(as.character(ID_HDC_G0_2) %in% names(city_list))

# World boundaries
world <- ne_countries(scale = "medium",returnclass = "sf")

#world map
p_world <- ggplot() +
  geom_sf(data = world,fill = "grey94",colour = "grey70",linewidth = 0.2 ) +
  geom_sf(data = cities_map,colour = "steelblue4",size = 1.1,alpha = 0.8) +
  coord_sf(xlim = c(-180, 180),ylim = c(-60, 85),expand = FALSE) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white",colour = NA ) )

#europe inset because too many cities close together
p_europe <- ggplot() +
  geom_sf(data = world,fill = "grey94",colour = "grey65",linewidth = 0.2) +
  geom_sf(data = cities_map,colour = "steelblue4",size = 1.5,alpha = 0.85) +
  coord_sf(xlim = c(-12, 35),ylim = c(35, 72),expand = FALSE) +
  theme_void() +
  theme(panel.border = element_rect(colour = "grey40",fill = NA,linewidth = 0.7),
    panel.background = element_rect(fill = "white",colour = NA) )

#combine
p_map <- p_world +
  inset_element(p_europe,left = 0.03,bottom = 0.04,right = 0.37,top = 0.48 )

p_map #doesn't really help, will just plot full world

p_world

ggsave("out/FS1_global_city_map.pdf",p_world,width = 11,height = 6,device = cairo_pdf)
