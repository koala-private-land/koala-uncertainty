library(dplyr)
library(ggplot2)
library(sf)
library(scales)
library(cowplot)
library(ggrepel)
library(ggspatial)

load('plots/shapes.RData')
town_centres <- read_sf('data/UCL_2021_AUST_GDA2020_SHP/UCL_2021_AUST_GDA2020.shp')

aus_plot <- ggplot() +
  #geom_sf(data = aus_border_union, fill = NA) +
  geom_sf(data = aus_border, fill = 'gray80') +
  geom_sf(data = nsw_lga_union, size = 0.5, fill = '#FFD580') +
  geom_rect(aes(xmin = nsw_bbox$xmin - bbox_buffer, 
                xmax = nsw_bbox$xmax + bbox_buffer, 
                ymin = nsw_bbox$ymin - bbox_buffer, 
                ymax = nsw_bbox$ymax + bbox_buffer), 
            fill = NA, color = 'gray50') +
  coord_sf(xlim = aus_xlim, ylim = aus_ylim) +
  theme_void()+
  theme(panel.background = element_rect(fill = 'white', colour = 'gray50'))

town_centre_df <- town_centres %>% 
  dplyr::filter(UCL_NAME21 %in% c('Sydney', 'Newcastle', 'Wollongong', 'Coffs Harbour', 'Port Macquarie', 'Armidale', 'Tamworth', 'Dubbo', 'Bathurst', 'Albury')) %>% 
  st_centroid() %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2])

main_map <- ggplot() +
  geom_sf(data = aus_border) +
  geom_sf(data = nsw_lga, size = 0.5, fill = 'antiquewhite1', color = 'gray50') +
  geom_sf(data = nsw_lga_union, size = 1, color = '#FFD580', fill = NA) +
  geom_sf(data = town_centre_df) +
  geom_text_repel(data = town_centre_df, aes(x = lon, y = lat, label = UCL_NAME21)) +
  coord_sf(xlim = c(nsw_bbox$xmin, nsw_bbox$xmax+3), ylim = c(nsw_bbox$ymin, nsw_bbox$ymax)) +
  theme_void()+
  ggspatial::annotation_scale(location = "tr", width_hint = 0.4) +
  ggspatial::annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering)+
  ggtitle("Map of study area") +
  theme(panel.background = element_rect(fill = 'aliceblue'),
        panel.border = element_rect(colour = 'gray50', linewidth = 1))


final_plot <- main_map %>%
  ggdraw() +
  draw_plot(aus_plot,
            x = .6, 
            y = 0,
            width = .33, 
            height = .33)

ggsave(final_plot, filename = 'plots/nsw_map.png')

