
rm(list=ls())

#Install and load required packages
packages <- c('ggplot2', 'sf', 'units', 'rnaturalearth')
install.packages(setdiff(packages, rownames(installed.packages())))
library(ggplot2)
library(units)
library(sf)
library(rnaturalearth)
library(ggpubr)

# Load world data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Create your dataframe
df <- data.frame(
  labels = c('Madonna', 'Muir',
             'Crocker', 'Cold Creek',
             'Bolinger', 'Gunstock',
             'Fort Columbia State Park', 'Gravel Pit',
             'Tenmile Creek', 'Lake Sylvia State Park',
             'Scott Lake National Forest', 'Wildboy Creek'),
  latitude = c(37.003992, 37.862, 
          38.770492, 34.09232115,
          37.812139, 37.836388,
          46.256075, 44.640211,
          44.10186, 47.002435,
          44.21835, 45.673955516464),
  longitude = c(-121.681136, -122.573239, 
           -122.970651, -118.6477477,
           -122.044156, -122.54455,
           -123.917127, -123.28802,
           -124.06072, -123.589958,
           -121.89018, -122.21913495478067),
  species = c('Taricha torosa', 'Taricha torosa',
              'Taricha torosa', 'Taricha torosa',
              'Taricha torosa', 'Taricha torosa',
              'Taricha granulosa', 'Taricha granulosa',
              'Taricha granulosa', 'Taricha granulosa',
              'Taricha granulosa', 'Taricha granulosa')
  )


tt_data <- sf_df[sf_df$species == 'Taricha torosa', ]
tg_data <- sf_df[sf_df$species == 'Taricha granulosa', ]

# Convert to sf object
sf_df <- st_as_sf(df, 
                  coords = c('longitude', 
                             'latitude'), 
                  crs = 4326)



# Plot the world map
plot(world['geometry'])

# Plot the data points
plot(sf_df, pch = 16, col = "red", add = TRUE)


northam <- ne_states(country = 'united states of america', 
                     returnclass = 'sf')
states <- northam[northam$name %in% c('California', 
                                      'Oregon', 
                                      'Washington'),]
ggplot() +
  geom_sf(data = states, 
          fill = 'lightgray', 
          color = 'black') +
  geom_sf(data = states, 
          fill = 'white', 
          color = 'black') +
  geom_sf(data = sf_df, 
          color = 'red', 
          size = 3)+
  theme_minimal()
  
bbox_sf_bay_area <- c(-122.6, 37.3, -121.7, 37.9)
bbox_sf_bay_area_tight <- c(-122.7, 37.5, -121.9, 38.1)  # xmin, ymin, xmax, ymax
corner_coords <- as.data.frame(expand.grid(x = c(-122.7, -121.9), y = c(37.5, 38.1)))


df$shape_factor <- as.factor(df$labels)

# Assign unique integers to each label in the factor
levels(df$shape_factor) <- seq_along(levels(df$shape_factor))


full_map <- ggplot() +
  geom_sf(data = states, fill = 'lightgray', color = 'black') +
  geom_sf(data = states, fill = 'white', color = 'black') +
  geom_sf(data = tt_data, aes(color = "Taricha torosa", 
                              shape = factor(shape)), size = 3) +
  geom_sf(data = tg_data, aes(color = "Taricha granulosa", 
                              shape = factor(shape)), size = 3) +
  scale_shape_manual(values = c("14" = 14, "15" = 15, "16" = 16, 
                                "17" = 17, "18" = 18, "19" = 19)) +
  scale_color_manual(values = c("Taricha torosa" = 'red', "Taricha granulosa" = 'blue')) +
  guides(color = guide_legend(title = "Species", 
                              ncol = 2, 
                              keywidth = 1, 
                              keyheight = 1,
                              override.aes = list(shape = c(16, 16))))
full_map



half_map
ggarrange(half_map, full_map, ncol = 2, widths = c(1,2))
plot_main_with_lines <- plot_main +
  geom_segment(data = corner_coords,
               aes(x = x, y = y, xend = bbox_sf_bay_area[c(1, 3)], yend = bbox_sf_bay_area[c(2, 4)]),
               color = "red",
               size = 1)

plot_main_with_lines

# Define the bounding box for the San Francisco Bay Area
bbox_sf_bay_area <- c(-122.6, 37.3, -121.7, 37.9)  # xmin, ymin, xmax, ymax

# Plot the overall map
full_map <- ggplot() +
  geom_sf(data = states, fill = 'lightgray', color = 'black') +
  geom_sf(data = states, fill = 'white', color = 'black') +
  geom_sf(data = tt_data, aes(color = "Taricha torosa"), size = 3) +
  geom_sf(data = tg_data, aes(color = "Taricha granulosa"), size = 3) +
  scale_color_manual(values = c("Taricha torosa" = 'red', "Taricha granulosa" = 'blue')) +
  guides(color = guide_legend(title = "Species", 
                              ncol = 1, 
                              keywidth = 1, 
                              keyheight = 1,
                              override.aes = list(shape = c(16, 16)),
                              title.theme = element_text(size = 30,
                                                         hjust = 0.5),
                              label.theme = element_text(face = 'italic',
                                                         size = 25))) +
  theme(legend.text = element_text(face = 'italic',
                                   size = 25),
        legend.title = element_text(face = 'bold',
                                    size = 20))+ 
  geom_rect(
    aes(xmin = bbox_sf_bay_area[1] - 0.1, xmax = bbox_sf_bay_area[3] + 0.1,
        ymin = bbox_sf_bay_area[2] - 0.1, ymax = bbox_sf_bay_area[4] + 0.1),
    fill = 'white', color = "black"  # Add a black outline
  )+
  theme_minimal()
full_map

# Plot the zoomed area
  half_map <- ggplot() +
    geom_sf(data = states, fill = 'lightgray', color = 'black') +
    geom_sf(data = states, fill = 'white', color = 'black') +
    geom_sf(data = tt_data, aes(color = "Taricha torosa"), size = 3) +
    geom_sf(data = tg_data, aes(color = "Taricha granulosa"), size = 3) +
    scale_color_manual(values = c("Taricha torosa" = 'red', "Taricha granulosa" = 'blue')) +
    guides(color = FALSE) +
    coord_sf(xlim = c(-122.7, -121.9), ylim = c(37.5, 38.1)) +  # Adjust the coordinate limits
    theme_minimal()+
    theme(panel.border = element_rect(color = 'black',
                                      size = 1,
                                      fill = NA))
half_map



full_map_w_lines <- full_map +
  geom_segment(
    x = c(-122.6, -121.7), y = c(37.9, 37.9),  # Coordinates for the lines
    xend = c(-122.6, -122.6), yend = c(37.9, 37.3),  # End coordinates for the lines
    color = "black", size = 0.5  # Line color and size
  ) +
  geom_segment(
    x = c(-121.7, -121.7), y = c(37.9, 37.3),  # Coordinates for the lines
    xend = c(-122.6, -121.7), yend = c(37.3, 37.3),  # End coordinates for the lines
    color = "black", size = 0.5  # Line color and size
  )




ggarrange(half_map, full_map, ncol = 2, widths = c(1,3))

# Arrange plots side by side
plot_main / plot_zoom
