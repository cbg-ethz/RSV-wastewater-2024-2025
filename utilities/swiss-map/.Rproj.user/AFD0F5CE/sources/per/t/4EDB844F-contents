library(sf)
library(ggplot2)
library(ggspatial)
library(tidyverse)
library(data.table)
library(patchwork)
library(scatterpie)
library(dplyr)
library(ggforce)
library(tidyr)
library(scales)

source('R/load_spatial_data.R')

rsv_ab_assay <- read.csv(file = '../../RSV_results/data_analysis/data/concentrations/rsv_ab_assay/251023_RSVAB.csv',
                         check.names = FALSE) %>%
  select(-c(1)) %>%
  rename(Date = collection_date) %>%
  rename(location = ara_name) %>%
  select(Date, pathogen, location, load)


conc_subtypes <- rsv_ab_assay %>%
  
  pivot_wider(
    names_from = pathogen,
    values_from = load
  )

conc_subtypes$location <- gsub("Werdhölzli", "Zurich", conc_subtypes$location)
conc_subtypes$location <- gsub("Aire", "Geneva", conc_subtypes$location)
conc_subtypes$location <- gsub("ARA Basel", "Basel", conc_subtypes$location)
conc_subtypes$location <- gsub("CDA Lugano", "Lugano", conc_subtypes$location)
conc_subtypes$location <- gsub("ARA Sensetal", "Laupen", conc_subtypes$location)
conc_subtypes$location <- gsub("ARA Chur", "Chur", conc_subtypes$location)

conc_subtypes$Date <- as.character(as.Date(as.character(conc_subtypes$Date),"%Y-%m-%d"))
conc_subtypes <- conc_subtypes[conc_subtypes$location %in% c("Zurich", "Geneva", "Basel", "Lugano", "Laupen", "Chur"), ]

rsv_ab_assay <- conc_subtypes %>% 
  rename("RSV_A" = "RSV-A") %>%
  rename("RSV_B" = "RSV-B")


rsv_ab_assay_summary <- rsv_ab_assay %>% 
  mutate(
    RSV_A_prop = RSV_A / (RSV_A + RSV_B),
    RSV_B_prop = RSV_B / (RSV_A + RSV_B)
    
  ) %>%
  group_by(location) %>%
  summarise(
    RSV_A_median_prop = median(RSV_A_prop, na.rm = TRUE),
    RSV_B_median_prop = median(RSV_B_prop, na.rm = TRUE),
    
    RSV_A_median_prop_IQR = IQR(RSV_A_prop, na.rm = TRUE),
    RSV_B_median_prop_IQR = IQR(RSV_B_prop, na.rm = TRUE)
  )
rsv_ab_assay_summary$ID <- c(270101, 390101, 664301, 66700, 515100, 26101)
rsv_ab_assay_summary$total_population <- c(262917, 54729, 405083, 74363, 117963, 475198)


catchments = load_catchments()
catchments = st_transform(catchments, 25830)
catchment_centroids = st_centroid(catchments)

plz_locs = load_plz_centroids(crs_required = 25830)

canton_boundaries = st_read('data/spatial/swissboundaries3d_2023-01_2056_5728.shp/swissBOUNDARIES3D_1_4_TLM_KANTONSGEBIET.shp')

aras_to_plot = rsv_ab_assay_summary$ID 

cat_locs_sites = catchment_centroids %>% filter(ara_id %in% aras_to_plot)

cat_locs_sites = cat_locs_sites %>% rename(ID = ara_id) 


cat_locs_sites_data = merge(cat_locs_sites, rsv_ab_assay_summary, by= c('ID'))

catch_poly_sites = catchments %>% filter(ara_id %in% aras_to_plot) 


cat_locs_sites_data$RSV_A_prop <- as.numeric(sapply(cat_locs_sites_data$RSV_A_median_prop, `[`, 1))
cat_locs_sites_data$RSV_B_prop <- as.numeric(sapply(cat_locs_sites_data$RSV_B_median_prop, `[`, 1))
cat_locs_sites_data$total <- as.numeric(cat_locs_sites_data$total)

# Extract numeric coordinates from geometry
coords <- st_coordinates(cat_locs_sites_data$geometry)[, 1:2]
cat_locs_sites_data <- cat_locs_sites_data %>%
  mutate(x = coords[,1], y = coords[,2])



# -------------------------------
# 1. TRANSFORM SPATIAL DATA TO PROJECTED CRS (meters)
# -------------------------------
canton_boundaries   <- st_transform(canton_boundaries, 2056)
catch_poly_sites    <- st_transform(catch_poly_sites, 2056)
cat_locs_sites_data <- st_transform(cat_locs_sites_data, 2056)

# -------------------------------
# 2. EXTRACT COORDINATES & CLEAN DATA
# -------------------------------
coords <- st_coordinates(cat_locs_sites_data$geometry)[, 1:2]

df <- cat_locs_sites_data %>%
  st_drop_geometry() %>%   # remove geometry for pie calculation
  mutate(
    x = coords[,1],
    y = coords[,2],
    RSV_A = as.numeric(RSV_A_prop),
    RSV_B = as.numeric(RSV_B_prop),
    total_population = as.numeric(total_population)
  ) %>%
  replace_na(list(RSV_A = 0, RSV_B = 0))  # replace NA with 0

# -------------------------------
# 3. COMPUTE PIE SLICES (long format + angles)
# -------------------------------
df_pies <- df %>%
  pivot_longer(cols = c("RSV_A", "RSV_B"), names_to = "variant", values_to = "fraction") %>%
  mutate(fraction = as.numeric(fraction)) %>%
  select(location, variant, fraction, x, y)

df_pies <- df_pies %>%
  group_by(location) %>%
  mutate(
    total = sum(fraction, na.rm = TRUE),
    prop  = fraction / total,
    start = c(0, head(cumsum(prop), -1)) * 2*pi,
    end   = cumsum(prop) * 2*pi
  ) %>%
  ungroup() %>%
  select(location, variant, fraction, x, y, start, end)

# -------------------------------

# ---- 4. CREATE BASE PIE RADIUS AND SCALE FOR FIT ----
max_pop <- max(df$total_population, na.rm = TRUE)

# initial radius based on population, 0.25 power transform to reduce extremes
df <- df %>%
  mutate(
    pop_r = (total_population / max_pop)^(1/4) * 70000  # initial max 70 km
  )

# then rescale slightly smaller for map fit
df <- df %>%
  mutate(
    pop_r_scaled = scales::rescale(pop_r, to = c(9000, 32000)) # 9–32 km
  )
df_pies <- df_pies %>%
  left_join(df %>% select(location, pop_r_scaled), by = "location")

# -------------------------------
# 5. COMPUTE MAP LIMITS WITH PADDING
# -------------------------------
x_range <- range(df$x)
y_range <- range(df$y)
x_pad <- diff(x_range) * 0.4 # padding
y_pad <- diff(y_range) * 0.5

# -------------------------------
# 6. BUILD FINAL MAP
# -------------------------------
pop_legend_values <- c(100000, 300000, 500000)  # example populations
pop_legend_radius <- (pop_legend_values / max_pop)^(1/4) * 70000
pop_legend_radius <- scales::rescale(pop_legend_radius, to = c(9000, 32000))
# Space them manually (top to bottom)
legend_x <- x_range[1] - x_pad * 0.6   # further right of map
# Desired vertical gap between edges
gap <- 8000  # adjust as needed

# Compute y positions from top to bottom
legend_y <- numeric(length(pop_legend_values))
legend_y[1] <- y_range[2] - pop_legend_radius[1]  # topmost circle

for(i in 2:length(pop_legend_values)){
  # previous bottom edge - gap - current radius
  legend_y[i] <- legend_y[i-1] - pop_legend_radius[i-1] - gap - pop_legend_radius[i]
}

# Create legend dataframe
legend_df <- data.frame(
  x = x_range[1] - x_pad*0.6,  # right of map
  y = legend_y + 50000,
  r = pop_legend_radius,
  label = paste0(pop_legend_values/1000, "k")
)

plot <- ggplot() +
  
  # Base map layers
  geom_sf(data = canton_boundaries, fill = NA, linewidth = 0.3, color= "grey") +
  geom_sf(data = catch_poly_sites, fill = NA, linewidth = 0.3, color= "gray42") +
  
  # Draw cake charts
  geom_arc_bar(
    data = df_pies,
    aes(
      x0 = x, y0 = y, r0 = 0, r = pop_r_scaled,
      start = start, end = end, fill = variant
    ),
    alpha = 0.6,
    color = "black",
    inherit.aes = FALSE
  ) +
  
  # Slice colors
  scale_fill_manual(values = c("RSV_A" = "blue", "RSV_B" = "green3"), name = "Subtype") +
  
  # Labels above pies
  geom_text(
    data = df,
    aes(x = x, y = y + pop_r_scaled + 7000, label = location),
    size = 5,
    inherit.aes = FALSE
  ) +
  
  geom_circle(
    data = legend_df,
    aes(x0 = x, y0 = y, r = r),
    color = "grey",
    fill = "grey90",
    linewidth = 0.5
  ) +
  geom_text(
    data = legend_df,
    aes(x = x, y = y, label = label),
    hjust = 0,
    vjust = 0.5,
    size = 4
  )+
  
  # Map limits with padding
  coord_sf(
    xlim = c(x_range[1] - x_pad, x_range[2] + x_pad),  # extra space for legend
    ylim = c(y_range[1] - y_pad, y_range[2] + y_pad),
    expand = FALSE
  ) +
  
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = c(0.05, 0.2),  # left (x=0.05), lower (y=0.2)
    
    legend.text = element_text(size = 14),  # increase size of legend labels
    legend.title = element_text(size = 14)
  )

ggsave('RSV_plot_cake.pdf')

print(plot)


