library(ggplot2)
library(dplyr)
library(ggplot2)
library(robustbase)
library(stringr)
library(tidyr)


# RSV-A analysis
coverage_a_20250307_2418653583 <- read.csv(file = '/Users/arimaite/Documents/GitHub/RSV_wastewater/regular_monitoring_2024_2025/data/all_data/20250307_2418539124_and_20250307_2418653583_experiment/collected_rsv_coverage_20250307_2418653583_EPI_ISL_412866.tsv')
coverage_a_20250307_2418539124 <- read.csv(file = '/Users/arimaite/Documents/GitHub/RSV_wastewater/regular_monitoring_2024_2025/data/all_data/20250307_2418539124_and_20250307_2418653583_experiment/collected_rsv_coverage_20250307_2418539124_EPI_ISL_412866.tsv')
coverage_a_20250307_2418653583$batch <- "1:10:10"
coverage_a_20250307_2418539124$batch <- "1:1:1"
coverage_a <- rbind(coverage_a_20250307_2418653583, coverage_a_20250307_2418539124)

coverage_a$Date = str_extract(coverage_a$sample, '\\d{4}_\\d{2}_\\d{2}')
coverage_a$Date <- gsub("_", "-", coverage_a$Date)
coverage_a$location <- str_extract(coverage_a$sample,"\\d{2}")
print(head(coverage_a))
coverage_a$location <- gsub("10", "Zurich", coverage_a$location)
coverage_a$location <- gsub("16", "Geneva", coverage_a$location)
coverage_a$location <- gsub("15", "Basel", coverage_a$location)
coverage_a$location <- gsub("05", "Lugano", coverage_a$location)
coverage_a$location <- gsub("25", "Laupen", coverage_a$location)
coverage_a$location <- gsub("17", "Chur", coverage_a$location)
coverage_a <- coverage_a[c("pos", "coverage", "Date", "location", "batch")]


length <- max(coverage_a$pos)
coverage_a$Date <- as.Date(coverage_a$Date)  # Convert to Date format

coverage_a <- coverage_a %>%
  group_by(location, Date, batch) %>%
  mutate(covered = sum(coverage >= 10) / length) %>%
  mutate(mean_pos_cov = mean(coverage)) %>%
  distinct(location, Date, batch, .keep_all = TRUE)

# Create the plot
plot1 <- ggplot(coverage_a, aes(x = Date, y = covered, color = batch)) +
  geom_point() +
  geom_line()+
  labs(
    title = "RSV-A full genome",
    x = "Date",
    y = "Covered >=10x"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels by 45 degrees
  )+
  facet_wrap(~ location)  # Create facets by location

#print(plot1)



# Create the plot
plot2 <- ggplot(coverage_a, aes(x = Date, y = mean_pos_cov, color = batch)) +
  geom_point() +
  geom_line()+
  labs(
    title = "RSV-A full genome",
    x = "Date",
    y = "Mean position-wise coverage"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels by 45 degrees
  )+
  facet_wrap(~ location)  # Create facets by location

#print(plot2)
############################################### BAR PLOTS
library(ggplot2)
library(dplyr)
library(tidyr)

# Ensure Date is a factor to control bar spacing
coverage_a$Date <- as.factor(coverage_a$Date)

# Color palette (customizable)
my_colors <- c("1:1:1" = "#1f77b4", "1:10:10" = "#ff7f0e")  # Blue and orange

plot_bar_covered <- ggplot(coverage_a, aes(x = Date, y = covered, fill = batch)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  scale_fill_manual(values = my_colors, name = "Mix Ratio") +
  labs(
    title = "RSV-A Genome Coverage",
    subtitle = "Fraction of genome covered ≥10x",
    x = "Date",
    y = "Genome Fraction Covered (≥10x)"
  ) +
  facet_wrap(~ location, ncol = 3) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 13)
  )

print(plot_bar_covered)

# Ensure Date is a factor to control bar spacing
coverage_a$Date <- as.factor(coverage_a$Date)

# Color palette (customizable)
my_colors <- c("1:1:1" = "#1f77b4", "1:10:10" = "#ff7f0e")  # Blue and orange

plot_bar_positionwise_coverage <- ggplot(coverage_a, aes(x = Date, y = mean_pos_cov, fill = batch)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  scale_fill_manual(values = my_colors, name = "Mix Ratio") +
  labs(
    title = "RSV-A Genome Coverage",
    subtitle = "Mean position-wise coverage",
    x = "Date",
    y = "Mean position-wise coverage"
  ) +
  facet_wrap(~ location, ncol = 3) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 13)
  )

print(plot_bar_positionwise_coverage)


# Clean the environeent
rm(list = ls())


# RSV-B analysis
coverage_b_20250307_2418653583 <- read.csv(file = '/Users/arimaite/Documents/GitHub/RSV_wastewater/regular_monitoring_2024_2025/data/all_data/20250307_2418539124_and_20250307_2418653583_experiment/collected_rsv_coverage_20250307_2418653583_EPI_ISL_1653999.tsv')
coverage_b_20250307_2418539124 <- read.csv(file = '/Users/arimaite/Documents/GitHub/RSV_wastewater/regular_monitoring_2024_2025/data/all_data/20250307_2418539124_and_20250307_2418653583_experiment/collected_rsv_coverage_20250307_2418539124_EPI_ISL_1653999.tsv')
coverage_b_20250307_2418653583$batch <- "1:10:10"
coverage_b_20250307_2418539124$batch <- "1:1:1"
coverage_b <- rbind(coverage_b_20250307_2418653583, coverage_b_20250307_2418539124)

coverage_b$Date = coverage_b$sample#str_extract(coverage_b$sample, '\\d{4}_\\d{2}_\\d{2}')
#coverage_b$Date <- gsub("_", "-", coverage_b$Date)
coverage_b$location <- str_extract(coverage_b$sample,"\\d{2}")
print(head(coverage_b))
coverage_b$location <- gsub("10", "Zurich", coverage_b$location)
coverage_b$location <- gsub("16", "Geneva", coverage_b$location)
coverage_b$location <- gsub("15", "Basel", coverage_b$location)
coverage_b$location <- gsub("05", "Lugano", coverage_b$location)
coverage_b$location <- gsub("25", "Laupen", coverage_b$location)
coverage_b$location <- gsub("17", "Chur", coverage_b$location)
coverage_b <- coverage_b[c("pos", "coverage", "Date", "location", "batch")]


length <- max(coverage_b$pos)
#coverage_b$Date <- as.Date(coverage_b$Date)  # Convert to Date format
#coverage_b$Date <- format(coverage_b$Date, "%b %d")
coverage_b <- coverage_b %>%
  group_by(location, Date, batch) %>%
  mutate(covered = sum(coverage >= 10) / length) %>%
  mutate(mean_pos_cov = mean(coverage)) %>%
  distinct(location, Date, batch, .keep_all = TRUE)

# Create the plot
plot1 <- ggplot(coverage_b, aes(x = Date, y = covered, color = batch)) +
  geom_point() +
  geom_line()+
  labs(
    title = "RSV-B full genome",
    x = "Date",
    y = "Covered >=10x"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels by 45 degrees
  )+
  facet_wrap(~ location)  # Create facets by location

#print(plot1)



# Create the plot
plot2 <- ggplot(coverage_b, aes(x = Date, y = mean_pos_cov, color = batch)) +
  geom_point() +
  geom_line()+
  labs(
    title = "RSV-B full genome",
    x = "Date",
    y = "Mean position-wise coverage"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels by 45 degrees
  )+
  facet_wrap(~ location)  # Create facets by location

#print(plot2)
############################################### BAR PLOTS

# Ensure Date is a factor to control bar spacing
#coverage_b$Date <- as.factor(coverage_b$Date)

# Color palette (customizable)
my_colors <- c("1:1:1" = "#1f77b4", "1:10:10" = "#ff7f0e")  # Blue and orange

plot_bar_covered <- ggplot(coverage_b, aes(x = Date, y = covered, fill = batch)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  scale_fill_manual(values = my_colors, name = "Mix Ratio") +
  labs(
    title = "RSV-B Genome Coverage",
    subtitle = "Fraction of genome covered ≥10x",
    x = "Date",
    y = "Genome Fraction Covered (≥10x)"
  ) +
  facet_wrap(~ location, ncol = 3) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 13)
  )

print(plot_bar_covered)



# Color palette (customizable)
my_colors <- c("1:1:1" = "#1f77b4", "1:10:10" = "#ff7f0e")  # Blue and orange

plot_bar_positionwise_coverage <- ggplot(coverage_b, aes(x = Date, y = mean_pos_cov, fill = batch)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  scale_fill_manual(values = my_colors, name = "Mix Ratio") +
  labs(
    title = "RSV-B Genome Coverage",
    subtitle = "Mean position-wise coverage",
    x = "Date",
    y = "Mean position-wise coverage"
  ) +
  facet_wrap(~ location, ncol = 3) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 13)
  )

print(plot_bar_positionwise_coverage)



