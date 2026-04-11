library(ggplot2)
library(dplyr)
library(ggplot2)
library(robustbase)
library(stringr)
library(tidyr)



# RSV-B analysis
# List of the dates at which the samples were sequenced for each treatment plant
samples <- read.csv("../../RSV_results/data_analysis/data/samples_tsv/samples_corrected_v2.tsv", 
                    sep = '\t',
                    header = TRUE,
                    colClasses = c(
                      "location" = "character"
                    ))


samples$location <- gsub("10", "Zurich", samples$location)
samples$location <- gsub("16", "Geneva", samples$location)
samples$location <- gsub("15", "Basel", samples$location)
samples$location <- gsub("05", "Lugano", samples$location)
samples$location <- gsub("25", "Laupen", samples$location)
samples$location <- gsub("17", "Chur", samples$location)
samples$date <- gsub("_", "-", samples$date)
samples$date <- as.character(samples$date)

# Have a list of samples that were sequenced for RSV-A (remove only B):
rsvb_samples <- samples %>% 
  filter(!grepl("A", samples$protocol))
rsvb_samples$date <- as.character(rsvb_samples$date)  # Convert to Date format



# Coverage
coverage_b <- read.csv('../../RSV_results/data_analysis/data/RSVB_2024_2025/coverage/collected_rsv_coverage_all_2024_2025_EPI_ISL_1653999.tsv')
coverage_b$Date = str_extract(coverage_b$sample, '\\d{4}_\\d{2}_\\d{2}')
coverage_b$Date <- gsub("_", "-", coverage_b$Date)



coverage_b$location <- str_extract(coverage_b$sample,"\\d{2}")
print(head(coverage_b))
coverage_b$location <- gsub("10", "Zurich", coverage_b$location)
coverage_b$location <- gsub("16", "Geneva", coverage_b$location)
coverage_b$location <- gsub("15", "Basel", coverage_b$location)
coverage_b$location <- gsub("05", "Lugano", coverage_b$location)
coverage_b$location <- gsub("25", "Laupen", coverage_b$location)
coverage_b$location <- gsub("17", "Chur", coverage_b$location)
coverage_b <- coverage_b[c("pos", "coverage", "Date", "location")]

length_rsvb_genome <- max(coverage_b$pos) #total number of genome bases
F_coordinates <- c(5676, 7400)
F_length = F_coordinates[2] - F_coordinates[1] + 1
coverage_b_plotting <- coverage_b %>%
  group_by(location, Date) %>%
  mutate(covered = sum(coverage >= 10) / length_rsvb_genome) %>%
  mutate(mean_pos_cov = mean(coverage)) %>%
  mutate(F_covered = sum(coverage >= 10 & pos >= F_coordinates[1] & pos <= F_coordinates[2]) / F_length) %>%
  distinct(location, Date, .keep_all = TRUE)  %>%
  ungroup() %>% select(-c("coverage", "pos"))


#print(coverage_b %>%
#        group_by(location, Date) %>%
#        summarize(count = n()), n = 200)


################################
# If sample was sequenced but because of 0 coverage it was not included into coverage file, include as zeros.
coverage_b_plotting <- coverage_b_plotting %>%
  full_join(
    rsvb_samples %>%  rename(Date = date),
    by = c("location", "Date")
  ) %>%
  replace_na(list(
    covered = 0,
    mean_pos_cov = 0,
    F_covered = 0
  ))
################################

# SUBTYPE-SPECIFIC Concentration


color_palette <- c(
  "Zurich" = "#1f78b4",
  "Geneva" = "#ff7f00",
  "Basel" = "#33a02c",
  "Lugano" = "#e31a1c",
  "Laupen" = "#6a3d9a",
  "Chur" = "#b15928"
)
# Get the concentrations from RSV-A/B assay



# Get the concentrations from RSV-A/B assay
rsv_ab_assay <- read.csv(file = '../../RSV_results/data_analysis/data/concentrations/rsv_ab_assay/251023_RSVAB.csv',
                         check.names = FALSE) %>%
  select(-c(1)) %>%
  rename(Date = collection_date) %>%
  rename(location = ara_name) %>%
  select(Date, pathogen, location, conc) %>%
  mutate(conc = conc / 1000) # convert gc/l -> gc/ml


concentration_subtypes <- rsv_ab_assay %>%
  
  pivot_wider(
    names_from = pathogen,
    values_from = conc
  ) %>%
  rename("RSV_A" = "RSV-A") %>%
  rename("RSV_B" = "RSV-B")

concentration_subtypes$location <- gsub("Werdhölzli", "Zurich", concentration_subtypes$location)
concentration_subtypes$location <- gsub("Aire", "Geneva", concentration_subtypes$location)
concentration_subtypes$location <- gsub("ARA Basel", "Basel", concentration_subtypes$location)
concentration_subtypes$location <- gsub("CDA Lugano", "Lugano", concentration_subtypes$location)
concentration_subtypes$location <- gsub("ARA Sensetal", "Laupen", concentration_subtypes$location)
concentration_subtypes$location <- gsub("ARA Chur", "Chur", concentration_subtypes$location)

concentration_subtypes$Date <- as.character(as.Date(as.character(concentration_subtypes$Date),"%Y-%m-%d"))
concentration_subtypes <- concentration_subtypes[concentration_subtypes$location %in% c("Zurich", "Geneva", "Basel", "Lugano", "Laupen", "Chur"), ]


concentration_subtypes <- concentration_subtypes[,c('location', 'Date', 'RSV_A', 'RSV_B')]

concentration_subtypes$Date <- as.character(concentration_subtypes$Date)


coverage_subtype_conc <- inner_join(coverage_b_plotting, concentration_subtypes, by = c("location", "Date"))[c('Date', 'location','covered','F_covered', 'mean_pos_cov', 'RSV_B')]
coverage_subtype_conc$location <- as.factor(coverage_subtype_conc$location)





coverage_subtype_conc$ln_covered_modified <- log(1/(1-coverage_subtype_conc$covered))

robust_exponential_models <- coverage_subtype_conc %>%
  group_by(location) %>%
  summarise(
    model = list(lmrob(ln_covered_modified ~ 0 + RSV_B)),
    slope = coef(model[[1]])[1],
    spearman_rho = cor.test(RSV_B, covered, method = 'spearman')$estimate,
    p_value = cor.test(RSV_B, covered, method = 'spearman')$p.value,
    max_RSVB = max(RSV_B, na.rm = TRUE)
  )



exp_curve_data <- robust_exponential_models %>%
  group_by(location) %>%
  tidyr::expand(cc_dat = seq(0, max_RSVB, max_RSVB/20)) %>%  # Expand grid for smooth curve
  left_join(robust_exponential_models, by = "location") %>%
  mutate(predicted_y = 1 - exp(-slope * cc_dat))





annotation_df <- robust_exponential_models %>%
  select(location, slope, spearman_rho, p_value) %>%
  mutate(
    x = max(coverage_subtype_conc$RSV_B),  # or dynamically max(coverage_subtype_conc$`RSV-A`)
    y1 = 0.1,
    y2 = 0.02,
    #label_slope = paste0("italic(k)==", round(slope, 3)),
    label_rho   = paste0("italic(rho)==", round(spearman_rho, 3)),
    label_p     = paste0("italic(p)==", signif(p_value, 3))
  )



plot_covered <- ggplot(coverage_subtype_conc, aes(x = RSV_B, y = covered, color = location)) +
  geom_point() +
  ylim(c(0,1)) +
  
  geom_point(aes(x = RSV_B, y = covered, color = location, shape = "Whole genome"), size = 2.5) +
  geom_point(aes(x = RSV_B, y = F_covered, shape = "F gene"), size = 2.5, color='black', alpha = 0.3) +
  geom_line(data = exp_curve_data, aes(x = cc_dat, y = predicted_y, color = location), linetype = "dashed") +
  
  
  labs(
    title = "", 
    x = "RSV-B, gc/mL wastewater", #"gc/mL wastewater", 
    y = bquote("Covered " >= 10 ~ "×")
  ) +
  theme_minimal()  + 
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5), # bigger title
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 16),
    strip.text = element_text(size = 18), # facet label
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  ) +
  facet_wrap(~ location)  + 
  
  geom_text(
    data = annotation_df,
    aes(x = x, y = y1, label = label_rho),
    hjust = 1, parse = TRUE, inherit.aes = FALSE, size =5
  ) +
  geom_text(
    data = annotation_df,
    aes(x = x, y = y2, label = label_p),
    hjust = 1, parse = TRUE, inherit.aes = FALSE, size =5
  ) +
  scale_color_manual(
    name = "Location",
    values = color_palette
  ) +
  scale_shape_manual(
    name = "",
    values = c("Whole genome" = 16, "F gene" = 17)
  )

print(plot_covered)

# Estimate correlations
# Estimate correlation between concentrations estimated with RSV-AB subtyping assay and the proportion of genome covered
correlation_all_locations <-  cor.test(coverage_subtype_conc$RSV_B, 
                                       coverage_subtype_conc$covered,
                                       method = 'spearman')
print(correlation_all_locations)
# coverage_b_plotting contains all RSV-B samples, whereas coverage_subtype_conc contains only the samples for which we have RSV-AB assay measurements
# Estimate median and range of proportion of genome covered
coverage_subtype_conc_nonzero_conc <- coverage_subtype_conc %>% 
  filter(
    RSV_B > 0
  )
coverage_subtype_conc_nonzero_conc <- coverage_subtype_conc_nonzero_conc %>%
  summarise(
    median_covered = median(covered, na.rm = TRUE),
    min_covered = min(covered),
    max_covered = max(covered),
  )
print(coverage_subtype_conc_nonzero_conc)

correlation_all_locations_F_GENE <-  cor.test(coverage_subtype_conc$RSV_B, 
                                              coverage_subtype_conc$F_covered,
                                              method = 'spearman')
print(correlation_all_locations_F_GENE)

# Estimate median coverage: RSV-AB ASSAY FILTERING
#***with non-zero RSV-B subtype concentration
#########################################
coverage_subtype_conc_filtered <- coverage_subtype_conc %>% 
  filter(
    RSV_B >= 6
  )
coverage_subtype_conc_filtered <- coverage_subtype_conc_filtered %>%

    summarise(
    median_covered = median(covered, na.rm = TRUE),
    min_covered = min(covered),
    max_covered = max(covered),
  )
print(coverage_subtype_conc_filtered)

# Compare locations - no filtering
print("Compare locations")
comp_loc_coverage <- coverage_subtype_conc %>%
  group_by(location) %>% 
  summarise(
    median_covered = median(covered, na.rm = TRUE),
    min_covered = min(covered),
    max_covered = max(covered),
  )
print(comp_loc_coverage)

