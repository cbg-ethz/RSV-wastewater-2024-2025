library(lme4)
library(ggplot2)
library(dplyr)
library(robustbase)
library(stringr)
library(tidyr)
library(ggpubr)



# load the data with genetic diversity estimates per sample
rsva_diversity_matrix <- read.csv("../../../RSV/data_analysis/results/RSVA_2024_2025/genomic_diversity/Switzerland-RSV-A-2024-2025_v1.csv")
rsva_diversity_matrix <- rsva_diversity_matrix[, c("date", 
                                                   "WWTP", 
                                                   "shannon_entropy", 
                                                   "covered_genome",
                                                   "covered_F",
                                                   "covered_G",
                                                   "genome_location"
                                                   )]

rsva_diversity_matrix$Date <- as.Date(rsva_diversity_matrix$date,format = "%Y_%m_%d")
rsva_diversity_matrix <- rsva_diversity_matrix %>%
  select(-date) %>%
  mutate (subtype = "A")

rsvb_diversity_matrix <- read.csv("../../../RSV/data_analysis/results/RSVB_2024_2025/genomic_diversity/Switzerland-RSV-B-2024-2025_v1.csv")
rsvb_diversity_matrix <- rsvb_diversity_matrix[, c("date", 
                                                   "WWTP", 
                                                   "shannon_entropy", 
                                                   "covered_genome",
                                                   "covered_F",
                                                   "covered_G",
                                                   "genome_location"
)]

rsvb_diversity_matrix$Date <- as.Date(rsvb_diversity_matrix$date,format = "%Y_%m_%d")
rsvb_diversity_matrix <- rsvb_diversity_matrix %>%
  select(-date) %>% 
  mutate(subtype = "B")

# merge genetic diversity dataframes
rsv_div_matrix <- rbind(rsva_diversity_matrix, rsvb_diversity_matrix)
# center the date
rsv_div_matrix$Date_cent = scale(as.numeric(as.Date(rsv_div_matrix$Date)), center = TRUE, scale = FALSE)[,1]
rsv_div_matrix <- rsv_div_matrix %>%
  mutate(WWTP = factor(WWTP)) %>%
  mutate(genome_location = factor(genome_location))
rsv_div_matrix$genome_location <- relevel(
  factor(rsv_div_matrix$genome_location), 
  ref = "full"
)


rsv_div_matrix <- rsv_div_matrix %>% filter(
  (covered_genome >= 0.5 & genome_location == 'full')
)




model <- lm(shannon_entropy ~ genome_location + subtype + Date_cent + WWTP,
              data = rsv_div_matrix)

summary(model)



# Get predictions from the model
#rsv_diversity_matrix$predictions <- predict(model)





# Shannon Entropy vs Time
ggplot(rsv_div_matrix[rsv_div_matrix$genome_location=="full",], aes(x = Date, y = shannon_entropy, color=subtype)) +
  geom_point(alpha=0.7, size=2) +                                # Add points
  geom_smooth(method = "lm", se = FALSE, linewidth=1) +  # Add a linear regression line
  #facet_wrap(~ WWTP, scales = "free_y") + # Facet by genome_location with separate y-axes
  theme_minimal(base_size = 14) +                             # Use a clean theme
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y")+  # Set x-axis ticks every month
  scale_color_manual(values = c("A" = "blue", "B" = "green")) +

  labs(
    title = expression("Shannon Entropy over time. Coverage " >= 50*"%"),
    x = "Time",
    y = "Mean Positional Shannon Entropy",
    color = "Subtype"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

#plot(g)





p <- ggplot(rsv_div_matrix, aes(x = genome_location, y = shannon_entropy, fill = subtype)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, linewidth = 0.7) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("A" = "blue", "B" = "green")) +
  scale_fill_manual(values = c("A" = "blue", "B" = "green")) +
  theme_minimal() +
  labs(
    x = "Genome region",
    y = "Mean Positional Shannon Entropy",
    title = expression("Shannon Entropy by Genome Location and Subtype"),
    color = "Subtype"
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 1.5),
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10)  # add top margin to prevent clipping
  )+
  facet_wrap(~ WWTP, scales = "free_y")  # Facet by genome_location with separate y-axes
  

plot(p)

# Shannon Entropy by RSV Subtype
# k <- ggplot(rsv_div_matrix[rsv_div_matrix$genome_location=="full",], aes(x = subtype, y = shannon_entropy, fill = subtype)) +
#   geom_boxplot(alpha = 0.8) +
#   geom_jitter(width = 0.2, alpha = 0.5) +
#   scale_fill_manual(values = c("A" = "blue", "B" = "green")) +
#   theme_minimal(base_size = 14) +
#   labs(title = expression("Mean Positional Shannon Entropy by RSV Subtype"), y = "Mean Shannon Entropy", x="Subtype") +
#   theme(
#     legend.position = "right",
#     legend.title = element_text(size = 14, face = "bold"),
#     legend.text = element_text(size = 14),
#     axis.title = element_text(size = 16, face = "bold"),
#     axis.text = element_text(size = 14),
#     plot.title = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 1.5),
#     plot.margin = margin(t = 20, r = 10, b = 10, l = 10)  # add top margin to prevent clipping
#   )
#plot(k)
###############################################################################

# Median Shannon Entropy: Calculate the median of the mean Shannon entropy values
# across all samples in that season. The median is less sensitive to outliers,
# which may be relevant if sample sizes or entropy distributions vary.

# Compute measures of variability (e.g., standard deviation, interquartile range) 
# to assess the spread of entropy values within each season.
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
# Quality control step - filter samples below 50% coverage


rsv_div_matrix <- rsv_div_matrix %>% filter(
  (covered_genome >= 0.5 & genome_location == 'full') |
  (covered_F >= 0.5 & genome_location == "F") |
  (covered_G >= 0.5 & genome_location == "G")
  
)


p <- ggplot(rsv_div_matrix, aes(x = genome_location, y = shannon_entropy, fill = subtype)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, linewidth = 0.7) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("A" = "blue", "B" = "green")) +
  scale_fill_manual(values = c("A" = "blue", "B" = "green")) +
  theme_minimal() +
  labs(
    x = "Genome region",
    y = "Mean Positional Shannon Entropy",
    title = expression("Shannon Entropy by Genome Location and Subtype. Coverage" >= 50*"%"),
    color = "Subtype"
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 1.5),
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10)  # add top margin to prevent clipping
  )+
facet_wrap(~ WWTP, scales = "free_y")  # Facet by genome_location with separate y-axes

plot(p)

# Shannon Entropy by RSV Subtype
k <- ggplot(rsv_div_matrix[rsv_div_matrix$genome_location=="full",], aes(x = subtype, y = shannon_entropy, fill = subtype)) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_fill_manual(values = c("A" = "blue", "B" = "green")) +
  theme_minimal(base_size = 14) +
  labs(title = expression("Mean Shannon Entropy by RSV Subtype. Coverage" >= 50*"%"), 
       y = "Mean Positional Shannon Entropy",
       x="Subtype") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 1.5),
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10)  # add top margin to prevent clipping
  )
#plot(k)

head(rsv_div_matrix)






# -------------------------------------------------
# RSV Shannon Entropy – FINAL ANALYSIS
# -------------------------------------------------
library(dplyr); library(emmeans); library(ggplot2)

df <- rsv_div_matrix %>%
  mutate(
    subtype         = factor(subtype),
    genome_location = factor(genome_location, levels = c("F","G","full")),
    WWTP            = factor(WWTP),
    Date_c          = as.numeric(Date - mean(Date)) / 7
  )

# OLS model (no random effects needed)
ols_mod <- lm(shannon_entropy ~ subtype * genome_location + Date_c, data = df)

# ANOVA table
anova(ols_mod) 

# Pairwise comparisons
cat("\n--- Gene differences within subtype ---\n")
print(pairs(emmeans(ols_mod, ~ genome_location | subtype), adjust="tukey"))

cat("\n--- Subtype differences within gene ---\n")
print(pairs(emmeans(ols_mod, ~ subtype | genome_location), adjust="tukey"))

# Plot
emm_df <- as.data.frame(emmeans(ols_mod, ~ genome_location | subtype))
ggplot(emm_df, aes(x=genome_location, y=emmean, colour=subtype)) +
  geom_point(position=position_dodge(0.3), size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL),
                width=.2, position=position_dodge(0.3)) +
  labs(y="Shannon Entropy (mean ± 95% CI)",
       title="RSV Diversity by Gene and Subtype") +
  theme_minimal()