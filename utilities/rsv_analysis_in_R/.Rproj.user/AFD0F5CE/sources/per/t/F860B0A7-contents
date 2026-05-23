library(ggplot2)
library(dplyr)
library(ggplot2)
library(robustbase)
library(stringr)
library(tidyr)
library(slider)


# ------------------------ List of the dates at which the samples were sequenced for each treatment plant ------------------------
samples <- read.csv("../../RSV_results/data_analysis/data/samples_tsv/samples_corrected_v2.tsv", 
                    sep = '\t',
                    header = TRUE)
samples$date <- as.character(samples$date)

samples$location <- gsub("10", "Zurich", samples$location)
samples$location <- gsub("16", "Geneva", samples$location)
samples$location <- gsub("15", "Basel", samples$location)
samples$location <- gsub("^5", "05", samples$location)
samples$location <- gsub("05", "Lugano", samples$location)
samples$location <- gsub("25", "Laupen", samples$location)
samples$location <- gsub("17", "Chur", samples$location)

samples <- samples %>%
  mutate(SampleName=str_extract(ID, "\\d{2}_\\d{4}_\\d{2}_\\d{2}"))
samples$date <- gsub("_", "-", samples$date)

# Have a list of samples that were sequenced for RSV-A:
rsva_samples <- samples %>%
  filter(!grepl("B", samples$protocol)) %>%
mutate(SampleName=str_extract(ID, "\\d{2}_\\d{4}_\\d{2}_\\d{2}"))

# Have a list of samples that were sequenced for RSV-B:
rsvb_samples <- samples %>%
  filter(!grepl("A", samples$protocol)) %>%
  mutate(SampleName=str_extract(ID, "\\d{2}_\\d{4}_\\d{2}_\\d{2}"))


# ------------------------ Concentration RSV-AB assay ------------------------


# Get the loads from file
rsv_ab_assay <- read.csv(file = '../../RSV_results/data_analysis/data/concentrations/rsv_ab_assay/251023_RSVAB.csv',
                                   check.names = FALSE) %>%
  select(-c(1)) %>%
  rename(Date = collection_date) %>%
  rename(location = ara_name) %>%
  select(Date, pathogen, location, load)


normalized_load_subtypes <- rsv_ab_assay %>%

  pivot_wider(
    names_from = pathogen,
    values_from = load
  )

normalized_load_subtypes$location <- gsub("Werdhölzli", "Zurich", normalized_load_subtypes$location)
normalized_load_subtypes$location <- gsub("Aire", "Geneva", normalized_load_subtypes$location)
normalized_load_subtypes$location <- gsub("ARA Basel", "Basel", normalized_load_subtypes$location)
normalized_load_subtypes$location <- gsub("CDA Lugano", "Lugano", normalized_load_subtypes$location)
normalized_load_subtypes$location <- gsub("ARA Sensetal", "Laupen", normalized_load_subtypes$location)
normalized_load_subtypes$location <- gsub("ARA Chur", "Chur", normalized_load_subtypes$location)

normalized_load_subtypes$Date <- as.character(as.Date(as.character(normalized_load_subtypes$Date),"%Y-%m-%d"))
normalized_load_subtypes <- normalized_load_subtypes[normalized_load_subtypes$location %in% c("Zurich", "Geneva", "Basel", "Lugano", "Laupen", "Chur"), ]


# ---------------- imputing missing RSV-AB normalized loads values if sample is sequenced but not sub-typed --------------
## these values are used for stacked plots of lineage relative abundances 
## scaled by subtype concentrations (imputed values are needed for the days with missing subtyping information)
# for each location

for (i in 1:nrow(samples)){
  date = samples[i, "date"]
  loc = samples[i, "location"]
  sample_id <- samples[i, "SampleName"]

  
  normalized_load_subtypes[(normalized_load_subtypes$Date == date & normalized_load_subtypes$location==loc), 'ID'] <- sample_id
}

# Filter samples that were not sequenced
normalized_load_subtypes <- normalized_load_subtypes %>%
  filter(!is.na(ID))



# add missing values for samples that were sequenced but not sub-typed
# for each location
for (i in 1:nrow(samples)){
  date = samples[i, "date"]
  loc = samples[i, "location"]
  sample_id <- samples[i, "SampleName"]
  
  if (!any(normalized_load_subtypes$Date==date & normalized_load_subtypes$location==loc)){
    
    # If not present, create a new row
    new_row <- data.frame(
      Date = date,
      location = loc,
      'RSV-A' = NA,
      'RSV-B' = NA,
      ID = sample_id,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    # Append new row to rsv_ab_assay
    normalized_load_subtypes <- rbind(normalized_load_subtypes, new_row)
  }
}


# Smooth the values of RSV/AB assay time series estimates and impute missing data

# 1. group by location and type
# RSV-A 
normalized_load_subtypes$Date <- as.Date(normalized_load_subtypes$Date, "%Y-%m-%d" )
normalized_load_subtypes_imputed <- normalized_load_subtypes %>%
  arrange(location, Date)  %>%
  group_by(location) %>%
  mutate(RSV_A_imputed = slide_index_dbl(
    .x = `RSV-A`,
    .i = Date,
    .f = ~median(.x, na.rm=TRUE),
    .before = lubridate::days(3),
    .after = lubridate::days(3),
    .complete = FALSE
  )) %>% ungroup()

# RSV-B 

normalized_load_subtypes_imputed <- normalized_load_subtypes_imputed %>%
  arrange(location, Date)  %>%
  group_by(location) %>%
  mutate(RSV_B_imputed = slide_index_dbl(
    .x = `RSV-B`,
    .i = Date,
    .f = ~median(.x, na.rm=TRUE),
    .before = lubridate::days(3),
    .after = lubridate::days(3),
    .complete = FALSE
  )) %>% ungroup()


# estimate proportions
normalized_load_subtypes_imputed <- normalized_load_subtypes_imputed %>%
  rename("RSV_A" = "RSV-A") %>%
  rename("RSV_B" = "RSV-B") %>% 
  mutate(RSVA_proportion = RSV_A_imputed/(RSV_A_imputed + RSV_B_imputed + 10e-10)) %>% # add a small constant to avoid division by zero
  mutate(RSVB_proportion = RSV_B_imputed/(RSV_A_imputed + RSV_B_imputed + 10e-10))
  


ggplot(data = normalized_load_subtypes_imputed, aes(as.Date(Date), y=RSV_A_imputed)) +
  geom_line()+
  geom_point()+
  facet_wrap(~location, scales = 'free_y')

ggplot(data = normalized_load_subtypes_imputed, aes(as.Date(Date), y=RSV_B_imputed)) +
  geom_line()+
  geom_point()+
  facet_wrap(~location, scales = 'free_y')



write.csv(normalized_load_subtypes_imputed, '../../RSV_results/data_analysis/data/concentrations/rsv_ab_assay/viral_loads_imputated_rsvab.csv')



# plot imputed
long_df <- normalized_load_subtypes_imputed %>%
  pivot_longer(
    cols = c(RSV_A_imputed, RSV_B_imputed),
    names_to = "type",
    values_to = "Concentration"
  )

ggplot(long_df, aes(x = as.Date(Date), y = Concentration, color=type)) +
  geom_line() +
  geom_point() +
  facet_wrap(~location, scales = "free_y") +
  labs(x = "Date", y = "RSV Concentration", title = "RSVA and RSVB Over Time") +
  theme_minimal()


# # # plot not imputed
# 
# long_df <- merged_df_normalized_pan %>%
#   pivot_longer(
#     cols = c(RSVA, RSVB),
#     names_to = "type",
#     values_to = "Concentration"
#   )
# 
# ggplot(long_df, aes(x = as.Date(Date), y = Concentration, color=type)) +
#   geom_line() +
#   geom_point() +
#   facet_wrap(~location, scales = "free_y") +
#   labs(x = "Date", y = "RSV Concentration", title = "RSVA and RSVB Over Time") +
#   theme_minimal()
