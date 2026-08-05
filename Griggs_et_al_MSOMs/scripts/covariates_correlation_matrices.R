##############################################
# Check MSOM covariate correlations
##############################################

# Purpose:
# 1. Check Pearson correlations among site-level occupancy covariates.
# 2. Check Pearson correlations among numeric detection covariates.
# 3. Check model-level collinearity among detection covariates using GVIF.
#
# Author: Johanna Griggs
# Last updated: 2026-07-26

# Clear environment (if needed)
rm(list = ls())

# Set working directory to location of 'Griggs_et_al_MSOMs' folder
setwd("C:/Users/johan/OneDrive - University of Victoria/MSc/HT Cam Project/Data/Griggs_et_al_MSOMs")
getwd()

# Load packages 
library(tidyverse)
library(corrplot)
library(car)

### LOAD DATASETS ### -------------------------------------------

# Run 3.0_prepare_MSOM_covariates.R to create the required model covariate files.

# For site covariates
occu_site_covs <- read.csv("inputs/occu_site_covs_biweekly_MSOM.csv")
site_year_list <- read.csv("inputs/site_year_list_biweekly.csv")

# For detection covariates
det_cov_trail_strength <- read.csv("inputs/det_cov_trail_strength.csv")
summer_effort_wide <- read.csv("inputs/summer_biweekly_cam_effort.csv")
winter_effort_wide <- read.csv("inputs/winter_biweekly_cam_effort.csv")
deer_summer_prop_det_df <- read.csv("inputs/det_cov_deer_summer_prop_det_biweekly.csv")
deer_winter_prop_det_df <- read.csv("inputs/det_cov_deer_winter_prop_det_biweekly.csv")
wolf_summer_prop_det_df <- read.csv("inputs/det_cov_wolf_summer_prop_det_biweekly.csv")
wolf_winter_prop_det_df <- read.csv("inputs/det_cov_wolf_winter_prop_det_biweekly.csv")
deer_binary_summer <- read.csv("inputs/det_cov_deer_summer_binary.csv")
deer_binary_winter <- read.csv("inputs/det_cov_deer_winter_binary.csv")
wolf_binary_summer <- read.csv("inputs/det_cov_wolf_summer_binary.csv")
wolf_binary_winter <- read.csv("inputs/det_cov_wolf_winter_binary.csv")

# Define effort column names
summer_eff_cols <- setdiff(names(summer_effort_wide),"SiteYear")
winter_eff_cols <- setdiff(names(winter_effort_wide),"SiteYear")

# QA/QC
stopifnot(identical(summer_eff_cols, paste0("E", 1:13)),
  identical(winter_eff_cols, paste0("E", c(1:6, 8:13))))

###################################
# CHECK CORRELATIONS BETWEEN SITE COVARIATES  
###################################

# Select covariates to be used in MSOMs
final_model_covs <- c(
  "recent_harvest_prop_buff_z",
  "second_growth_prop_buff_z",
  "old_growth_prop_buff_z",
  "coast_decay_z",
  "roads_decay_z",
  "elevation_z",
  "mean_tri_500_z")

# Confirm all expected covariates exist and SiteIDs are unique
stopifnot(
  all(final_model_covs %in% names(occu_site_covs)),
  "SiteYear" %in% names(occu_site_covs),
  !anyDuplicated(site_year_list$SiteYear))

candidate_occu_covs <- occu_site_covs %>%
  left_join(
    site_year_list %>% select(SiteYear, SiteID),
    by = "SiteYear") %>%
  distinct(SiteID, .keep_all = TRUE) %>%
  select(all_of(final_model_covs))

### CALCULATE PEARSON CORRELATIONS ### -----------------------------------------

candidate_occu_cor <- cor(
  candidate_occu_covs,
  method = "pearson",
  use = "complete.obs")

print(round(candidate_occu_cor, 2))

### PLOT CORRELATION MATRIX ### ------------------------------------------------

corrplot(
  candidate_occu_cor,
  method = "color",
  type = "upper",
  tl.cex = 0.9,
  tl.col = "black",
  addCoef.col = "black",
  number.cex = 1.5,
  diag = FALSE)

###################################
# CHECK CORRELATIONS BETWEEN DETECTION COVARIATES 
###################################

# BUILD DETECTION-COVARIATE DATAFRAMES

# Summer: binary covariates
summer_df_binary <- summer_effort_wide %>%
  pivot_longer(all_of(summer_eff_cols), names_to = "Occasion",
               values_to = "effort_days") %>%
  filter(effort_days > 0) %>%
  left_join(deer_binary_summer, by = "SiteYear") %>%
  left_join(wolf_binary_summer, by = "SiteYear") %>%
  left_join(det_cov_trail_strength, by = "SiteYear") %>%
  rename(deer_binary = deer_binary_summer,
         wolf_binary = wolf_binary_summer,
         trail = trail_strength) %>%
  mutate(trail = factor(trail))

# Summer: proportion-detected covariates
summer_df_prop <- summer_effort_wide %>%
  pivot_longer(all_of(summer_eff_cols), names_to = "Occasion",
               values_to = "effort_days") %>%
  filter(effort_days > 0) %>%
  left_join(deer_summer_prop_det_df, by = "SiteYear") %>%
  left_join(wolf_summer_prop_det_df, by = "SiteYear") %>%
  left_join(det_cov_trail_strength, by = "SiteYear") %>%
  rename(deer_prop_det = deer_summer_prop_det,
         wolf_prop_det = wolf_summer_prop_det,
         trail = trail_strength) %>%
  mutate(trail = factor(trail))

# Winter: binary covariates
winter_df_binary <- winter_effort_wide %>%
  pivot_longer(all_of(winter_eff_cols), names_to = "Occasion",
               values_to = "effort_days") %>%
  filter(effort_days > 0) %>%
  left_join(deer_binary_winter, by = "SiteYear") %>%
  left_join(wolf_binary_winter, by = "SiteYear") %>%
  left_join(det_cov_trail_strength, by = "SiteYear") %>%
  rename(deer_binary = deer_binary_winter,
         wolf_binary = wolf_binary_winter,
         trail = trail_strength) %>%
  mutate(trail = factor(trail))

# Winter: proportion-detected covariates
winter_df_prop <- winter_effort_wide %>%
  pivot_longer(all_of(winter_eff_cols), names_to = "Occasion",
               values_to = "effort_days") %>%
  filter(effort_days > 0) %>%
  left_join(deer_winter_prop_det_df, by = "SiteYear") %>%
  left_join(wolf_winter_prop_det_df, by = "SiteYear") %>%
  left_join(det_cov_trail_strength, by = "SiteYear") %>%
  rename(deer_prop_det = deer_winter_prop_det,
         wolf_prop_det = wolf_winter_prop_det,
         trail = trail_strength) %>%
  mutate(trail = factor(trail))

# Create one corrplot with both binary and proportion detection covariates for simplicity
    # Note: Binary detection and proportion detected are alternative formulations of
    # the same cross-species detection covariate and will not be included together
    # in the same candidate detection model.

# We are using:
  # Pearson correlations for numeric covariates
  # GVIF for all covariates including trail strength (categorical)

# Summer
summer_cor_combined_dat <- summer_df_binary %>%
  select(SiteYear, Occasion, effort_days, deer_binary, wolf_binary) %>%
  left_join(summer_df_prop %>% select(SiteYear, Occasion, deer_prop_det, wolf_prop_det),
            by = c("SiteYear", "Occasion")) %>%
  select(effort_days, deer_binary, deer_prop_det, wolf_binary, wolf_prop_det) %>%
  drop_na()

summer_cor_combined <- cor(summer_cor_combined_dat, method = "pearson")
print(round(summer_cor_combined, 2))

corrplot(summer_cor_combined, method = "color", type = "upper",
         addCoef.col = "black", number.cex = 1.5, tl.col = "black", tl.srt = 45,
         diag = FALSE, mar = c(0, 0, 1, 0),
         title = "Summer numeric detection covariates")

# Winter
winter_cor_combined_dat <- winter_df_binary %>%
  select(SiteYear, Occasion, effort_days, deer_binary, wolf_binary) %>%
  left_join(winter_df_prop %>% select(SiteYear, Occasion, deer_prop_det, wolf_prop_det),
            by = c("SiteYear", "Occasion")) %>%
  select(effort_days, deer_binary, deer_prop_det, wolf_binary, wolf_prop_det) %>%
  drop_na()

winter_cor_combined <- cor(winter_cor_combined_dat, method = "pearson")
print(round(winter_cor_combined, 2))

corrplot(winter_cor_combined, method = "color", type = "upper",
         addCoef.col = "black", number.cex = 1.5, tl.col = "black", tl.srt = 45,
         diag = FALSE, mar = c(0, 0, 1, 0),
         title = "Winter numeric detection covariates")

# GVIF model-level collinearity diagnostics

# GVIF checks model-level collinearity among covariates that could actually be included together in the same species-specific detection submodel.
# Deer detection models include wolf detection covariates.
# Wolf detection models include deer detection covariates.

# Summer: binary detection covariates

# Deer detection: effort_days + wolf_binary + trail
summer_vif_deer_binary_dat <- summer_df_binary %>%
  select(effort_days, wolf_binary, trail) %>%
  drop_na() %>%
  mutate(trail = droplevels(trail), dummy = seq_len(n()))

summer_vif_deer_binary <- car::vif(
  lm(dummy ~ effort_days + wolf_binary + trail, data = summer_vif_deer_binary_dat))

summer_vif_deer_binary_table <- as.data.frame(summer_vif_deer_binary) %>%
  rownames_to_column("term") %>%
  rename(GVIF_adj = `GVIF^(1/(2*Df))`) %>%
  mutate(Dataset = "Summer binary", Detection_model = "Deer detection") %>%
  select(Dataset, Detection_model, term, GVIF, Df, GVIF_adj)

summer_vif_deer_binary_table

# Wolf detection: effort_days + deer_binary + trail
summer_vif_wolf_binary_dat <- summer_df_binary %>%
  select(effort_days, deer_binary, trail) %>%
  drop_na() %>%
  mutate(trail = droplevels(trail), dummy = seq_len(n()))

summer_vif_wolf_binary <- car::vif(
  lm(dummy ~ effort_days + deer_binary + trail, data = summer_vif_wolf_binary_dat))

summer_vif_wolf_binary_table <- as.data.frame(summer_vif_wolf_binary) %>%
  rownames_to_column("term") %>%
  rename(GVIF_adj = `GVIF^(1/(2*Df))`) %>%
  mutate(Dataset = "Summer binary", Detection_model = "Wolf detection") %>%
  select(Dataset, Detection_model, term, GVIF, Df, GVIF_adj)

summer_vif_wolf_binary_table

# Summer: proportion-detected covariates

# Deer detection: effort_days + wolf_prop_det + trail
summer_vif_deer_prop_dat <- summer_df_prop %>%
  select(effort_days, wolf_prop_det, trail) %>%
  drop_na() %>%
  mutate(trail = droplevels(trail), dummy = seq_len(n()))

summer_vif_deer_prop <- car::vif(
  lm(dummy ~ effort_days + wolf_prop_det + trail, data = summer_vif_deer_prop_dat))

summer_vif_deer_prop_table <- as.data.frame(summer_vif_deer_prop) %>%
  rownames_to_column("term") %>%
  rename(GVIF_adj = `GVIF^(1/(2*Df))`) %>%
  mutate(Dataset = "Summer proportion", Detection_model = "Deer detection") %>%
  select(Dataset, Detection_model, term, GVIF, Df, GVIF_adj)

summer_vif_deer_prop_table

# Wolf detection: effort_days + deer_prop_det + trail
summer_vif_wolf_prop_dat <- summer_df_prop %>%
  select(effort_days, deer_prop_det, trail) %>%
  drop_na() %>%
  mutate(trail = droplevels(trail), dummy = seq_len(n()))

summer_vif_wolf_prop <- car::vif(
  lm(dummy ~ effort_days + deer_prop_det + trail, data = summer_vif_wolf_prop_dat))

summer_vif_wolf_prop_table <- as.data.frame(summer_vif_wolf_prop) %>%
  rownames_to_column("term") %>%
  rename(GVIF_adj = `GVIF^(1/(2*Df))`) %>%
  mutate(Dataset = "Summer proportion", Detection_model = "Wolf detection") %>%
  select(Dataset, Detection_model, term, GVIF, Df, GVIF_adj)

summer_vif_wolf_prop_table

# Winter: binary detection covariates

# Deer detection: effort_days + wolf_binary + trail
winter_vif_deer_binary_dat <- winter_df_binary %>%
  select(effort_days, wolf_binary, trail) %>%
  drop_na() %>%
  mutate(trail = droplevels(trail), dummy = seq_len(n()))

winter_vif_deer_binary <- car::vif(
  lm(dummy ~ effort_days + wolf_binary + trail, data = winter_vif_deer_binary_dat))

winter_vif_deer_binary_table <- as.data.frame(winter_vif_deer_binary) %>%
  rownames_to_column("term") %>%
  rename(GVIF_adj = `GVIF^(1/(2*Df))`) %>%
  mutate(Dataset = "Winter binary", Detection_model = "Deer detection") %>%
  select(Dataset, Detection_model, term, GVIF, Df, GVIF_adj)

winter_vif_deer_binary_table

# Wolf detection: effort_days + deer_binary + trail
winter_vif_wolf_binary_dat <- winter_df_binary %>%
  select(effort_days, deer_binary, trail) %>%
  drop_na() %>%
  mutate(trail = droplevels(trail), dummy = seq_len(n()))

winter_vif_wolf_binary <- car::vif(
  lm(dummy ~ effort_days + deer_binary + trail, data = winter_vif_wolf_binary_dat))

winter_vif_wolf_binary_table <- as.data.frame(winter_vif_wolf_binary) %>%
  rownames_to_column("term") %>%
  rename(GVIF_adj = `GVIF^(1/(2*Df))`) %>%
  mutate(Dataset = "Winter binary", Detection_model = "Wolf detection") %>%
  select(Dataset, Detection_model, term, GVIF, Df, GVIF_adj)

winter_vif_wolf_binary_table

# Winter: proportion-detected covariates

# Deer detection: effort_days + wolf_prop_det + trail
winter_vif_deer_prop_dat <- winter_df_prop %>%
  select(effort_days, wolf_prop_det, trail) %>%
  drop_na() %>%
  mutate(trail = droplevels(trail), dummy = seq_len(n()))

winter_vif_deer_prop <- car::vif(
  lm(dummy ~ effort_days + wolf_prop_det + trail, data = winter_vif_deer_prop_dat))

winter_vif_deer_prop_table <- as.data.frame(winter_vif_deer_prop) %>%
  rownames_to_column("term") %>%
  rename(GVIF_adj = `GVIF^(1/(2*Df))`) %>%
  mutate(Dataset = "Winter proportion", Detection_model = "Deer detection") %>%
  select(Dataset, Detection_model, term, GVIF, Df, GVIF_adj)

winter_vif_deer_prop_table

# Wolf detection: effort_days + deer_prop_det + trail
winter_vif_wolf_prop_dat <- winter_df_prop %>%
  select(effort_days, deer_prop_det, trail) %>%
  drop_na() %>%
  mutate(trail = droplevels(trail), dummy = seq_len(n()))

winter_vif_wolf_prop <- car::vif(
  lm(dummy ~ effort_days + deer_prop_det + trail, data = winter_vif_wolf_prop_dat))

winter_vif_wolf_prop_table <- as.data.frame(winter_vif_wolf_prop) %>%
  rownames_to_column("term") %>%
  rename(GVIF_adj = `GVIF^(1/(2*Df))`) %>%
  mutate(Dataset = "Winter proportion", Detection_model = "Wolf detection") %>%
  select(Dataset, Detection_model, term, GVIF, Df, GVIF_adj)

winter_vif_wolf_prop_table

# # Combine GVIF results
gvif_diagnostics <- bind_rows(
  summer_vif_deer_binary_table,
  summer_vif_wolf_binary_table,
  summer_vif_deer_prop_table,
  summer_vif_wolf_prop_table,
  winter_vif_deer_binary_table,
  winter_vif_wolf_binary_table,
  winter_vif_deer_prop_table,
  winter_vif_wolf_prop_table) %>%
  mutate(
    collinearity_interpretation = case_when(
      GVIF_adj < 2 ~ "low",
      GVIF_adj < 3 ~ "moderate",
      GVIF_adj >= 3 ~ "potentially high"))

gvif_diagnostics

