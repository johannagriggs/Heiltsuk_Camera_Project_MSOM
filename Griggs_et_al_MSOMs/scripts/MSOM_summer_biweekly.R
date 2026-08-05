##############################################
# Run Summer MSOM (with biweekly sampling occasions, dropping any partial occasions)
##############################################

# Purpose of script: Running a single-season (summer), multi-species occupancy model 
  # Code is adapted from code from Dyck et al. 2022 (https://github.com/marissadyck/Multi_occ_LWW/tree/main)
# Author: Johanna Griggs
# Date Created: 2026-05-05

# Data structure:
# We defined:
  # Summer as April-September; with dropped occasions the final dates are April 1st to September 29th

# Remove all existing data in R session (if needed)
rm(list = ls())

###################################
# 1.0 Load data
###################################

# # Optional: install "unmarked" Version 1.4.1 (newer versions may use different syntax)
# install.packages("remotes")
# library(remotes)
# remotes::install_version("unmarked",version="1.4.1")

# # Load packages
library(unmarked)
library(tidyverse)
library(magrittr)
library(ggplot2)
library(boot)
library(AICcmodavg)
library(patchwork)

# Set working directory to location of 'Griggs_et_al_MSOMs' folder
setwd("C:/Users/johan/OneDrive - University of Victoria/MSc/HT Cam Project/Data/Griggs_et_al_MSOMs")
getwd()

# Source code provided from Ken Kellner via Dyck et al. to fix package bugs related to the predict function
source("inputs/om_predict_fix.R")

# # Load data 

# Deer detection history 
deer_summer_det_hist <- read.csv("inputs/deer_summer_det_hist_biweekly.csv", row.names = "SiteYear")
head(deer_summer_det_hist) # ensure first column is not interpreted as data

# Wolf detection history 
wolf_summer_det_hist <- read.csv("inputs/wolf_summer_det_hist_biweekly.csv", row.names = "SiteYear")
head(wolf_summer_det_hist)

# Convert detection histories to matrices 
deer_sum_det_matrix <- as.matrix(deer_summer_det_hist)
wolf_sum_det_matrix <- as.matrix(wolf_summer_det_hist)
class(deer_sum_det_matrix)
class(wolf_sum_det_matrix)

# Load occupancy covariates
occu_site_covs <- read.csv("inputs/occu_site_covs_biweekly_MSOM.csv", row.names = "SiteYear")
colnames(occu_site_covs) 

# Load untransformed covariates (for backtransforming axes later)
orig_site_covs <- read.csv("inputs/orig_site_covs_biweekly_MSOM.csv", row.names = "SiteYear")

# Trail_strength (FAINT, MODERATE, STRONG levels)
trail_strength <- read.csv("inputs/det_cov_trail_strength.csv",  row.names = "SiteYear")

# Camera effort 
effort_biweekly <- read.csv("inputs/summer_biweekly_cam_effort.csv",  row.names = "SiteYear")

# Species-binary 
deer_binary_summer <-  read.csv("inputs/det_cov_deer_summer_binary.csv",  row.names = "SiteYear")
wolf_binary_summer <-  read.csv("inputs/det_cov_wolf_summer_binary.csv",  row.names = "SiteYear")

# Species proportion detected 
deer_summer_prop_det <- read.csv("inputs/det_cov_deer_summer_prop_det_biweekly.csv",row.names = "SiteYear")
wolf_summer_prop_det <- read.csv("inputs/det_cov_wolf_summer_prop_det_biweekly.csv",row.names = "SiteYear")

# Ensure SiteYear rows are in the same order
stopifnot(
  identical(rownames(deer_summer_det_hist), rownames(wolf_summer_det_hist)),
  identical(rownames(deer_summer_det_hist), rownames(occu_site_covs)),
  identical(rownames(deer_summer_det_hist), rownames(orig_site_covs)),
  identical(rownames(deer_summer_det_hist), rownames(trail_strength)),
  identical(rownames(deer_summer_det_hist), rownames(effort_biweekly)),
  identical(rownames(deer_summer_det_hist), rownames(deer_binary_summer)),
  identical(rownames(deer_summer_det_hist), rownames(wolf_binary_summer)),
  identical(rownames(deer_summer_det_hist), rownames(deer_summer_prop_det)),
  identical(rownames(deer_summer_det_hist), rownames(wolf_summer_prop_det)))

stopifnot(
  identical(dim(deer_sum_det_matrix), dim(wolf_sum_det_matrix)),
  nrow(occu_site_covs) == nrow(deer_sum_det_matrix),
  identical(dim(as.matrix(effort_biweekly)), dim(deer_sum_det_matrix)))

###################################
# 2.0 Prepare data
###################################

# Detection covariates need to match dimensions of detection histories
n_surveys <- ncol(deer_summer_det_hist)
n_sites <- nrow(deer_summer_det_hist)

trail_strength_matrix <- matrix(trail_strength[,1], nrow = n_sites, ncol = n_surveys)
deer_binary_summer_matrix <- matrix(deer_binary_summer[,1], nrow = n_sites, ncol = n_surveys)
wolf_binary_summer_matrix <- matrix(wolf_binary_summer[,1], nrow = n_sites, ncol = n_surveys)
deer_summer_prop_det_matrix <- matrix(deer_summer_prop_det[,1], nrow = n_sites, ncol = n_surveys)
wolf_summer_prop_det_matrix <- matrix(wolf_summer_prop_det[,1], nrow = n_sites, ncol = n_surveys)

# # Create OccuMulti Frame
# Create detection history list 
ylist <- list(deer=deer_sum_det_matrix, wolf=wolf_sum_det_matrix)
lapply(ylist, head)

# List of detection covariates 
det_list <- list(trail=trail_strength_matrix, effort = effort_biweekly, 
                 deer_binary=deer_binary_summer_matrix, wolf_binary=wolf_binary_summer_matrix,
                 deer=deer_summer_prop_det_matrix, wolf=wolf_summer_prop_det_matrix )
lapply(det_list, head)

# Construct an unmarkedFrameOccuMulti object
msom_data <- unmarkedFrameOccuMulti(y=ylist, 
                                    siteCovs=occu_site_covs, 
                                    obsCovs=det_list)

summary(msom_data)

###################################
# 3.0 Test GOF for global model
###################################

# Using deer/wolf proportion detected instead of deer/wolf binary (deer/wolf prop is more complex)

# # Define global model for summer

# Detection structure
DW_det_global <- c(
  "~ effort + trail + wolf",
  "~ effort + trail + deer")

# Full marginal occupancy formulas
deer_state_global <- paste(
  "~ coast_decay_z + elevation_z + mean_tri_500_z +",
  "old_growth_prop_buff_z + recent_harvest_prop_buff_z +",
  "second_growth_prop_buff_z + roads_decay_z")

wolf_state_global <- paste(
  "~ coast_decay_z + elevation_z + mean_tri_500_z +",
  "old_growth_prop_buff_z + recent_harvest_prop_buff_z +",
  "second_growth_prop_buff_z + roads_decay_z")

# Full pairwise interaction formula
interaction_state_global <- paste(
  "~ coast_decay_z + elevation_z + mean_tri_500_z +",
  "old_growth_prop_buff_z + recent_harvest_prop_buff_z +",
  "second_growth_prop_buff_z + roads_decay_z")

# Fit global model
summer_global <- occuMulti(
  detformulas = DW_det_global,
  stateformulas = c(
    deer_state_global,
    wolf_state_global,
    interaction_state_global),
  maxOrder = 2,
  data = msom_data)

summary(summer_global)

# GOF function (raw detection-history cell-level SSE, Pearson chi-square, and Freeman-Tukey)
  # Adapted from Appendix 2 Code from Andrade-Ponce, G.P., Mandujano, S., Dáttilo, W., Farías-González, V., Jiménez, J., Velásquez-C, K. and Zavaleta, A. (2022),
  # The original code was modified for an unmarked::occuMulti object by matching observed, fitted, and residual matrices by species.

fitstats_GOF <- function(fm) {
  
  # Observed detection histories for each species
  # Patricio et al. used: observed <- getY(fm@data)
  # For occuMulti, ylist preserves species-specific detection matrices.
  observed <- fm@data@ylist
  
  # Fitted values and residuals
  expected <- fitted(fm)
  resids <- residuals(fm)
  
  # Match species names across observed, expected, and residual outputs
  species_to_use <- Reduce(
    intersect,
    list(names(observed), names(expected), names(resids))
  )
  
  # Sum of squared errors
  sse <- sum(sapply(species_to_use, FUN = function(x) {
    sum(resids[[x]]^2, na.rm = TRUE)
  }))
  
  # Pearson chi-square
  chisq <- sum(sapply(species_to_use, FUN = function(x) {
    sum((observed[[x]] - expected[[x]])^2 / expected[[x]], na.rm = TRUE)
  }))
  
  # Freeman-Tukey discrepancy
  freeTuke <- sum(sapply(species_to_use, FUN = function(x) {
    sum((sqrt(observed[[x]]) - sqrt(expected[[x]]))^2, na.rm = TRUE)
  }))
  
  # Return statistics
  out <- c(SSE = sse, Chisq = chisq, freemanTukey = freeTuke)
  
  return(out)}

# Parametric bootstrap GOF test
pb_file <- "outputs/summer_mod/pb_summer_global_1000_2026-07-30.rds" 

if (file.exists(pb_file)) {
  pb_summer_global <- readRDS(pb_file)
} else {
  set.seed(123)
  
  pb_summer_global <- parboot(
    object = summer_global,
    statistic = fitstats_GOF,
    nsim = 1000,
    parallel = FALSE)
  
  saveRDS(pb_summer_global, pb_file)}

pb_summer_global

# Calculate c-hat from chi squared result
summer_global_chat <- pb_summer_global@t0["Chisq"] / mean(pb_summer_global@t.star[, "Chisq"])
summer_global_chat #1.281104 

# # GOF Plots 
# Extract observed and simulated statistics
t0_summer <- pb_summer_global@t0
tstar_summer <- as.data.frame(pb_summer_global@t.star)
colnames(tstar_summer) <- names(t0_summer)

# Calculate bootstrap p-values
pval_summer_sse <- mean(tstar_summer$SSE >= t0_summer["SSE"], na.rm = TRUE)
pval_summer_chisq <- mean(tstar_summer$Chisq >= t0_summer["Chisq"], na.rm = TRUE)
pval_summer_ft <- mean(tstar_summer$freemanTukey >= t0_summer["freemanTukey"], na.rm = TRUE)

# Format p-values for plot titles
pval_summer_sse_label <- formatC(pval_summer_sse, digits = 3, format = "f")
pval_summer_chisq_label <- formatC(pval_summer_chisq, digits = 3, format = "f")
pval_summer_ft_label <- formatC(pval_summer_ft, digits = 3, format = "f")

# Calculate c-hat values
chat_summer_sse <- t0_summer["SSE"] / mean(tstar_summer$SSE, na.rm = TRUE)
chat_summer_chisq <- t0_summer["Chisq"] / mean(tstar_summer$Chisq, na.rm = TRUE)
chat_summer_ft <- t0_summer["freemanTukey"] / mean(tstar_summer$freemanTukey, na.rm = TRUE)

# Format c-hat values for plot titles
chat_summer <- chat_summer_chisq
chat_summer_chisq_label <- formatC(chat_summer_chisq, digits = 2, format = "f")

# Separate dataframes
sse_summer_df <- tibble(t_B = tstar_summer$SSE)
chisq_summer_df <- tibble(t_B = tstar_summer$Chisq)
ft_summer_df <- tibble(t_B = tstar_summer$freemanTukey)

# SSE 
plot_summer_sse <- ggplot(sse_summer_df, aes(x = t_B)) +
  geom_histogram(bins = 35, fill = "grey80", color = "grey40") +
  geom_vline(
    xintercept = t0_summer["SSE"],
    linetype = "dashed",
    linewidth = 1,
    color = "red") +
  theme_classic() +
  labs(
    title = paste0("Multi-species Bootstrapped SSE statistic\np = ", pval_summer_sse_label),
    x = "Simulated statistic",
    y = "Frequency") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 11),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9))

plot_summer_sse

# Chi-square 
plot_summer_chisq <- ggplot(chisq_summer_df, aes(x = t_B)) +
  geom_histogram(bins = 35, fill = "grey80", color = "grey40") +
  geom_vline(
    xintercept = t0_summer["Chisq"],
    linetype = "dashed",
    linewidth = 1,
    color = "red") +
  theme_classic() +
  labs(
    title = paste0(
      "Multi-species Bootstrapped Chisq statistic\np = ",
      pval_summer_chisq_label,
      "; c-hat = ",
      chat_summer_chisq_label),
    x = "Simulated statistic",
    y = "Frequency") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 11),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9))

plot_summer_chisq

# Freeman-Tukey
plot_summer_ft <- ggplot(ft_summer_df, aes(x = t_B)) +
  geom_histogram(bins = 35, fill = "grey80", color = "grey40") +
  geom_vline(
    xintercept = t0_summer["freemanTukey"],
    linetype = "dashed",
    linewidth = 1,
    color = "red") +
  theme_classic() +
  labs(
    title = paste0("Multi-species Bootstrapped Freeman-Tukey statistic\np = ", pval_summer_ft_label),
    x = "Simulated statistic",
    y = "Frequency") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 11),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9))

plot_summer_ft

# # Combine plots with patchwork

plot_summer_sse <- plot_summer_sse +
  scale_x_continuous(expand = expansion(mult = c(0.10, 0.12)))

plot_summer_chisq <- plot_summer_chisq +
  scale_x_continuous(
    expand = expansion(mult = c(0.05, 0.20))) +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(t = 5, r = 20, b = 5, l = 5))

plot_summer_ft <- plot_summer_ft +
  scale_x_continuous(expand = expansion(mult = c(0.10, 0.12)))

gof_summer_combined <- plot_summer_sse + plot_summer_chisq + plot_summer_ft +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(face = "bold", size = 12)) & 
  theme(
    # Axis values / tick labels
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    # Axis labels / titles
    axis.title.x = element_text(size = 16, margin = margin(t = 8)))

gof_summer_combined

# Create summer GOF summary table

summer_t0 <- tibble(
  Statistic = names(pb_summer_global@t0),
  observed_t0 = as.numeric(pb_summer_global@t0))

summer_gof_summary <- pb_summer_global@t.star %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Statistic", values_to = "tB") %>%
  left_join(summer_t0, by = "Statistic") %>%
  group_by(Statistic, observed_t0) %>%
  summarise(
    mean_tB = mean(tB),
    median_tB = median(tB),
    sd_tB = sd(tB),
    max_tB = max(tB),
    p_value = mean(tB > observed_t0),
    .groups = "drop") %>%
  mutate(
    season = "Summer",
    c_hat_mean = observed_t0 / mean_tB,
    c_hat_median = observed_t0 / median_tB,
    Statistic = case_when(
      Statistic == "Chisq" ~ "Pearson chi-square",
      Statistic == "freemanTukey" ~ "Freeman-Tukey",
      TRUE ~ Statistic)) %>%
  select(
    season, Statistic, observed_t0, mean_tB, median_tB,
    sd_tB, max_tB, p_value, c_hat_mean, c_hat_median)

summer_gof_summary

write.csv(summer_gof_summary,
          "outputs/summer_mod/MSOM_summer_GOF_summary_table.csv",
          row.names = FALSE)

###################################
# 4.0 Detection model selection 
###################################

# null model 
d0 <- occuMulti(
  detformulas = c('~1','~1'),
  stateformulas = c('~1','~1','~1'),
  maxOrder = 2,
  data = msom_data)

summary(d0)

# effort only 
d1 <- occuMulti(
  detformulas = c('~effort','~effort'),
  stateformulas = c('~1','~1','~1'),
  maxOrder = 2,
  data = msom_data)

summary(d1)

# trail strength only 
d2 <- occuMulti(
  detformulas = c('~trail','~trail'),
  stateformulas = c('~1','~1','~1'),
  maxOrder = 2,
  data = msom_data)

summary(d2)

# binary only
d3 <- occuMulti(
  detformulas = c('~wolf_binary','~deer_binary'),
  stateformulas = c('~1','~1','~1'),
  maxOrder = 2,
  data = msom_data)

summary(d3)

# detection history only
d4 <- occuMulti(
  detformulas = c('~wolf','~deer'),
  stateformulas = c('~1','~1','~1'),
  maxOrder = 2,
  data = msom_data)

summary(d4)

# effort + trail
d5 <- occuMulti(
  detformulas = c('~effort + trail','~effort + trail'),
  stateformulas = c('~1','~1','~1'),
  maxOrder = 2,
  data = msom_data)

summary(d5)

# # effort + binary
d6 <- occuMulti(
  detformulas = c('~effort + wolf_binary','~effort + deer_binary'),
  stateformulas = c('~1','~1','~1'),
  maxOrder = 2,
  data = msom_data)

summary(d6)

# effort + detection history
d7 <- occuMulti(
  detformulas = c('~effort + wolf','~effort + deer'),
  stateformulas = c('~1','~1','~1'),
  maxOrder = 2,
  data = msom_data)

summary(d7)

# # trail strength + binary
d8 <- occuMulti(
  detformulas = c('~trail + wolf_binary','~trail + deer_binary'),
  stateformulas = c('~1','~1','~1'),
  maxOrder = 2,
  data = msom_data)

summary(d8)

# # trail strength + detection history
d9 <- occuMulti(
  detformulas = c('~trail + wolf','~trail + deer'),
  stateformulas = c('~1','~1','~1'),
  maxOrder = 2,
  data = msom_data)

summary(d9)

# # effort + trail + binary
d10 <- occuMulti(
  detformulas = c('~effort + trail + wolf_binary',
                  '~effort + trail + deer_binary'),
  stateformulas = c('~1','~1','~1'),
  maxOrder = 2,
  data = msom_data)

summary(d10)

# # effort + trail + detection history
d11 <- occuMulti(
  detformulas = c('~effort + trail + wolf',
                  '~effort + trail + deer'),
  stateformulas = c('~1','~1','~1'),
  maxOrder = 2,
  data = msom_data)

summary(d11)

# Global-model c-hat for summer MSOM (to be used in model selection)
chat_summer <- as.numeric(summer_global_chat)

# Number of site-years
n_siteyears <- nrow(getY(msom_data))

# Create fitlist for model selection table 
det.mods <- list(
  d0 = d0, d1 = d1, d2 = d2, d3 = d3, d4 = d4, d5 = d5,
  d6 = d6, d7 = d7, d8 = d8, d9 = d9, d10 = d10, d11 = d11)

# Run QAICc model selection
det_qaicc_table <- aictab(
  cand.set = det.mods,
  modnames = names(det.mods),
  second.ord = TRUE,
  nobs = n_siteyears,
  c.hat = chat_summer)

det_qaicc_table

# Save AIC table
write.csv(as.data.frame(det_qaicc_table),"outputs/summer_mod/MSOM_biweekly_summer_det_QAICc_table.csv",row.names = FALSE)

# Define detection function
DW_det <-c('~effort + trail + wolf','~effort + trail + deer')

###################################
# 5.0 Marginal Occupancy Model Selection 
###################################

# # Deer

# null model 
deer0 <- occuMulti(
  detformulas = DW_det,
  stateformulas = c('~1','~1','~1'),
  maxOrder = 2,
  data = msom_data)

summary(deer0)

# environment only 
deer1 <- occuMulti(
  detformulas = DW_det,
  stateformulas = c('~coast_decay_z + elevation_z + mean_tri_500_z + old_growth_prop_buff_z','~1','~1'),
  maxOrder = 2,
  data = msom_data)

summary(deer1)

# forestry only 
deer2 <- occuMulti(
  detformulas = DW_det,
  stateformulas = c('~recent_harvest_prop_buff_z + second_growth_prop_buff_z + old_growth_prop_buff_z + roads_decay_z','~1','~1'),
  maxOrder = 2,
  data = msom_data)

summary(deer2)

# environment + roads
deer3 <- occuMulti(
  detformulas = DW_det,
  stateformulas = c('~coast_decay_z + elevation_z + mean_tri_500_z + old_growth_prop_buff_z + roads_decay_z','~1','~1'),
  maxOrder = 2,
  data = msom_data)

summary(deer3)

# environment + forestry
deer4 <- occuMulti(
  detformulas = DW_det,
  stateformulas = c('~coast_decay_z + elevation_z + mean_tri_500_z + old_growth_prop_buff_z + recent_harvest_prop_buff_z + second_growth_prop_buff_z + roads_decay_z','~1','~1'),
  maxOrder = 2,
  data = msom_data)

summary(deer4)

# create fitlist for model selection table 
deer.mods <- list(deer0 = deer0, deer1 = deer1, deer2 = deer2, deer3 = deer3, deer4 = deer4)

# Run QAICc model selection
deer_qaicc_table <- aictab(
  cand.set = deer.mods,
  modnames = names(deer.mods),
  second.ord = TRUE,
  nobs = n_siteyears,
  c.hat = chat_summer)

deer_qaicc_table

# Save AIC table
write.csv(as.data.frame(deer_qaicc_table), "outputs/summer_mod/MSOM_biweekly_summer_deer_QAICc_table.csv", row.names = FALSE)

# # Wolf

# null model 
wolf0 <- occuMulti(
  detformulas = DW_det,
  stateformulas = c('~1','~1','~1'),
  maxOrder = 2,
  data = msom_data)

summary(wolf0)

# environment only 
wolf1 <- occuMulti(
  detformulas = DW_det,
  stateformulas = c('~1','~coast_decay_z + elevation_z + mean_tri_500_z + old_growth_prop_buff_z','~1'),
  maxOrder = 2,
  data = msom_data)

summary(wolf1)

# forestry only 
wolf2 <- occuMulti(
  detformulas = DW_det,
  stateformulas = c('~1','~recent_harvest_prop_buff_z + second_growth_prop_buff_z + old_growth_prop_buff_z + roads_decay_z','~1'),
  maxOrder = 2,
  data = msom_data)

summary(wolf2)

# environment + roads
wolf3 <- occuMulti(
  detformulas = DW_det,
  stateformulas = c('~1','~coast_decay_z + elevation_z + mean_tri_500_z + old_growth_prop_buff_z + roads_decay_z','~1'),
  maxOrder = 2,
  data = msom_data)

summary(wolf3)

# environment + forestry 
wolf4 <- occuMulti(
  detformulas = DW_det,
  stateformulas = c('~1','~recent_harvest_prop_buff_z + second_growth_prop_buff_z + old_growth_prop_buff_z + roads_decay_z + coast_decay_z + elevation_z + mean_tri_500_z','~1'),
  maxOrder = 2,
  data = msom_data)

summary(wolf4)

# Create fitlist
wolf.mods <- list(wolf0 = wolf0, wolf1 = wolf1, wolf2 = wolf2, wolf3 = wolf3, wolf4 = wolf4)

wolf_qaicc_table <- aictab(
  cand.set = wolf.mods,
  modnames = names(wolf.mods),
  second.ord = TRUE,
  nobs = n_siteyears,
  c.hat = chat_summer)

wolf_qaicc_table

# Save AIC table
write.csv(as.data.frame(wolf_qaicc_table),"outputs/summer_mod/MSOM_biweekly_summer_wolf_QAICc_table.csv",row.names = FALSE)

###################################
# 6.0 Co-Occupancy Model Selection 
###################################

# null model
null <- occuMulti(
  detformulas = DW_det,
  stateformulas = c('~coast_decay_z + elevation_z + mean_tri_500_z + old_growth_prop_buff_z + recent_harvest_prop_buff_z + second_growth_prop_buff_z + roads_decay_z',
                    '~coast_decay_z + elevation_z + mean_tri_500_z + old_growth_prop_buff_z',
                    '~1'),
  maxOrder = 2,
  data = msom_data)

summary(null)

# distance to coast (decay)
coast <- occuMulti(
  detformulas = DW_det,
  stateformulas = c('~coast_decay_z + elevation_z + mean_tri_500_z + old_growth_prop_buff_z + recent_harvest_prop_buff_z + second_growth_prop_buff_z + roads_decay_z',
                    '~coast_decay_z + elevation_z + mean_tri_500_z + old_growth_prop_buff_z',
                    '~coast_decay_z'),
  maxOrder = 2,
  data = msom_data)

summary(coast)

# elevation
elevation <- occuMulti(
  detformulas = DW_det,
  stateformulas = c('~coast_decay_z + elevation_z + mean_tri_500_z + old_growth_prop_buff_z + recent_harvest_prop_buff_z + second_growth_prop_buff_z + roads_decay_z',
                    '~coast_decay_z + elevation_z + mean_tri_500_z + old_growth_prop_buff_z',
                    '~elevation_z'),
  maxOrder = 2,
  data = msom_data)

summary(elevation)

# terrain ruggedness index (mean over 500 m buffer around site)
tri <- occuMulti(
  detformulas = DW_det,
  stateformulas = c('~coast_decay_z + elevation_z + mean_tri_500_z + old_growth_prop_buff_z + recent_harvest_prop_buff_z + second_growth_prop_buff_z + roads_decay_z',
                    '~coast_decay_z + elevation_z + mean_tri_500_z + old_growth_prop_buff_z',
                    '~mean_tri_500_z'),
  maxOrder = 2,
  data = msom_data)

summary(tri)

# proportion of recent harvest (500 m buffer around site)
recentharvest <- occuMulti(
  detformulas = DW_det,
  stateformulas = c('~coast_decay_z + elevation_z + mean_tri_500_z + old_growth_prop_buff_z + recent_harvest_prop_buff_z + second_growth_prop_buff_z + roads_decay_z',
                    '~coast_decay_z + elevation_z + mean_tri_500_z + old_growth_prop_buff_z',
                    '~recent_harvest_prop_buff_z'),
  maxOrder = 2,
  data = msom_data)

summary(recentharvest)

# proportion of second growth (500 m buffer around site)
secondgrowth <- occuMulti(
  detformulas = DW_det,
  stateformulas = c('~coast_decay_z + elevation_z + mean_tri_500_z + old_growth_prop_buff_z + recent_harvest_prop_buff_z + second_growth_prop_buff_z + roads_decay_z',
                    '~coast_decay_z + elevation_z + mean_tri_500_z + old_growth_prop_buff_z',
                    '~second_growth_prop_buff_z'),
  maxOrder = 2,
  data = msom_data)

summary(secondgrowth)

# proportion of old growth (500 m buffer around site)
oldgrowth <- occuMulti(
  detformulas = DW_det,
  stateformulas = c('~coast_decay_z + elevation_z + mean_tri_500_z + old_growth_prop_buff_z + recent_harvest_prop_buff_z + second_growth_prop_buff_z + roads_decay_z',
                    '~coast_decay_z + elevation_z + mean_tri_500_z + old_growth_prop_buff_z',
                    '~old_growth_prop_buff_z'),
  maxOrder = 2,
  data = msom_data)

summary(oldgrowth)

# distance to nearest road (decay)
roads <- occuMulti(
  detformulas = DW_det,
  stateformulas = c('~coast_decay_z + elevation_z + mean_tri_500_z + old_growth_prop_buff_z + recent_harvest_prop_buff_z + second_growth_prop_buff_z + roads_decay_z',
                    '~coast_decay_z + elevation_z + mean_tri_500_z + old_growth_prop_buff_z',
                    '~roads_decay_z'),
  maxOrder = 2,
  data = msom_data)

summary(roads)

# Create fitList
cond.mods <- list(null = null, coast = coast, elevation = elevation, tri = tri, recentharvest = recentharvest,
                  secondgrowth = secondgrowth, oldgrowth = oldgrowth, roads = roads)

# Run model selection 
cond_qaicc_table <- aictab(
  cand.set = cond.mods,
  modnames = names(cond.mods),
  second.ord = TRUE,
  nobs = n_siteyears,
  c.hat = chat_summer)

cond_qaicc_table

# Save AICc table 
write.csv(as.data.frame(cond_qaicc_table),"outputs/summer_mod/MSOM_biweekly_summer_cond_QAICc_table.csv",row.names = FALSE)

summary(null) #top model

##################################
# Adjust summary results with c-hat-adjusted SEs/CIs (using global model c-hat)
##################################

# Extract estimates and original SEs
est <- coef(null)
se_original <- sqrt(diag(vcov(null)))

# Adjust SEs using sqrt(c-hat)
se_adjusted <- se_original * sqrt(summer_global_chat)

# Recalculate z and p values
z_adjusted <- est / se_adjusted
p_adjusted <- 2 * pnorm(abs(z_adjusted), lower.tail = FALSE)

# Create coefficient table similar to summary(null)
summary_chat <- cbind(
  Estimate = est,
  SE = se_adjusted,
  z = z_adjusted,
  `P(>|z|)` = p_adjusted)

# Split occupancy and detection parameters
state_rows <- grepl("^psi", rownames(summary_chat))
det_rows <- grepl("^p", rownames(summary_chat))

# # Print adjusted summary

cat("\nCall:\n")
print(null@call)

cat("\nOverdispersion adjustment:\n")
cat("c-hat =", summer_global_chat, "\n")
cat("SE multiplier = sqrt(c-hat) =", sqrt(summer_global_chat), "\n\n")

cat("State model coefficients:\n")
printCoefmat(
  summary_chat[state_rows, , drop = FALSE],
  digits = 3,
  signif.stars = TRUE,
  P.values = TRUE,
  has.Pvalue = TRUE)

cat("\nDetection model coefficients:\n")
printCoefmat(
  summary_chat[det_rows, , drop = FALSE],
  digits = 3,
  signif.stars = TRUE,
  P.values = TRUE,
  has.Pvalue = TRUE)

summary_chat_df <- data.frame(
  parameter = rownames(summary_chat),
  Estimate = summary_chat[, "Estimate"],
  SE_adjusted = summary_chat[, "SE"],
  z_adjusted = summary_chat[, "z"],
  p_adjusted = summary_chat[, "P(>|z|)"],
  row.names = NULL)

summary_chat_df

###################################
# 7.0 Assess Top Model GOF (using previous GOF function)
###################################

pb_file <- "outputs/summer_mod/pb_summer_topmod_1000_2026-07-30.rds" 

if (file.exists(pb_file)) {
  pb_topmod <- readRDS(pb_file)
} else {
  set.seed(123)
  
  pb_topmod <- parboot(
    object = null,
    statistic = fitstats_GOF,
    nsim = 1000,
    parallel = FALSE)
  
  saveRDS(pb_topmod, pb_file)}

pb_topmod

# Calculate c-hat from chi squared statistic result 
    # Use global-model c-hat for model selection, coefficient SE inflation,
    # and prediction uncertainty. Top-model c-hat is calculated for reference only.
topmod_chat <- pb_topmod@t0["Chisq"] / mean(pb_topmod@t.star[, "Chisq"])
topmod_chat #1.31

# # GOF Plots 

# Replace with your actual parboot object name
pb <- pb_topmod

# Extract observed and simulated statistics
t0 <- pb@t0
tstar <- as.data.frame(pb@t.star)
colnames(tstar) <- names(t0)

# Calculate bootstrap p-values
pval_sse <- mean(tstar$SSE >= t0["SSE"], na.rm = TRUE)
pval_chisq <- mean(tstar$Chisq >= t0["Chisq"], na.rm = TRUE)
pval_ft <- mean(tstar$freemanTukey >= t0["freemanTukey"], na.rm = TRUE)

# Format p-values for plot titles
pval_sse_label <- formatC(pval_sse, digits = 3, format = "f")
pval_chisq_label <- formatC(pval_chisq, digits = 3, format = "f")
pval_ft_label <- formatC(pval_ft, digits = 3, format = "f")

# Calculate c-hat values
chat_sse <- t0["SSE"] / mean(tstar$SSE, na.rm = TRUE)
chat_chisq <- t0["Chisq"] / mean(tstar$Chisq, na.rm = TRUE)
chat_ft <- t0["freemanTukey"] / mean(tstar$freemanTukey, na.rm = TRUE)

# Format c-hat values for plot titles
chat_chisq_label <- formatC(chat_chisq, digits = 2, format = "f")

# Separate dataframes
sse_df <- tibble(t_B = tstar$SSE)
chisq_df <- tibble(t_B = tstar$Chisq)
ft_df <- tibble(t_B = tstar$freemanTukey)

# SSE 
plot_sse <- ggplot(sse_df, aes(x = t_B)) +
  geom_histogram(bins = 35, fill = "grey80", color = "grey40") +
  geom_vline(
    xintercept = t0["SSE"],
    linetype = "dashed",
    linewidth = 1,
    color = "red") +
  theme_classic() +
  labs(
    title = paste0("Multi-species Bootstrapped SSE statistic\np = ", pval_sse_label),
    x = "Simulated statistic",
    y = "Frequency") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 11),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9))

plot_sse

# Chi-square 
plot_chisq <- ggplot(chisq_df, aes(x = t_B)) +
  geom_histogram(bins = 35, fill = "grey80", color = "grey40") +
  geom_vline(
    xintercept = t0["Chisq"],
    linetype = "dashed",
    linewidth = 1,
    color = "red") +
  theme_classic() +
  labs(
    title = paste0(
      "Multi-species Bootstrapped Chisq statistic\np = ",
      pval_chisq_label,
      "; c-hat = ",
      chat_chisq_label),
    x = "Simulated statistic",
    y = "Frequency") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 11),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9))

plot_chisq

# Freeman-Tukey
plot_ft <- ggplot(ft_df, aes(x = t_B)) +
  geom_histogram(bins = 35, fill = "grey80", color = "grey40") +
  geom_vline(
    xintercept = t0["freemanTukey"],
    linetype = "dashed",
    linewidth = 1,
    color = "red") +
  theme_classic() +
  labs(
    title = paste0("Multi-species Bootstrapped Freeman-Tukey statistic\np = ", pval_ft_label),
    x = "Simulated statistic",
    y = "Frequency") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 11),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9))

plot_ft

# # Combine plots with patchwork
plot_sse <- plot_sse +
  scale_x_continuous(
    expand = expansion(mult = c(0.10, 0.12)))

plot_chisq <- plot_chisq +
  scale_x_continuous(
    expand = expansion(mult = c(0.05, 0.20))) +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(t = 5, r = 20, b = 5, l = 5))

plot_ft <- plot_ft +
  scale_x_continuous(
    expand = expansion(mult = c(0.10, 0.12)))

gof_summer_topmod_combined <-
  plot_sse + plot_chisq + plot_ft +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 16, margin = margin(t = 8)))

gof_summer_topmod_combined

###################################
# 8.0 Overdispersion-adjusted predictions using AICcmodavg
###################################

# Single-model candidate set for top-model predictions
top_mod <- list(null = null)

# State covariates for summary occupancy predictions
state_cov_names <- c(
  "coast_decay_z",
  "elevation_z",
  "mean_tri_500_z",
  "old_growth_prop_buff_z",
  "recent_harvest_prop_buff_z",
  "second_growth_prop_buff_z",
  "roads_decay_z")

state_newdata <- occu_site_covs[, state_cov_names]

# Detection covariates for summary detection predictions
trail_levels <- c("FAINT", "MODERATE", "STRONG")

det_newdata <- data.frame(
  effort = as.numeric(as.matrix(det_list$effort)),
  trail = factor(as.vector(as.matrix(det_list$trail)), levels = trail_levels),
  wolf = as.numeric(as.matrix(det_list$wolf)),
  deer = as.numeric(as.matrix(det_list$deer)))

det_newdata <- det_newdata[complete.cases(det_newdata), ]

###################################
# 9.0 Plot mean predicted marginal occupancy and detection probabilities
#     with averaged c-hat-adjusted lower and upper confidence limits
###################################

# # Calculate mean marginal occupancy probability

# Deer marginal occupancy
null_mo_deer_raw <- modavgPred(
  cand.set = top_mod,
  modnames = names(top_mod),
  newdata = state_newdata,
  parm.type = "psi",
  type = "response",
  c.hat = summer_global_chat,
  nobs = n_siteyears,
  uncond.se = "revised",
  species = "deer")

null_mo_deer <- data.frame(
  Predicted = as.numeric(null_mo_deer_raw$mod.avg.pred),
  SE = as.numeric(null_mo_deer_raw$uncond.se),
  lower = as.numeric(null_mo_deer_raw$lower.CL),
  upper = as.numeric(null_mo_deer_raw$upper.CL))

null_d_mo <- mean(null_mo_deer$Predicted, na.rm = TRUE)
null_d_lower <- mean(pmax(null_mo_deer$lower, 0), na.rm = TRUE)
null_d_upper <- mean(pmin(null_mo_deer$upper, 1), na.rm = TRUE)

# Wolf marginal occupancy
null_mo_wolf_raw <- modavgPred(
  cand.set = top_mod,
  modnames = names(top_mod),
  newdata = state_newdata,
  parm.type = "psi",
  type = "response",
  c.hat = summer_global_chat,
  nobs = n_siteyears,
  uncond.se = "revised",
  species = "wolf")

null_mo_wolf <- data.frame(
  Predicted = as.numeric(null_mo_wolf_raw$mod.avg.pred),
  SE = as.numeric(null_mo_wolf_raw$uncond.se),
  lower = as.numeric(null_mo_wolf_raw$lower.CL),
  upper = as.numeric(null_mo_wolf_raw$upper.CL))

null_w_mo <- mean(null_mo_wolf$Predicted, na.rm = TRUE)
null_w_lower <- mean(pmax(null_mo_wolf$lower, 0), na.rm = TRUE)
null_w_upper <- mean(pmin(null_mo_wolf$upper, 1), na.rm = TRUE)

# Save mean predicted values and mean lower/upper limits
null_occ_pred <- c(null_d_mo, null_w_mo)
null_occ_lower <- c(null_d_lower, null_w_lower)
null_occ_upper <- c(null_d_upper, null_w_upper)

# Combine into a dataframe
null_occ_pred.df <- data.frame(
  species = c("Deer", "Wolf"),
  pred = null_occ_pred,
  lower = null_occ_lower,
  upper = null_occ_upper)

# Plot predicted marginal occupancy probabilities 
p_occ <- ggplot(
  data = null_occ_pred.df, 
  aes(x = species, y = pred)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = pmax(lower, 0), ymax = pmin(upper, 1)), width = 0.1, color = "black") +
  scale_y_continuous(
    breaks = seq(0, 1, 0.1),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  labs(x = "Species", y = "Predicted Marginal Occupancy") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15))

# # Calculate detection probability
# Deer detection
null_det_deer_raw <- modavgPred(
  cand.set = top_mod,
  modnames = names(top_mod),
  newdata = det_newdata,
  parm.type = "detect",
  type = "response",
  c.hat = summer_global_chat,
  nobs = n_siteyears,
  uncond.se = "revised",
  species = "deer")

null_det_deer <- data.frame(
  Predicted = as.numeric(null_det_deer_raw$mod.avg.pred),
  SE = as.numeric(null_det_deer_raw$uncond.se),
  lower = as.numeric(null_det_deer_raw$lower.CL),
  upper = as.numeric(null_det_deer_raw$upper.CL))

null_d_det <- mean(null_det_deer$Predicted, na.rm = TRUE)
null_d_det_lower <- mean(pmax(null_det_deer$lower, 0), na.rm = TRUE)
null_d_det_upper <- mean(pmin(null_det_deer$upper, 1), na.rm = TRUE)

# Wolf detection
null_det_wolf_raw <- modavgPred(
  cand.set = top_mod,
  modnames = names(top_mod),
  newdata = det_newdata,
  parm.type = "detect",
  type = "response",
  c.hat = summer_global_chat,
  nobs = n_siteyears,
  uncond.se = "revised",
  species = "wolf")

null_det_wolf <- data.frame(
  Predicted = as.numeric(null_det_wolf_raw$mod.avg.pred),
  SE = as.numeric(null_det_wolf_raw$uncond.se),
  lower = as.numeric(null_det_wolf_raw$lower.CL),
  upper = as.numeric(null_det_wolf_raw$upper.CL))

null_w_det <- mean(null_det_wolf$Predicted, na.rm = TRUE)
null_w_det_lower <- mean(pmax(null_det_wolf$lower, 0), na.rm = TRUE)
null_w_det_upper <- mean(pmin(null_det_wolf$upper, 1), na.rm = TRUE)

# Save predicted values
det_pred <- c(null_d_det, null_w_det)
det_lower <- c(null_d_det_lower, null_w_det_lower)
det_upper <- c(null_d_det_upper, null_w_det_upper)

# Combine into dataframe 
det_pred.df <- data.frame(
  species = c("Deer", "Wolf"),
  pred = c(null_d_det, null_w_det),
  lower = c(null_d_det_lower, null_w_det_lower),
  upper = c(null_d_det_upper, null_w_det_upper))

# Plot predicted detection probabilities
p_det <- ggplot(
  data = det_pred.df,
  aes(x = species, y = pred)) +
  geom_point(size = 4, color = "black") +
  geom_errorbar(aes(ymin = pmax(lower, 0), ymax = pmin(upper, 1)), width = 0.1, color = "black") +
  scale_y_continuous(
    breaks = seq(0, 1, 0.1),
    limits = c(0, 1),
    expand = c(0, 0)) +
  labs(x = "Species", y = "Predicted Detection Probability") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15))

# # Combine marginal occupancy and detection plots
p_occ + p_det

###################################
# 10.0 Predict detection probabilities with c-hat-adjusted SEs/CIs
###################################

# # Create dataframes for detection covariates 
trail_levels <- c("FAINT", "MODERATE", "STRONG")

effort_vals <- as.numeric(as.matrix(det_list$effort))
wolf_vals <- as.numeric(as.matrix(det_list$wolf))
deer_vals <- as.numeric(as.matrix(det_list$deer))

mean_effort <- mean(effort_vals, na.rm = TRUE)
mean_wolf <- mean(wolf_vals, na.rm = TRUE)
mean_deer <- mean(deer_vals, na.rm = TRUE)

# Build a dataframe for trail strength
trail_df <- data.frame(
  trail = factor(trail_levels, levels = trail_levels),
  effort = mean_effort,
  wolf = mean_wolf,
  deer = mean_deer)

# for wolf prop_det effects on deer detection
deer_wolf_df <- data.frame(
  trail = factor("FAINT", levels = trail_levels),
  effort = mean_effort,
  wolf = seq(min(wolf_vals, na.rm = TRUE), max(wolf_vals, na.rm = TRUE), length.out = 100),
  deer = mean_deer)

# for deer prop_det effects on wolf detection
wolf_deer_df <- data.frame(
  trail = factor("FAINT", levels = trail_levels),
  effort = mean_effort,
  wolf = mean_wolf,
  deer = seq(min(deer_vals, na.rm = TRUE), max(deer_vals, na.rm = TRUE), length.out = 100))

# for effort effect on deer detection
deer_effort_df <- data.frame(
  trail = factor("FAINT", levels = trail_levels),
  effort = seq(min(effort_vals, na.rm = TRUE), max(effort_vals, na.rm = TRUE), length.out = 100),
  wolf = mean_wolf,
  deer = mean_deer)

# for effort effect on wolf detection
wolf_effort_df <- data.frame(
  trail = factor("FAINT", levels = trail_levels),
  effort = seq(min(effort_vals, na.rm = TRUE), max(effort_vals, na.rm = TRUE), length.out = 100),
  wolf = mean_wolf,
  deer = mean_deer)

# Predictions 
deer_trail_raw <- modavgPred(
  cand.set = top_mod,
  modnames = names(top_mod),
  newdata = trail_df,
  parm.type = "detect",
  type = "response",
  c.hat = summer_global_chat,
  nobs = n_siteyears,
  uncond.se = "revised",
  species = "deer")

deer_trail_pred <- cbind(
  trail_df,
  data.frame(
    Predicted = as.numeric(deer_trail_raw$mod.avg.pred),
    SE = as.numeric(deer_trail_raw$uncond.se),
    lower = as.numeric(deer_trail_raw$lower.CL),
    upper = as.numeric(deer_trail_raw$upper.CL)))

wolf_trail_raw <- modavgPred(
  cand.set = top_mod,
  modnames = names(top_mod),
  newdata = trail_df,
  parm.type = "detect",
  type = "response",
  c.hat = summer_global_chat,
  nobs = n_siteyears,
  uncond.se = "revised",
  species = "wolf")

wolf_trail_pred <- cbind(
  trail_df,
  data.frame(
    Predicted = as.numeric(wolf_trail_raw$mod.avg.pred),
    SE = as.numeric(wolf_trail_raw$uncond.se),
    lower = as.numeric(wolf_trail_raw$lower.CL),
    upper = as.numeric(wolf_trail_raw$upper.CL)))

deer_wolf_raw <- modavgPred(
  cand.set = top_mod,
  modnames = names(top_mod),
  newdata = deer_wolf_df,
  parm.type = "detect",
  type = "response",
  c.hat = summer_global_chat,
  nobs = n_siteyears,
  uncond.se = "revised",
  species = "deer")

deer_wolf_pred <- cbind(
  deer_wolf_df,
  data.frame(
    Predicted = as.numeric(deer_wolf_raw$mod.avg.pred),
    SE = as.numeric(deer_wolf_raw$uncond.se),
    lower = as.numeric(deer_wolf_raw$lower.CL),
    upper = as.numeric(deer_wolf_raw$upper.CL)))

wolf_deer_raw <- modavgPred(
  cand.set = top_mod,
  modnames = names(top_mod),
  newdata = wolf_deer_df,
  parm.type = "detect",
  type = "response",
  c.hat = summer_global_chat,
  nobs = n_siteyears,
  uncond.se = "revised",
  species = "wolf")

wolf_deer_pred <- cbind(
  wolf_deer_df,
  data.frame(
    Predicted = as.numeric(wolf_deer_raw$mod.avg.pred),
    SE = as.numeric(wolf_deer_raw$uncond.se),
    lower = as.numeric(wolf_deer_raw$lower.CL),
    upper = as.numeric(wolf_deer_raw$upper.CL)))

deer_effort_raw <- modavgPred(
  cand.set = top_mod,
  modnames = names(top_mod),
  newdata = deer_effort_df,
  parm.type = "detect",
  type = "response",
  c.hat = summer_global_chat,
  nobs = n_siteyears,
  uncond.se = "revised",
  species = "deer")

deer_effort_pred <- cbind(
  deer_effort_df,
  data.frame(
    Predicted = as.numeric(deer_effort_raw$mod.avg.pred),
    SE = as.numeric(deer_effort_raw$uncond.se),
    lower = as.numeric(deer_effort_raw$lower.CL),
    upper = as.numeric(deer_effort_raw$upper.CL)))

wolf_effort_raw <- modavgPred(
  cand.set = top_mod,
  modnames = names(top_mod),
  newdata = wolf_effort_df,
  parm.type = "detect",
  type = "response",
  c.hat = summer_global_chat,
  nobs = n_siteyears,
  uncond.se = "revised",
  species = "wolf")

wolf_effort_pred <- cbind(
  wolf_effort_df,
  data.frame(
    Predicted = as.numeric(wolf_effort_raw$mod.avg.pred),
    SE = as.numeric(wolf_effort_raw$uncond.se),
    lower = as.numeric(wolf_effort_raw$lower.CL),
    upper = as.numeric(wolf_effort_raw$upper.CL)))

# # Plots

# Deer detection vs trail strength
deer_trail_plot <- ggplot(deer_trail_pred, aes(x = trail, y = Predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = pmax(lower, 0), ymax = pmin(upper, 1)), width = 0.2) +
  coord_cartesian(ylim = c(0, 0.6)) +
  labs(x = "Trail Strength", y = "Deer Detection Probability") +
  theme_classic() +
  theme(axis.title.y = element_text(size = 16, angle = 90, margin = margin(r = 8)))

# Deer detection vs wolf detected
deer_wolf_plot <- ggplot(deer_wolf_pred, aes(x = wolf, y = Predicted)) +
  geom_line(linewidth = 1, colour="goldenrod3") +
  geom_ribbon(aes(ymin = pmax(lower, 0), ymax = pmin(upper, 1)), fill="goldenrod3",alpha = 0.30, colour=NA) +
  coord_cartesian(ylim = c(0, 0.6)) +
  labs(x = "Proportion of sampling occasions with a wolf detection", y = "Deer Detection Probability") +
  theme_classic() +
  theme(axis.title.y = element_blank()) 

# Wolf detection vs trail strength
wolf_trail_plot <- ggplot(wolf_trail_pred, aes(x = trail, y = Predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = pmax(lower, 0), ymax = pmin(upper, 1)), width = 0.2) +
  coord_cartesian(ylim = c(0, 0.6)) +
  labs(x = "Trail Strength", y = "Wolf Detection Probability") +
  theme_classic() +
  theme(axis.title.y = element_text(size = 16, angle = 90, margin = margin(r = 8)))


# Wolf detection vs deer detected
wolf_deer_plot <- ggplot(wolf_deer_pred, aes(x = deer, y = Predicted)) +
  geom_line(linewidth = 1, colour="goldenrod3") +
  geom_ribbon(aes(ymin = pmax(lower, 0), ymax = pmin(upper, 1)), fill="goldenrod3",alpha = 0.30, colour=NA) +
  coord_cartesian(ylim = c(0, 0.6)) +
  labs(x = "Proportion of sampling occasions with a deer detection", y = "Wolf Detection Probability") +
  theme_classic() +
  theme(axis.title.y = element_blank())

# Deer detection vs effort
deer_effort_plot <- ggplot(deer_effort_pred, aes(x = effort, y = Predicted)) +
  geom_line(linewidth = 1, colour="royalblue4") +
  geom_ribbon(aes(ymin = pmax(lower, 0), ymax = pmin(upper, 1)), fill="royalblue4",alpha = 0.30, colour=NA) +
  coord_cartesian(ylim = c(0, 0.6)) +
  labs(
    x = "Effort (number of active camera days per bi-weekly occasion)",
    y = "Deer Detection Probability"
  ) +
  theme_classic() +
  theme(axis.title.y = element_blank())

# Wolf detection vs effort
wolf_effort_plot <- ggplot(wolf_effort_pred, aes(x = effort, y = Predicted)) +
  geom_line(linewidth = 1, colour="royalblue4") +
  geom_ribbon(aes(ymin = pmax(lower, 0), ymax = pmin(upper, 1)), fill="royalblue4",alpha = 0.30, colour=NA) +
  coord_cartesian(ylim = c(0, 0.6)) +
  labs(
    x = "Effort (number of active camera days per bi-weekly occasion)",
    y = "Wolf Detection Probability"
  ) +
  theme_classic() +
  theme(axis.title.y = element_blank())

# Combine plots with patchwork 
((deer_trail_plot | deer_wolf_plot | deer_effort_plot) /
    (wolf_trail_plot | wolf_deer_plot | wolf_effort_plot)) & 
  theme(
    # Axis values / tick labels
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    # Axis labels / titles
    axis.title.x = element_text(size = 16, margin = margin(t = 8)))

###################################
# 11.0 Plot Predicted Marginal Occupancy with c-hat-adjusted SEs/CIs
###################################

# # Build dataframes for each environmental covariate (hold all other covs at mean)

# Calculate mean for each covariate (to be held constant)
mean_recent <- mean(occu_site_covs$recent_harvest_prop_buff_z)
mean_second <- mean(occu_site_covs$second_growth_prop_buff_z)
mean_og     <- mean(occu_site_covs$old_growth_prop_buff_z)
mean_roads  <- mean(occu_site_covs$roads_decay_z)
mean_coast  <- mean(occu_site_covs$coast_decay_z)
mean_elev   <- mean(occu_site_covs$elevation_z)
mean_tri    <- mean(occu_site_covs$mean_tri_500_z)

# Recent Harvest dataframe
recent_df <- expand.grid(
  recent_harvest_prop_buff_z = seq(min(occu_site_covs$recent_harvest_prop_buff_z),
                                   max(occu_site_covs$recent_harvest_prop_buff_z),
                                   0.01),
  second_growth_prop_buff_z = mean_second,
  old_growth_prop_buff_z = mean_og,
  roads_decay_z = mean_roads,
  coast_decay_z = mean_coast,
  elevation_z = mean_elev,
  mean_tri_500_z = mean_tri)

# back-transform for x-axis
recent_df$recent_harvest_prop_buff <-
  recent_df$recent_harvest_prop_buff_z *
  sd(orig_site_covs$recent_harvest_prop_buff, na.rm = TRUE) +
  mean(orig_site_covs$recent_harvest_prop_buff, na.rm = TRUE)

# Second growth dataframe
second_df <- expand.grid(
  second_growth_prop_buff_z = seq(min(occu_site_covs$second_growth_prop_buff_z),
                                  max(occu_site_covs$second_growth_prop_buff_z),
                                  0.01),
  recent_harvest_prop_buff_z = mean_recent,
  old_growth_prop_buff_z = mean_og,
  roads_decay_z = mean_roads,
  coast_decay_z = mean_coast,
  elevation_z = mean_elev,
  mean_tri_500_z = mean_tri)

second_df$second_growth_prop_buff <-
  second_df$second_growth_prop_buff_z *
  sd(orig_site_covs$second_growth_prop_buff, na.rm = TRUE) +
  mean(orig_site_covs$second_growth_prop_buff, na.rm = TRUE)

summary(second_df$second_growth_prop_buff_z)
summary(second_df$second_growth_prop_buff)

range(second_df$second_growth_prop_buff_z, na.rm = TRUE)
range(second_df$second_growth_prop_buff, na.rm = TRUE)

# Old growth dataframe
og_df <- expand.grid(
  old_growth_prop_buff_z = seq(min(occu_site_covs$old_growth_prop_buff_z),
                               max(occu_site_covs$old_growth_prop_buff_z),
                               0.01),
  recent_harvest_prop_buff_z = mean_recent,
  second_growth_prop_buff_z = mean_second,
  roads_decay_z = mean_roads,
  coast_decay_z = mean_coast,
  elevation_z = mean_elev,
  mean_tri_500_z = mean_tri)

og_df$old_growth_prop_buff <-
  og_df$old_growth_prop_buff_z *
  sd(orig_site_covs$old_growth_prop_buff, na.rm = TRUE) +
  mean(orig_site_covs$old_growth_prop_buff, na.rm = TRUE)

# Roads dataframe
roads_df <- expand.grid(
  roads_decay_z = seq(min(occu_site_covs$roads_decay_z),
                      max(occu_site_covs$roads_decay_z),
                      0.01),
  recent_harvest_prop_buff_z = mean_recent,
  second_growth_prop_buff_z = mean_second,
  old_growth_prop_buff_z = mean_og,
  coast_decay_z = mean_coast,
  elevation_z = mean_elev,
  mean_tri_500_z = mean_tri)

roads_df$roads_decay <-
  roads_df$roads_decay_z *
  sd(orig_site_covs$roads_decay, na.rm = TRUE) +
  mean(orig_site_covs$roads_decay, na.rm = TRUE)

# Coast 
coast_df <- expand.grid(
  coast_decay_z = seq(min(occu_site_covs$coast_decay_z),
                      max(occu_site_covs$coast_decay_z),
                      0.01),
  elevation_z = mean_elev,
  mean_tri_500_z = mean_tri,
  recent_harvest_prop_buff_z = mean_recent,
  second_growth_prop_buff_z = mean_second,
  old_growth_prop_buff_z = mean_og,
  roads_decay_z = mean_roads)

coast_df$coast_decay <-
  coast_df$coast_decay_z *
  sd(orig_site_covs$coast_decay, na.rm = TRUE) +
  mean(orig_site_covs$coast_decay, na.rm = TRUE)

# Elevation dataframe
elev_df <- expand.grid(
  elevation_z = seq(min(occu_site_covs$elevation_z),
                    max(occu_site_covs$elevation_z),
                    0.01),
  coast_decay_z = mean_coast,
  mean_tri_500_z = mean_tri,
  recent_harvest_prop_buff_z = mean_recent,
  second_growth_prop_buff_z = mean_second,
  old_growth_prop_buff_z = mean_og,
  roads_decay_z = mean_roads)

elev_df$elevation <-
  elev_df$elevation_z *
  sd(orig_site_covs$elevation, na.rm = TRUE) +
  mean(orig_site_covs$elevation, na.rm = TRUE)

# TRI dataframe
tri_df <- expand.grid(
  mean_tri_500_z = seq(min(occu_site_covs$mean_tri_500_z),
                       max(occu_site_covs$mean_tri_500_z),
                       0.01),
  coast_decay_z = mean_coast,
  elevation_z = mean_elev,
  recent_harvest_prop_buff_z = mean_recent,
  second_growth_prop_buff_z = mean_second,
  old_growth_prop_buff_z = mean_og,
  roads_decay_z = mean_roads)

tri_df$mean_tri_500 <-
  tri_df$mean_tri_500_z *
  sd(orig_site_covs$mean_tri_500, na.rm = TRUE) +
  mean(orig_site_covs$mean_tri_500, na.rm = TRUE)

# # Predict deer marginal occupancy probabilities
deer_recent_raw <- modavgPred(
  cand.set = top_mod,
  modnames = names(top_mod),
  newdata = recent_df,
  parm.type = "psi",
  type = "response",
  c.hat = summer_global_chat,
  nobs = n_siteyears,
  uncond.se = "revised",
  species = "deer")

deer_recent_pred <- cbind(
  recent_df,
  data.frame(
    Predicted = as.numeric(deer_recent_raw$mod.avg.pred),
    SE = as.numeric(deer_recent_raw$uncond.se),
    lower = as.numeric(deer_recent_raw$lower.CL),
    upper = as.numeric(deer_recent_raw$upper.CL)))

deer_second_raw <- modavgPred(
  cand.set = top_mod,
  modnames = names(top_mod),
  newdata = second_df,
  parm.type = "psi",
  type = "response",
  c.hat = summer_global_chat,
  nobs = n_siteyears,
  uncond.se = "revised",
  species = "deer")

deer_second_pred <- cbind(
  second_df,
  data.frame(
    Predicted = as.numeric(deer_second_raw$mod.avg.pred),
    SE = as.numeric(deer_second_raw$uncond.se),
    lower = as.numeric(deer_second_raw$lower.CL),
    upper = as.numeric(deer_second_raw$upper.CL)))

deer_og_raw <- modavgPred(
  cand.set = top_mod,
  modnames = names(top_mod),
  newdata = og_df,
  parm.type = "psi",
  type = "response",
  c.hat = summer_global_chat,
  nobs = n_siteyears,
  uncond.se = "revised",
  species = "deer")

deer_og_pred <- cbind(
  og_df,
  data.frame(
    Predicted = as.numeric(deer_og_raw$mod.avg.pred),
    SE = as.numeric(deer_og_raw$uncond.se),
    lower = as.numeric(deer_og_raw$lower.CL),
    upper = as.numeric(deer_og_raw$upper.CL)))

deer_roads_raw <- modavgPred(
  cand.set = top_mod,
  modnames = names(top_mod),
  newdata = roads_df,
  parm.type = "psi",
  type = "response",
  c.hat = summer_global_chat,
  nobs = n_siteyears,
  uncond.se = "revised",
  species = "deer")

deer_roads_pred <- cbind(
  roads_df,
  data.frame(
    Predicted = as.numeric(deer_roads_raw$mod.avg.pred),
    SE = as.numeric(deer_roads_raw$uncond.se),
    lower = as.numeric(deer_roads_raw$lower.CL),
    upper = as.numeric(deer_roads_raw$upper.CL)))

deer_coast_raw <- modavgPred(
  cand.set = top_mod,
  modnames = names(top_mod),
  newdata = coast_df,
  parm.type = "psi",
  type = "response",
  c.hat = summer_global_chat,
  nobs = n_siteyears,
  uncond.se = "revised",
  species = "deer")

deer_coast_pred <- cbind(
  coast_df,
  data.frame(
    Predicted = as.numeric(deer_coast_raw$mod.avg.pred),
    SE = as.numeric(deer_coast_raw$uncond.se),
    lower = as.numeric(deer_coast_raw$lower.CL),
    upper = as.numeric(deer_coast_raw$upper.CL)))

deer_elev_raw <- modavgPred(
  cand.set = top_mod,
  modnames = names(top_mod),
  newdata = elev_df,
  parm.type = "psi",
  type = "response",
  c.hat = summer_global_chat,
  nobs = n_siteyears,
  uncond.se = "revised",
  species = "deer")

deer_elev_pred <- cbind(
  elev_df,
  data.frame(
    Predicted = as.numeric(deer_elev_raw$mod.avg.pred),
    SE = as.numeric(deer_elev_raw$uncond.se),
    lower = as.numeric(deer_elev_raw$lower.CL),
    upper = as.numeric(deer_elev_raw$upper.CL)))

deer_tri_raw <- modavgPred(
  cand.set = top_mod,
  modnames = names(top_mod),
  newdata = tri_df,
  parm.type = "psi",
  type = "response",
  c.hat = summer_global_chat,
  nobs = n_siteyears,
  uncond.se = "revised",
  species = "deer")

deer_tri_pred <- cbind(
  tri_df,
  data.frame(
    Predicted = as.numeric(deer_tri_raw$mod.avg.pred),
    SE = as.numeric(deer_tri_raw$uncond.se),
    lower = as.numeric(deer_tri_raw$lower.CL),
    upper = as.numeric(deer_tri_raw$upper.CL)))

# # Plot deer marginal occupancy probabilities

# Coast
deer_coast_pred_plot <- ggplot(deer_coast_pred, aes(x = coast_decay, y = Predicted)) +
  geom_line(linewidth = 1, color = "lightblue3") +
  geom_ribbon(aes(ymin = pmax(lower, 0), ymax = pmin(upper, 1)), fill = "lightblue3", alpha = 0.45) +
  xlab("Distance from Coast (Decay Rate)") +
  ylab("Deer Occupancy Probability") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic() +
  theme(axis.title.y = element_blank())

# Elevation
deer_elev_pred_plot <- ggplot(deer_elev_pred, aes(x = elevation, y = Predicted)) +
  geom_line(linewidth = 1, color = "thistle3") +
  geom_ribbon(aes(ymin = pmax(lower, 0), ymax = pmin(upper, 1)), fill = "thistle3", alpha = 0.45) +
  xlab("Elevation (m)") +
  ylab("Deer Occupancy Probability") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic() +
  theme(axis.title.y = element_blank())

# TRI
deer_tri_pred_plot <- ggplot(deer_tri_pred, aes(x = mean_tri_500, y = Predicted)) +
  geom_line(linewidth = 1, color = "lightgoldenrod") +
  geom_ribbon(aes(ymin = pmax(lower, 0), ymax = pmin(upper, 1)), fill = "lightgoldenrod", alpha = 0.45) +
  xlab("Mean Terrain Ruggedness Index (500 m buffer)") +
  ylab("Deer Occupancy Probability") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic() + 
  theme(axis.title.y = element_blank())

# Old Growth
deer_og_pred_plot <- ggplot(deer_og_pred, aes(x = old_growth_prop_buff, y = Predicted)) +
  geom_line(linewidth = 1, color = "darkolivegreen4") +
  geom_ribbon(aes(ymin = pmax(lower, 0), ymax = pmin(upper, 1)), fill = "darkolivegreen4", alpha = 0.45) +
  xlab("Proportion of Old Growth (500 m buffer)") +
  ylab("Deer Occupancy Probability") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic() + 
  theme(axis.title.y = element_blank())

# Recent Harvest
deer_recent_pred_plot <- ggplot(deer_recent_pred, aes(x = recent_harvest_prop_buff, y = Predicted)) +
  geom_line(linewidth = 1, color = "brown3") +
  geom_ribbon(aes(ymin = pmax(lower, 0), ymax = pmin(upper, 1)), fill = "brown3", alpha = 0.45) +
  xlab("Proportion of Recent Harvest (500 m buffer)") +
  ylab("Deer Occupancy Probability") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic() + 
  theme(axis.title.y = element_blank())

# Second growth
deer_second_pred_plot <- ggplot(deer_second_pred, aes(x = second_growth_prop_buff, y = Predicted)) +
  geom_line(linewidth = 1, color = "darkorange3") +
  geom_ribbon(aes(ymin = pmax(lower, 0), ymax = pmin(upper, 1)), fill = "darkorange3", alpha = 0.45) +
  xlab("Proportion of Second Growth (500 m buffer)") +
  ylab("Deer Occupancy Probability") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic()  + 
  theme(axis.title.y = element_blank())

# Roads
deer_roads_pred_plot <- ggplot(deer_roads_pred, aes(x = roads_decay, y = Predicted)) +
  geom_line(linewidth = 1, color = "tan") +
  geom_ribbon(aes(ymin = pmax(lower, 0), ymax = pmin(upper, 1)), fill = "tan", alpha = 0.45) +
  xlab("Distance from Nearest Road (Decay Rate)") +
  ylab("Deer Occupancy Probability") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic()  + 
  theme(axis.title.y = element_blank())

# Combine plots with patchwork
(deer_coast_pred_plot + deer_elev_pred_plot + deer_tri_pred_plot) &
  theme(
    # Axis values / tick labels
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    # Axis labels / titles
    axis.title.x = element_text(size = 16, margin = margin(t = 8)))

(deer_recent_pred_plot + deer_second_pred_plot) / (deer_og_pred_plot + deer_roads_pred_plot) & 
  theme(
    # Axis values / tick labels
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    # Axis labels / titles
    axis.title.x = element_text(size = 16, margin = margin(t = 8)))

# # Predict wolf marginal occupancy probabilities
wolf_coast_raw <- modavgPred(
  cand.set = top_mod,
  modnames = names(top_mod),
  newdata = coast_df,
  parm.type = "psi",
  type = "response",
  c.hat = summer_global_chat,
  nobs = n_siteyears,
  uncond.se = "revised",
  species = "wolf")

wolf_coast_pred <- cbind(
  coast_df,
  data.frame(
    Predicted = as.numeric(wolf_coast_raw$mod.avg.pred),
    SE = as.numeric(wolf_coast_raw$uncond.se),
    lower = as.numeric(wolf_coast_raw$lower.CL),
    upper = as.numeric(wolf_coast_raw$upper.CL)))

wolf_elev_raw <- modavgPred(
  cand.set = top_mod,
  modnames = names(top_mod),
  newdata = elev_df,
  parm.type = "psi",
  type = "response",
  c.hat = summer_global_chat,
  nobs = n_siteyears,
  uncond.se = "revised",
  species = "wolf")

wolf_elev_pred <- cbind(
  elev_df,
  data.frame(
    Predicted = as.numeric(wolf_elev_raw$mod.avg.pred),
    SE = as.numeric(wolf_elev_raw$uncond.se),
    lower = as.numeric(wolf_elev_raw$lower.CL),
    upper = as.numeric(wolf_elev_raw$upper.CL)))

wolf_tri_raw <- modavgPred(
  cand.set = top_mod,
  modnames = names(top_mod),
  newdata = tri_df,
  parm.type = "psi",
  type = "response",
  c.hat = summer_global_chat,
  nobs = n_siteyears,
  uncond.se = "revised",
  species = "wolf")

wolf_tri_pred <- cbind(
  tri_df,
  data.frame(
    Predicted = as.numeric(wolf_tri_raw$mod.avg.pred),
    SE = as.numeric(wolf_tri_raw$uncond.se),
    lower = as.numeric(wolf_tri_raw$lower.CL),
    upper = as.numeric(wolf_tri_raw$upper.CL)))

wolf_og_raw <- modavgPred(
  cand.set = top_mod,
  modnames = names(top_mod),
  newdata = og_df,
  parm.type = "psi",
  type = "response",
  c.hat = summer_global_chat,
  nobs = n_siteyears,
  uncond.se = "revised",
  species = "wolf")

wolf_og_pred <- cbind(
  og_df,
  data.frame(
    Predicted = as.numeric(wolf_og_raw$mod.avg.pred),
    SE = as.numeric(wolf_og_raw$uncond.se),
    lower = as.numeric(wolf_og_raw$lower.CL),
    upper = as.numeric(wolf_og_raw$upper.CL)))

# # Plot wolf marginal occupancy probabilities

# Coast
wolf_coast_pred_plot <- ggplot(wolf_coast_pred, aes(x = coast_decay, y = Predicted)) +
  geom_line(linewidth = 1, color = "lightblue3") +
  geom_ribbon(aes(ymin = pmax(lower, 0), ymax = pmin(upper, 1)), fill = "lightblue3", alpha = 0.45) +
  xlab("Distance from Coast (Decay Rate)") +
  ylab("Wolf Occupancy Probability") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic() +
  theme(axis.title.y = element_blank())

# Elevation
wolf_elev_pred_plot <- ggplot(wolf_elev_pred, aes(x = elevation, y = Predicted)) +
  geom_line(linewidth = 1, color = "thistle3") +
  geom_ribbon(aes(ymin = pmax(lower, 0), ymax = pmin(upper, 1)), fill = "thistle3", alpha = 0.45) +
  xlab("Elevation (m)") +
  ylab("Wolf Occupancy Probability") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic() +
  theme(axis.title.y = element_blank())

# TRI
wolf_tri_pred_plot <- ggplot(wolf_tri_pred, aes(x = mean_tri_500, y = Predicted)) +
  geom_line(linewidth = 1, color = "lightgoldenrod") +
  geom_ribbon(aes(ymin = pmax(lower, 0), ymax = pmin(upper, 1)), fill = "lightgoldenrod", alpha = 0.45) +
  xlab("Mean Terrain Ruggedness Index (500 m buffer)") +
  ylab("Wolf Occupancy Probability") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic() +
  theme(axis.title.y = element_blank())

# Old Growth
wolf_og_pred_plot <- ggplot(wolf_og_pred, aes(x = old_growth_prop_buff, y = Predicted)) +
  geom_line(linewidth = 1, color = "darkolivegreen4") +
  geom_ribbon(aes(ymin = pmax(lower, 0), ymax = pmin(upper, 1)), fill = "darkolivegreen4", alpha = 0.45) +
  xlab("Proportion of Old Growth (500 m buffer)") +
  ylab("Wolf Occupancy Probability") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic() +
  theme(axis.title.y = element_blank())

# Combine plots with patchwork
(wolf_coast_pred_plot + wolf_elev_pred_plot + wolf_tri_pred_plot + wolf_og_pred_plot) & 
  theme(
    # Axis values / tick labels
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    # Axis labels / titles
    axis.title.x = element_text(size = 16, margin = margin(t = 8)))

###################################
# 12.0 Create combined MSOM summary table with c-hat-adjusted SEs/CIs
###################################

# Deer conditional on wolf
deer_wolf_present_raw <- modavgPred(
  cand.set = top_mod,
  modnames = names(top_mod),
  newdata = state_newdata,
  parm.type = "psi",
  type = "response",
  c.hat = summer_global_chat,
  nobs = n_siteyears,
  uncond.se = "revised",
  species = "deer",
  cond = "wolf")

deer_wolf_present_summary <- data.frame(
  Predicted = as.numeric(deer_wolf_present_raw$mod.avg.pred),
  SE = as.numeric(deer_wolf_present_raw$uncond.se),
  lower = as.numeric(deer_wolf_present_raw$lower.CL),
  upper = as.numeric(deer_wolf_present_raw$upper.CL))

deer_wolf_absent_raw <- modavgPred(
  cand.set = top_mod,
  modnames = names(top_mod),
  newdata = state_newdata,
  parm.type = "psi",
  type = "response",
  c.hat = summer_global_chat,
  nobs = n_siteyears,
  uncond.se = "revised",
  species = "deer",
  cond = "-wolf")

deer_wolf_absent_summary <- data.frame(
  Predicted = as.numeric(deer_wolf_absent_raw$mod.avg.pred),
  SE = as.numeric(deer_wolf_absent_raw$uncond.se),
  lower = as.numeric(deer_wolf_absent_raw$lower.CL),
  upper = as.numeric(deer_wolf_absent_raw$upper.CL))

deer_cond_present   <- mean(deer_wolf_present_summary$Predicted, na.rm = TRUE)
deer_cond_present_se <- mean(deer_wolf_present_summary$SE, na.rm = TRUE)
deer_cond_present_l <- mean(pmax(deer_wolf_present_summary$lower, 0), na.rm = TRUE)
deer_cond_present_u <- mean(pmin(deer_wolf_present_summary$upper, 1), na.rm = TRUE)

deer_cond_absent   <- mean(deer_wolf_absent_summary$Predicted, na.rm = TRUE)
deer_cond_absent_se <- mean(deer_wolf_absent_summary$SE, na.rm = TRUE)
deer_cond_absent_l <- mean(pmax(deer_wolf_absent_summary$lower, 0), na.rm = TRUE)
deer_cond_absent_u <- mean(pmin(deer_wolf_absent_summary$upper, 1), na.rm = TRUE)

# Wolf conditional on deer
wolf_deer_present_raw <- modavgPred(
  cand.set = top_mod,
  modnames = names(top_mod),
  newdata = state_newdata,
  parm.type = "psi",
  type = "response",
  c.hat = summer_global_chat,
  nobs = n_siteyears,
  uncond.se = "revised",
  species = "wolf",
  cond = "deer")

wolf_deer_present_summary <- data.frame(
  Predicted = as.numeric(wolf_deer_present_raw$mod.avg.pred),
  SE = as.numeric(wolf_deer_present_raw$uncond.se),
  lower = as.numeric(wolf_deer_present_raw$lower.CL),
  upper = as.numeric(wolf_deer_present_raw$upper.CL))

wolf_deer_absent_raw <- modavgPred(
  cand.set = top_mod,
  modnames = names(top_mod),
  newdata = state_newdata,
  parm.type = "psi",
  type = "response",
  c.hat = summer_global_chat,
  nobs = n_siteyears,
  uncond.se = "revised",
  species = "wolf",
  cond = "-deer")

wolf_deer_absent_summary <- data.frame(
  Predicted = as.numeric(wolf_deer_absent_raw$mod.avg.pred),
  SE = as.numeric(wolf_deer_absent_raw$uncond.se),
  lower = as.numeric(wolf_deer_absent_raw$lower.CL),
  upper = as.numeric(wolf_deer_absent_raw$upper.CL))

wolf_cond_present   <- mean(wolf_deer_present_summary$Predicted, na.rm = TRUE)
wolf_cond_present_se <- mean(wolf_deer_present_summary$SE, na.rm = TRUE)
wolf_cond_present_l <- mean(pmax(wolf_deer_present_summary$lower, 0), na.rm = TRUE)
wolf_cond_present_u <- mean(pmin(wolf_deer_present_summary$upper, 1), na.rm = TRUE)

wolf_cond_absent   <- mean(wolf_deer_absent_summary$Predicted, na.rm = TRUE)
wolf_cond_absent_se <- mean(wolf_deer_absent_summary$SE, na.rm = TRUE)
wolf_cond_absent_l <- mean(pmax(wolf_deer_absent_summary$lower, 0), na.rm = TRUE)
wolf_cond_absent_u <- mean(pmin(wolf_deer_absent_summary$upper, 1), na.rm = TRUE)

# Combined MSOM summary table
msom_summary <- data.frame(
  species = c("Deer", "Wolf"),
  
  marginal_occupancy = c(null_d_mo, null_w_mo),
  mo_lower = c(null_d_lower, null_w_lower),
  mo_upper = c(null_d_upper, null_w_upper),
  
  marginal_detection = c(null_d_det, null_w_det),
  det_lower = c(null_d_det_lower, null_w_det_lower),
  det_upper = c(null_d_det_upper, null_w_det_upper),
  
  conditional_present = c(deer_cond_present, wolf_cond_present),
  cond_present_SE = c(deer_cond_present_se, wolf_cond_present_se),
  cond_present_lower = c(deer_cond_present_l, wolf_cond_present_l),
  cond_present_upper = c(deer_cond_present_u, wolf_cond_present_u),
  
  conditional_absent = c(deer_cond_absent, wolf_cond_absent),
  cond_absent_SE = c(deer_cond_absent_se, wolf_cond_absent_se),
  cond_absent_lower = c(deer_cond_absent_l, wolf_cond_absent_l),
  cond_absent_upper = c(deer_cond_absent_u, wolf_cond_absent_u))

msom_summary


