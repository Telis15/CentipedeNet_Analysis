source('scripts/00_AI_Export_Utils.R')
library(tidyverse)
library(rlang)

# Load pre-wrangled dataframes and remove what I don't need right now
load("data/WrangledData.RData")

#**CAREFUL**! Clear environment except listed objects
rm(list=setdiff(ls(), c("SampleData")))

# Source the global ggplot2 theme and aesthetic scales
source("scripts/0_PlotTheme.R")

# Boxplot Template Function
PlotBySite <- function(Metric, Label) {
  ggplot(SampleData, aes(x = .data[[Metric]], y = Site)) +
    geom_boxplot(fill = SiteColors(), alpha = 0.8) +
    labs(x = Label, y = "Site") +
    plot_theme
}

# 1. Effort by Gear & Site
labeller <- function(labels){
  new_labels <- c("Cast Net" = "Cast Net (Throws)", "Centipede Net" = "Centipede Net (Net Group Hours)", "Seine" = "Seine (Hauls)")
  return(new_labels[labels])
}

p_Effort <- ggplot(SampleData, aes(y = Site, x = Effort, fill = Gear)) +
  geom_boxplot(alpha = 0.65) +  
  scale_fill_manual(values = GearColors()) +
  facet_wrap(. ~ Gear, scales = "free", nrow = 3,  strip.position = "top", labeller = as_labeller(labeller)) +
  theme(panel.grid.major = element_line(colour = "grey90"), panel.grid.minor = element_line(colour = "grey95"), legend.position = "none") 


# Define numerical variables for standard boxplots
NumericVars <- tribble(
  ~Variable, ~Label,
  "Abundance", "Abundance (count)",
  "CPUE", "CPUE",
  "logCPUE", "log-CPUE (ln(x + 1))",
  "Richness", "Species Richness",
  "Shannon", "Shannon Diversity",
  "Simpson", "Simpson Diversity",
  "BPUE", "Biomass Per Unit Effort (g)",
  "DaylightHrs", "Hours of Daylight Sampling",
  "DaylightPercent", "Percent of Daylight Sampling",
  "SecchiDepth", "Secchi Depth (cm)",
  "SecchiDepth_SD", "Secchi Depth Variability (Standard Deviation)",
  "pH", "pH",
  "DO", "Dissolved Oxygen (mg/L)",
  "PSU", "Salinity (PSU)",
  "PSU_SD", "Salinity Variability (Standard Deviation)",
  "Temperature", "Temperature (°C)"
)

# Generate list of numerical plots dynamically
Plots <- map2(NumericVars$Variable, NumericVars$Label, PlotBySite) |>
  set_names(NumericVars$Variable)

Plots$Effort <- p_Effort

# Site-Level Constant Metrics (Barplots)
SiteConstants <- SampleData |>
  select(Site, Occlusion, Occlusion_SD, MudDominant, Steepness) |>
  mutate(across(c(MudDominant, Steepness), as.numeric)) |>
  distinct()

Plots$Occlusion <- ggplot(SiteConstants, aes(y = Site, x = Occlusion, fill = Site)) +
  geom_bar(stat = "identity", alpha = .8) +
  scale_fill_manual(values = SiteColors()) +
  labs(x = "Percent Occlusion", y = "Site") + plot_theme

Plots$Occlusion_SD <- ggplot(SiteConstants, aes(y = Site, x = Occlusion_SD, fill = Site)) +
  geom_bar(stat = "identity", alpha = .8) +
  scale_fill_manual(values = SiteColors()) +
  labs(x = "Occlusion Variability (Standard Deviation)", y = "Site") + plot_theme

Plots$MudDominant <- ggplot(SiteConstants, aes(x = MudDominant, y = Site, fill = Site)) +
  geom_bar(stat = "identity", alpha = .8) +
  scale_fill_manual(values = SiteColors()) +
  scale_x_continuous(breaks = c(0,1), labels = c("Coarser Substrate", "Mud Dominant")) + 
  labs(x = "Dominant Substrate", y = "Site") + plot_theme

Plots$Steepness <- ggplot(SiteConstants, aes(x = Steepness, y = Site, fill = Site)) +
  geom_bar(stat = "identity", alpha = .8) +
  scale_fill_manual(values = SiteColors()) +
  scale_x_continuous(breaks = c(0,1,2), labels = c("Low", "Medium", "High")) + 
  labs(x = "Bank Steepness", y = "Site") + plot_theme

# Save all plots
if(!dir.exists("output/plots/appendix")) dir.create("output/plots/appendix", recursive = TRUE)

iwalk(Plots, ~ggsave(
  filename = paste0("output/plots/appendix/Appendix_", .y, ".eps"), 
  device = cairo_ps,
  plot = .x,
  width = 2.75,
  height = 3.66,
  units = "in"
))

iwalk(Plots, ~ggsave(
  filename = paste0("output/plots/appendix/Appendix_", .y, ".png"), 
  plot = .x,
  width = 2.75,
  height = 3.66,
  units = "in",
  dpi = 600
))


