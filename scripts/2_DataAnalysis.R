##**Confirmed to be functional 3/6/2026. Made minor adjustments to make more sharable on GitHub**
# This script does not modify any pre-existing dataframes from WrangledData.RData and is not necessary to run prior to "3_Model Selection.R"
# It primarily synthesizes, summarizes, and visualizes data created in "1_DataWrangling.R", with the exception of the iNEXT analysis that takes place in this script.

# Load Packages & Pre-Wrangled Data (From 1_DataWrangling.R) --------------
library(iNEXT)
library(vegan)
library(cowplot)
library(eulerr)
library(tidyverse)
# library(devtools)
library(openxlsx2)
library(extrafont)

# Load pre-wrangled dataframes and remove what I don't need right now
load("data/WrangledData.RData")
load("output/models/Rarefaction.RData")
# load("output/models/Rarefaction_New.RData")
rm(FamilyList, Sample_Community_Matrix)

# Remove likely invalid centipede net sample (Tivives visit 2) where nets were never fully submerged
MergedData <- MergedData %>% 
  filter(SampleID != 'Tiv2_Cent')

SampleData <- SampleData %>% 
  filter(SampleID != 'Tiv2_Cent')

# This matrix includes ALL valid samples and is used for calculating 
# the true total diversity (87 species), NOT for iNEXT.
StudyWide_Incidence_Matrix <- MergedData %>% 
  filter(Effort != 0) %>%
  pivot_wider(id_cols = Species, names_from = SampleID, values_from = Count, values_fn = length, values_fill = 0) %>%
  mutate(across(where(is.numeric), ~ replace(., . != 0, 1))) %>%
  filter(!is.na(Species)) %>%
  column_to_rownames(var = "Species")

# Source the global ggplot2 theme and aesthetic scales
source("scripts/0_PlotTheme.R")

# Load Windows fonts from device and set serif as default.
loadfonts(device = "win")
par(family = "serif")

# Something to fix margin issues with cowplot and serif font (overlapping margins in grid plots).
set_null_device(cairo_pdf)


# Summary Stats -----------------------------------------------------------
# 'Tables' is a growing list of dfs that will be exported to a .xlsx file. Each named entry will be a separate, labelled worksheet.
Tables <- list()

# Species List
Tables$`Species List` <- left_join(GearSpeciesGrid, SpeciesList[, c("Species", "FBname")], by = "Species") %>%
  relocate(Species, "Common Name" = FBname)

# Site Summary
Tables$`Site Summary` <- SiteData %>%
  filter(Gear == "Centipede Net") %>%
  select(Site, Visit, DaylightHrs:Steepness, OcclusionGroupAvg) %>%
  summarize(
    OcclusionSTDev = sd(OcclusionGroupAvg, na.rm = TRUE),
    DominantSubstrate = names(which.max(table(DominantSubstrate))),
    across(where(is.numeric) & !Visit, \(x) mean(x, na.rm = TRUE)), .by = Site
  ) %>%
  relocate(OcclusionSTDev, .after = OcclusionGroupAvg)

# Study Summary
Tables$`Study Summary` <- MergedData %>%
  summarise(
    "Total Catch (n)" = sum(!is.na(Species)),
    "Individual Measurements" = sum(!is.na(StandardLength_mm)),
    "Species Richness" = length(unique(na.omit(Species))),
    "Shannon Diversity" = exp(diversity(t(StudyWide_Incidence_Matrix), groups = "none", index = "shannon")),
    "Simpson Diversity" = diversity(t(StudyWide_Incidence_Matrix), groups = "none", index = "invsimpson"),
    "Total Biomass (kg)" = sum(Weight_g, AverageWeight_g, na.rm = T) / 1000,
    "Minimum Weight (g)" = min(Weight_g, AverageWeight_g, na.rm = T),
    "Maximum Weight (g)" = max(Weight_g, AverageWeight_g, na.rm = T),
    "Minimum Standard Length (mm)" = min(StandardLength_mm, na.rm = T),
    "Maximum Standard Length (mm)" = max(StandardLength_mm, na.rm = T)
  )


# Visit duration info:
SampleData %>% 
  select(SampleID, Site, Visit, StartTime, EndTime) %>% 
  group_by(Site, Visit) %>% 
  summarize(Duration = as.numeric(max(EndTime, na.rm = T) - min(StartTime, na.rm = T))) %>% 
  summary(Duration)


# # iNEXT Analysis
# # Species Accumulation Curves & Rarefaction with iNEXT **These take time! (~30 minutes) --------------------
# Gear_Out_Raw <- iNEXT(Incidence_Matrices, q = c(0,1,2), datatype = "incidence_raw", nboot = 1000)
# Gear_Out_Single_Raw <- iNEXT(Incidence_Matrices[1:3],q = 0, datatype = "incidence_raw", nboot = 1000)
# Gear_Out_Combos_Raw <- iNEXT(Incidence_Matrices[4:7],q = 0, datatype = "incidence_raw", nboot = 1000)
# 
# # Creating individual iNEXT objects for each combo. !!These take ~15 minutes each!!
# Cast <- iNEXT(Incidence_Matrices['Cast Net'], q = 0, datatype = "incidence_raw", nboot = 1000, knots = length(Incidence_Matrices$`Cast Net`))
# Cent <- iNEXT(Incidence_Matrices['Centipede Net'], q = 0, datatype = "incidence_raw", nboot = 1000, knots = length(Incidence_Matrices$`Centipede Net`))
# Seine <- iNEXT(Incidence_Matrices["Seine"], q = 0, datatype = "incidence_raw", nboot = 1000, knots = length(Incidence_Matrices$Seine))
# Cast_Cent <- iNEXT(Incidence_Matrices['Cast Net & Centipede Net'], q = 0, datatype = "incidence_raw", nboot = 1000, knots = length(Incidence_Matrices$`Cast Net & Centipede Net`))
# Cast_Seine <- iNEXT(Incidence_Matrices['Cast Net & Seine'], q = 0, datatype = "incidence_raw", nboot = 1000, knots = length(Incidence_Matrices$`Cast Net & Seine`))
# Cent_Seine <- iNEXT(Incidence_Matrices['Centipede Net & Seine'], q = 0, datatype = "incidence_raw", nboot = 1000, knots = length(Incidence_Matrices$`Centipede Net & Seine`))
# All_Gears <- iNEXT(Incidence_Matrices['All Gears'], q = 0, datatype = "incidence_raw", nboot = 1000, knots = length(Incidence_Matrices$`All Gears`))
# 
# #How many samples would be required to attain set sample coverage? (Reported based on nboot = 1000).
# SampleEst <-
#   estimateD(Incidence_Matrices, q = 0, datatype = "incidence_raw", base = "coverage", level = c(0.75, .95), nboot = 1000) %>%
#   rename("Samples (n)" = t) %>%
#   mutate(SC = paste(100 * round(SC, digits = 2), "% Coverage")) %>%
#   select(-Order.q) %>%
#   arrange(SC) %>%
#   add_row(.after = 7)
# 
# save(
#   list = c(
#     "Gear_Out_Raw", "Gear_Out_Single_Raw", "Gear_Out_Combos_Raw",
#     "Cast", "Cent", "Seine",
#     "Cast_Cent", "Cast_Seine", "Cent_Seine", "All_Gears",
#     "SampleEst"
#   ),
#   file = "output/models/Rarefaction.RData"
# )

## Sample Coverage Plots for Solo Gears & Combo Gears ----------------------
# #**Template**
# # This is a mess, but it allows us to customize pretty much everything about the ggiNEXT plots. Most importantly, using distinct linetypes.
# # Get the ggplot build object
# gb <- ggplot_build(ggiNEXT(GEAR_OUT_df_FROM_INEXT, type = 2, color.var = "Assemblage") +
#                      labs(x = "Sampling Units", y = "Sample Coverage") +
#                      scale_shape_manual(values = GearShapes[c("Cast Net", "Centipede Net", "Seine")]) +
#                      scale_color_manual(values = GearColorsAll()[c("Cast Net", "Centipede Net", "Seine")], aesthetics = c("colour", "fill")) +
#                      scale_x_continuous(expand = c(0.01, 0.01)) + # Reduce padding
#                      scale_y_continuous(expand = c(0.01, 0.01), labels = scales::percent) +
#                      # ggiNEXT seems to ignore most ggplot global settings
#                      theme_classic(base_family = "serif", base_size = 9) +
#                      theme(
#                        panel.background = element_rect(fill = "transparent"),
#                        legend.background = element_rect(fill = "transparent"),
#                        legend.title = element_blank(),
#                        legend.position = "inside",
#                        legend.justification.inside = c(1, 0), #(1,1) is upper right
#                        legend.key.spacing.y = unit(0.05, "in") # Vertical space between legend keys (legend.key.spacing should work, but doesn't)
#                      ) +
#                      # Remove all legend entries except for shape & linetype (within color for some reason); format appearance
#                      guides(
#                        linetype = guide_none(),
#                        fill = guide_none(),
#                        shape = guide_none(),
#                        color = guide_legend(override.aes = list(
#                          shape = GearShapes[c("Cast Net", "Centipede Net", "Seine")],
#                          size = 2,
#                          linetype = c(GearLineTypes$`Cast Net`, GearLineTypes$`Centipede Net`, GearLineTypes$Seine),
#                          fill = NA,
#                          linewidth = 0.5
#                        ))
#                      ))
# 
# gb$data[[1]] <- gb$data[[1]] %>% # Forcibly reduce size of shapes on plot
#   mutate(size = 3)
# 
# gb$data[[2]] <- gb$data[[2]] %>% # Forcibly reduce linewidth and assign linetypes by group (should be twice as many groups as assemblages. Low numbers are rarefied, high numbers are extrapolated). Can also assign as GearLineTypes[[1,2,etc.]]
#   mutate(
#     linewidth = 0.75,
#     linetype = case_when(
#       group < 4 ~ linetype,
#       group == 4 ~ GearLineTypes$`Cast Net`,
#       group == 5 ~ GearLineTypes$`Centipede Net`,
#       group == 6 ~ GearLineTypes$Seine
#     )
#   )
# # Control spacing between legend entries (keys) and key dimensions
# gb$plot$theme$legend.key.spacing.y <- unit(0.05, "in")
# gb$plot$theme$legend.key.height <- unit(0.15, "in")
# gb$plot$theme$legend.key.width <- unit(0.5, "in")
# 
# # Create the gtable
# Plot_Name <- ggplot_gtable(gb)
# #
# # # Display the plot
# plot(Plot_Name)
# If it's a grid you're lookin' for:
# Grid_Name <- plot_grid(Plot_Name1, Plot_Name2, nrow = 2, labels = Letters(), label_x = 0.11, align = "hv", label_fontface = "plain", label_size = 12, label_fontfamily = "serif")

# ggsave(plot = Plot_Name, "output/plots/Plot_Name.png", width = 6, height = 4, units = "in")


# Sample Coverage Plots
# --- Streamlining: Create a common theme for both plots ---
common_theme <- function() {
  theme_classic(base_family = "serif", base_size = 9) +
    theme(
      panel.background = element_rect(fill = "transparent"),
      legend.background = element_rect(fill = "transparent"),
      legend.title = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(1, 0),
      legend.justification = c(1, 0),
      legend.key.spacing.y = unit(0.05, "in")
    )
}

# --- Singles_Coverage_Plot (Refactored with Final Legend Fix) ---

# 1. Create the initial ggplot object | Sample Completeness Curve
p_singles <- ggiNEXT(Gear_Out_Single_Raw, type = 2, color.var = "Assemblage", se = TRUE) +
  labs(x = "Sampling Units", y = "Sample Coverage") +
  scale_shape_manual(values = GearShapes[c("Cast Net", "Centipede Net", "Seine")]) +
  scale_color_manual(values = GearColorsAll()[c("Cast Net", "Centipede Net", "Seine")], aesthetics = c("colour", "fill")) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01), labels = scales::percent) +
  common_theme() +
  # This scale overrides the default linetype legend from ggiNEXT and removes it
  scale_linetype_manual(values = c("rarefaction" = "solid", "extrapolation" = "dashed"), guide = "none") +
  guides(
    fill = "none",
    shape = "none",
    color = guide_legend(override.aes = list(
      shape = GearShapes[c("Cast Net", "Centipede Net", "Seine")],
      size = 2,
      stroke = 1,
      linewidth = 0.4,
      alpha = 0.7,
      linetype = unlist(GearLineTypes[c("Cast Net", "Centipede Net", "Seine")])
    ))
  )

# 2. Apply the fix to remove ribbon from legend
for (i in seq_along(p_singles$layers)) {
  if (inherits(p_singles$layers[[i]]$geom, "GeomRibbon")) {
    p_singles$layers[[i]]$show.legend <- FALSE
  }
}

# 3. Build the modified plot object
gb_singles <- ggplot_build(p_singles)

gb_singles$data[[1]] <- gb_singles$data[[1]] %>% mutate(size = 3)
gb_singles$data[[2]] <- gb_singles$data[[2]] %>%
  mutate(
    linewidth = 0.75,
    linetype = case_when(
      group >= 4 ~ "15",
      group == 1 ~ GearLineTypes$`Cast Net`,
      group == 2 ~ GearLineTypes$`Centipede Net`,
      group == 3 ~ GearLineTypes$Seine
    )
  )

gb_singles$plot$theme$legend.key.height <- unit(0.15, "in")
gb_singles$plot$theme$legend.key.width <- unit(0.35, "in")

Singles_Coverage_Plot <- ggplot_gtable(gb_singles)


# --- Combos_Coverage_Plot (Refactored with Final Legend Fix) ---

# 1. Create the initial ggplot object | Sample Completeness Curve
p_combos <- ggiNEXT(Gear_Out_Combos_Raw, type = 2, color.var = "Assemblage", se = TRUE) +
  labs(x = "Sampling Units", y = NULL) +
  scale_shape_manual(values = GearShapes[c("All Gears", "Cast Net & Centipede Net", "Cast Net & Seine", "Centipede Net & Seine")]) +
  scale_color_manual(values = GearColorsAll()[c("All Gears", "Cast Net & Centipede Net", "Cast Net & Seine", "Centipede Net & Seine")], aesthetics = c("colour", "fill")) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01), labels = scales::percent, limits = c(.20, NA)) +
  common_theme() +
  # This scale overrides the default linetype legend from ggiNEXT and removes it
  scale_linetype_manual(values = c("rarefaction" = "solid", "extrapolation" = "dashed"), guide = "none") +
  guides(
    fill = "none",
    shape = "none",
    color = guide_legend(override.aes = list(
      shape = GearShapes[c("All Gears", "Cast Net & Centipede Net", "Cast Net & Seine", "Centipede Net & Seine")],
      size = 2,
      stroke = 1,
      linewidth = 0.4,
      alpha = 0.7,
      linetype = unlist(GearLineTypes[c("All Gears", "Cast Net & Centipede Net", "Cast Net & Seine", "Centipede Net & Seine")])
    ))
  )

# 2. Apply the fix to remove ribbon from legend
for (i in seq_along(p_combos$layers)) {
  if (inherits(p_combos$layers[[i]]$geom, "GeomRibbon")) {
    p_combos$layers[[i]]$show.legend <- FALSE
  }
}

# 3. Build the modified plot object
gb_combos <- ggplot_build(p_combos)

gb_combos$data[[1]] <- gb_combos$data[[1]] %>% mutate(size = 3)
gb_combos$data[[2]] <- gb_combos$data[[2]] %>%
  mutate(
    linewidth = 0.75,
    linetype = case_when(
      group >= 5 ~ "15",
      group == 1 ~ GearLineTypes$`All Gears`,
      group == 2 ~ GearLineTypes$`Cast Net & Centipede Net`,
      group == 3 ~ GearLineTypes$`Cast Net & Seine`,
      group == 4 ~ GearLineTypes$`Centipede Net & Seine`
    )
  )

gb_combos$plot$theme$legend.key.height <- unit(0.15, "in")
gb_combos$plot$theme$legend.key.width <- unit(0.35, "in")

Combos_Coverage_Plot <- ggplot_gtable(gb_combos)


# --- Combine Plots ---
Coverage_Grid <-
  plot_grid(Singles_Coverage_Plot, Combos_Coverage_Plot, nrow = 1, labels = c("Individual Gears", "Gear Combinations"), label_x = 0.05, align = "hv", rel_widths = c(1.03, 1), label_fontface = "plain", label_size = 9, label_fontfamily = "serif")

# Display the final plot
plot(Coverage_Grid)

# --- Save the Final Plot for Publication ---
ggsave(plot = Coverage_Grid, 
       filename = "output/plots/Coverage_2Panel_Final_Fig6.eps", 
       device = cairo_ps, 
       width = 5.62, 
       height = 3.75, 
       units = "in")
ggsave(plot = Coverage_Grid, 
       filename = "output/plots/Coverage_2Panel_Final_Fig6.png", 
       width = 5.62, 
       height = 3.75, 
       units = "in",
       dpi = 600)


# All q Plot | Diversity Accumulation Curve --------------------------------------------------------------
# 1. Create the initial ggplot object
p <- ggiNEXT(Gear_Out_Raw, type = 1, se = TRUE, facet.var = "Order.q", color.var = "Assemblage") +
  labs(x = "Sampling Units", y = "Effective Species") +
  scale_shape_manual(values = GearShapes) +
  scale_color_manual(values = GearColorsAll(), aesthetics = c("colour", "fill")) +
  scale_x_continuous(expand = c(0.01, 0.01)) + # Remove padding
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme_classic(base_family = "serif", base_size = 9) +
  theme(
    panel.background = element_rect(fill = "transparent"),
    legend.background = element_rect(fill = "transparent"), 
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.justification = "center",
    strip.text = element_text(vjust = 0.5) # Vertically center facet labels
  ) +
  facet_wrap(vars(Order.q), scales="free_y", labeller = labeller(Order.q = c("0" = "Species Richness", "1" = "Shannon Diversity", "2" = "Simpson Diversity"))) +
  guides(
    linetype = "none",
    fill = "none",
    shape = "none",
    color = guide_legend(override.aes = list(
      shape = GearShapes,
      size = 2,
      stroke = 1,
      linewidth = 0.4,
      alpha = 1,
      linetype = unlist(GearLineTypes) # Add linetypes back to the legend
    ))
  )

# 2. Apply the definitive fix to remove ribbon from legend
for (i in seq_along(p$layers)) {
  if (inherits(p$layers[[i]]$geom, "GeomRibbon")) {
    p$layers[[i]]$show.legend <- FALSE
  }
}

# 3. Build the modified plot object
gb <- ggplot_build(p)

gb$data[[1]] <- gb$data[[1]] %>%
  mutate(size = 3)

gb$data[[2]] <- gb$data[[2]] %>%
  mutate(
    linewidth = 0.75,
    linetype = case_when(
      group >= 8 ~ "15",
      group == 1 ~ GearLineTypes[[1]],
      group == 2 ~ GearLineTypes[[2]],
      group == 3 ~ GearLineTypes[[3]],
      group == 4 ~ GearLineTypes[[4]],
      group == 5 ~ GearLineTypes[[5]],
      group == 6 ~ GearLineTypes[[6]],
      group == 7 ~ GearLineTypes[[7]]
    )
  )

gb$plot$theme$legend.key.height <- unit(0.15, "in")
gb$plot$theme$legend.key.width <- unit(0.5, "in")

# Create the gtable
All_q_Plot <- ggplot_gtable(gb)

# Display the plot
plot(All_q_Plot)

# --- Save the Final Plot for Publication --- Figure 9
ggsave(plot = All_q_Plot, 
       filename = "output/plots/All_q_Plot_Final.eps", 
       device = cairo_ps,
       width = 5.62, 
       height = 3.75, 
       units = "in")
ggsave(plot = All_q_Plot, 
       filename = "output/plots/All_q_Plot_Final.png", 
       width = 5.62, 
       height = 3.75, 
       units = "in",
       dpi = 600)

# ## Testing Type 3 Plot:
# # All q Plot | Diversity Accumulation Curve --------------------------------------------------------------
# # 1. Create the initial ggplot object
# p <- ggiNEXT(Gear_Out_Raw, type = 3, se = TRUE, facet.var = "Order.q", color.var = "Assemblage") +
#   labs(x = "Percent Coverage", y = "Effective Species") +
#   scale_shape_manual(values = GearShapes) +
#   scale_color_manual(values = GearColorsAll(), aesthetics = c("colour", "fill")) +
#   scale_x_continuous(expand = c(0.01, 0.01), labels = scales::percent) + # Remove padding
#   scale_y_continuous(expand = c(0.01, 0.01)) +
#   theme_classic(base_family = "serif", base_size = 9) +
#   theme(
#     panel.background = element_rect(fill = "transparent"),
#     legend.background = element_rect(fill = "transparent"), 
#     legend.title = element_blank(),
#     legend.position = "bottom",
#     legend.justification = "center",
#     strip.text = element_text(vjust = 0.5) # Vertically center facet labels
#   ) +
#   facet_wrap(vars(Order.q), scales="free_y", labeller = labeller(Order.q = c("0" = "Species Richness", "1" = "Shannon Diversity", "2" = "Simpson Diversity"))) +
#   guides(
#     linetype = "none",
#     fill = "none",
#     shape = "none",
#     color = guide_legend(override.aes = list(
#       shape = GearShapes,
#       size = 2,
#       stroke = 1,
#       linewidth = 0.4,
#       alpha = 1,
#       linetype = unlist(GearLineTypes) # Add linetypes back to the legend
#     ))
#   )
# 
# # 2. Apply the definitive fix to remove ribbon from legend
# for (i in seq_along(p$layers)) {
#   if (inherits(p$layers[[i]]$geom, "GeomRibbon")) {
#     p$layers[[i]]$show.legend <- FALSE
#   }
# }
# 
# # 3. Build the modified plot object
# gb <- ggplot_build(p)
# 
# gb$data[[1]] <- gb$data[[1]] %>%
#   mutate(size = 3)
# 
# gb$data[[2]] <- gb$data[[2]] %>%
#   mutate(
#     linewidth = 0.75,
#     linetype = case_when(
#       group >= 8 ~ "15",
#       group == 1 ~ GearLineTypes[[1]],
#       group == 2 ~ GearLineTypes[[2]],
#       group == 3 ~ GearLineTypes[[3]],
#       group == 4 ~ GearLineTypes[[4]],
#       group == 5 ~ GearLineTypes[[5]],
#       group == 6 ~ GearLineTypes[[6]],
#       group == 7 ~ GearLineTypes[[7]]
#     )
#   )
# 
# gb$plot$theme$legend.key.height <- unit(0.15, "in")
# gb$plot$theme$legend.key.width <- unit(0.5, "in")
# 
# # Create the gtable
# Type3_Curve <- ggplot_gtable(gb)
# 
# # Display the plot
# plot(Type3_Curve)
# 
# # --- Save the Final Plot for Publication --- Figure 9
# ggsave(plot = All_q_Plot, 
#        filename = "output/plots/Type3_Curve.eps", 
#        device = cairo_ps,
#        width = 6, 
#        height = 4, 
#        units = "in")
# ggsave(plot = All_q_Plot, 
#        filename = "output/plots/Type3_Curve.png", 
#        width = 6, 
#        height = 4, 
#        units = "in")

# Skip to here if not re-running iNEXT analysis ---------------------------
# Table 5 How many samples would be required to attain set sample coverage? (Reported based on nboot = 1000).
Tables$`iNEXT Coverage Projection` <- SampleEst

SampleEst_Wide <- SampleEst %>%
  filter(Assemblage %in% c("All Gears", "Cast Net & Seine")) %>%
  select(SC, Assemblage, "Samples (n)") %>%
  pivot_wider(names_from = Assemblage, values_from = "Samples (n)") %>%
  mutate(`% Improvement` = (`Cast Net & Seine` - `All Gears`) / `Cast Net & Seine` * 100)

Tables$`iNEXT Coverage` <-
  Gear_Out_Raw$DataInfo %>%
  select(1:2, 4:5) %>%
  rename("Number of Samples" = T, "Observed Richness" = S.obs, "Estimated Sample Coverage" = SC)

# Table 4 Hill Number Estimates:
Tables$`iNEXT Estimates` <-
  Gear_Out_Raw$AsyEst %>%
  rename(
    "Metric" = "Diversity",
    "Observed" = "Observed",
    "Estimated" = "Estimator"
  ) %>%
  arrange(match(Metric, c("Species richness", "Shannon diversity", "Simpson diversity")), match(Assemblage, c("Cast Net", "Centipede Net", "Seine", "Cast Net & Centipede Net", "Cast Net & Seine", "Centipede Net & Seine", "All Gears"))) %>%
  add_row(.after = 14) %>%
  add_row(.after = 7)

# Gear Combo Diversity
Gear_Diversity <- Gear_Out_Raw$AsyEst %>%
      as.data.frame() %>%
      filter(Diversity == "Species richness") %>%
      select(Assemblage, Observed) %>%
      rename(Gear = Assemblage, Richness = Observed) %>%
  left_join(
    Gear_Out_Raw$AsyEst %>%
      as.data.frame() %>%
      filter(Diversity == "Shannon diversity") %>%
      select(Assemblage, Observed) %>%
      rename(Gear = Assemblage, Shannon = Observed),
    by = "Gear"
  ) %>%
  left_join(
    Gear_Out_Raw$AsyEst %>%
      as.data.frame() %>%
      filter(Diversity == "Simpson diversity") %>%
      select(Assemblage, Observed) %>%
      rename(Gear = Assemblage, Simpson = Observed),
    by = "Gear"
  )

# Summary data by gear

# --- 1. Create a table of TOTAL stats (Catch, Richness, Diversity) ---
# This table calculates the TRUE study-wide totals for all 7 combinations
# by pooling all valid samples, regardless of site-visit matching.

# --- First, create all the correct (un-paired) incidence matrices for this table ---
Matrix_Cast <- MergedData %>% filter(Gear == "Cast Net" & Effort != 0) %>% pivot_wider(id_cols = Species, names_from = SampleID, values_from = Count, values_fn = length, values_fill = 0) %>% mutate(across(where(is.numeric), ~ replace(., . != 0, 1))) %>% filter(!is.na(Species)) %>% column_to_rownames(var = "Species")
Matrix_Cent <- MergedData %>% filter(Gear == "Centipede Net" & Effort != 0 & SampleID != 'Tiv2_Cent') %>% pivot_wider(id_cols = Species, names_from = SampleID, values_from = Count, values_fn = length, values_fill = 0) %>% mutate(across(where(is.numeric), ~ replace(., . != 0, 1))) %>% filter(!is.na(Species)) %>% column_to_rownames(var = "Species")
Matrix_Seine <- MergedData %>% filter(Gear == "Seine" & Effort != 0) %>% pivot_wider(id_cols = Species, names_from = SampleID, values_from = Count, values_fn = length, values_fill = 0) %>% mutate(across(where(is.numeric), ~ replace(., . != 0, 1))) %>% filter(!is.na(Species)) %>% column_to_rownames(var = "Species")
Matrix_Cast_Cent <- MergedData %>% filter((Gear == "Cast Net" | Gear == "Centipede Net") & Effort != 0 & SampleID != 'Tiv2_Cent') %>% pivot_wider(id_cols = Species, names_from = SampleID, values_from = Count, values_fn = length, values_fill = 0) %>% mutate(across(where(is.numeric), ~ replace(., . != 0, 1))) %>% filter(!is.na(Species)) %>% column_to_rownames(var = "Species")
Matrix_Cast_Seine <- MergedData %>% filter((Gear == "Cast Net" | Gear == "Seine") & Effort != 0) %>% pivot_wider(id_cols = Species, names_from = SampleID, values_from = Count, values_fn = length, values_fill = 0) %>% mutate(across(where(is.numeric), ~ replace(., . != 0, 1))) %>% filter(!is.na(Species)) %>% column_to_rownames(var = "Species")
Matrix_Cent_Seine <- MergedData %>% filter((Gear == "Centipede Net" | Gear == "Seine") & Effort != 0 & SampleID != 'Tiv2_Cent') %>% pivot_wider(id_cols = Species, names_from = SampleID, values_from = Count, values_fn = length, values_fill = 0) %>% mutate(across(where(is.numeric), ~ replace(., . != 0, 1))) %>% filter(!is.na(Species)) %>% column_to_rownames(var = "Species")
# Matrix_All_Gears is your existing 'StudyWide_Incidence_Matrix' (which is correct)

# --- Now, build the 7-row table of total stats ---
total_stats_final <- bind_rows(
  # Cast Net
  GearSpeciesGrid %>% select(`Cast Net`) %>% summarize("Total Catch (n)" = sum(`Cast Net`, na.rm = TRUE)) %>%
    mutate(Gear = "Cast Net", "Species Richness" = nrow(Matrix_Cast),
           "Shannon Diversity" = exp(diversity(t(Matrix_Cast), groups = "none", index = "shannon")),
           "Simpson Diversity" = diversity(t(Matrix_Cast), groups = "none", index = "invsimpson")),
  # Centipede Net
  GearSpeciesGrid %>% select(`Centipede Net`) %>% summarize("Total Catch (n)" = sum(`Centipede Net`, na.rm = TRUE)) %>%
    mutate(Gear = "Centipede Net", "Species Richness" = nrow(Matrix_Cent),
           "Shannon Diversity" = exp(diversity(t(Matrix_Cent), groups = "none", index = "shannon")),
           "Simpson Diversity" = diversity(t(Matrix_Cent), groups = "none", index = "invsimpson")),
  # Seine
  GearSpeciesGrid %>% select(Seine) %>% summarize("Total Catch (n)" = sum(Seine, na.rm = TRUE)) %>%
    mutate(Gear = "Seine", "Species Richness" = nrow(Matrix_Seine),
           "Shannon Diversity" = exp(diversity(t(Matrix_Seine), groups = "none", index = "shannon")),
           "Simpson Diversity" = diversity(t(Matrix_Seine), groups = "none", index = "invsimpson")),
  # Cast Net & Centipede Net
  GearSpeciesGrid %>% select(`Cast Net`, `Centipede Net`) %>% summarize("Total Catch (n)" = sum(., na.rm = TRUE)) %>%
    mutate(Gear = "Cast Net & Centipede Net", "Species Richness" = nrow(Matrix_Cast_Cent),
           "Shannon Diversity" = exp(diversity(t(Matrix_Cast_Cent), groups = "none", index = "shannon")),
           "Simpson Diversity" = diversity(t(Matrix_Cast_Cent), groups = "none", index = "invsimpson")),
  # Cast Net & Seine
  GearSpeciesGrid %>% select(`Cast Net`, Seine) %>% summarize("Total Catch (n)" = sum(., na.rm = TRUE)) %>%
    mutate(Gear = "Cast Net & Seine", "Species Richness" = nrow(Matrix_Cast_Seine),
           "Shannon Diversity" = exp(diversity(t(Matrix_Cast_Seine), groups = "none", index = "shannon")),
           "Simpson Diversity" = diversity(t(Matrix_Cast_Seine), groups = "none", index = "invsimpson")),
  # Centipede Net & Seine
  GearSpeciesGrid %>% select(`Centipede Net`, Seine) %>% summarize("Total Catch (n)" = sum(., na.rm = TRUE)) %>%
    mutate(Gear = "Centipede Net & Seine", "Species Richness" = nrow(Matrix_Cent_Seine),
           "Shannon Diversity" = exp(diversity(t(Matrix_Cent_Seine), groups = "none", index = "shannon")),
           "Simpson Diversity" = diversity(t(Matrix_Cent_Seine), groups = "none", index = "invsimpson")),
  # All Gears (using the true study-wide objects)
  tibble(
    Gear = "All Gears",
    "Total Catch (n)" = sum(GearSpeciesGrid[,c("Cast Net", "Centipede Net", "Seine")], na.rm = TRUE),
    "Species Richness" = length(unique(na.omit(MergedData$Species))),
    "Shannon Diversity" = exp(diversity(t(StudyWide_Incidence_Matrix), groups = "none", index = "shannon")),
    "Simpson Diversity" = diversity(t(StudyWide_Incidence_Matrix), groups = "none", index = "invsimpson")
  )
)

# --- 2. Create the final table by merging MEAN and TOTAL stats ---
Tables$`Gear Study Summary` <- full_join(
  bind_rows(
    SampleData %>%
      summarize(
        "Total Samples" = length(Richness), "Mean Effort" = mean(Effort),
        "Mean Catch (n)" = mean(Abundance, na.rm = T), "Mean CPUE" = mean(CPUE, na.rm = T),
        "Mean Biomass (g)" = mean(Biomass, na.rm = T), "Mean Richness" = mean(Richness, na.rm = T),
        "Mean Shannon Diversity" = mean(Shannon, na.rm = T), "Mean Simpson Diversity" = mean(Simpson, na.rm = T),
        .by = Gear
      ),
    SampleData %>%
      summarize(
        Gear = "All Gears", "Total Samples" = length(Richness), "Mean Effort" = NA_real_,
        "Mean Catch (n)" = mean(Abundance, na.rm = T), "Mean CPUE" = NA_real_,
        "Mean Biomass (g)" = mean(Biomass, na.rm = T), "Mean Richness" = mean(Richness, na.rm = T),
        "Mean Shannon Diversity" = mean(Shannon, na.rm = T), "Mean Simpson Diversity" = mean(Simpson, na.rm = T)
      )
  ),
  total_stats_final,
  by = "Gear"
) %>%
  mutate("Effort Unit" = c("Throws", "Net Group Hours", "Hauls", NA, NA, NA, NA), .after = "Mean Effort") %>% 
  arrange(match(Gear, c("Cast Net", "Centipede Net", "Seine", "Cast Net & Centipede Net", "Cast Net & Seine", "Centipede Net & Seine", "All Gears"))) %>%
  relocate(Gear, `Total Samples`, `Total Catch (n)`, `Species Richness`, `Shannon Diversity`, `Simpson Diversity`, `Mean Catch (n)`, `Mean Richness`, `Mean Shannon Diversity`, `Mean Simpson Diversity`, `Mean Biomass (g)`, `Mean CPUE`)


# Boxplot of lengths by gear
ggplot(MergedData %>% filter(!is.na(StandardLength_mm)), aes(x = Gear, y = StandardLength_mm, fill = Gear)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(outlier.shape = 16, outlier.size = 2, alpha = 0.75) +
  geom_text(
    data = MergedData %>% group_by(Gear) %>% summarize(StandardLength_mm = mean(StandardLength_mm, na.rm = TRUE)),
    aes(label = paste("Mean =", round(StandardLength_mm, 1)), y = StandardLength_mm, family = "serif"),
    vjust = c(1.2, 0.1, 0.2), color = "black"
  ) +
  scale_fill_manual(values = GearColors()[2:4]) +
  labs(y = "Standard Length (mm)") +
  theme_classic(base_size = 9, base_family = "serif") +
  theme(axis.title.x = element_blank(), legend.position = "none")
# ggsave("output/plots/Length_Gear_Boxplot.png", width = 6, height = 6)


# Boxplot of weights by gear **Includes 'AverageWeight_g' for batch-measured catch for samples with >30 of one species. That only occurred in Seines.
ggplot(MergedData %>% mutate(Weight_g = if_else(is.na(Weight_g), AverageWeight_g, Weight_g)) %>% select(Gear, Weight_g) %>% na.omit(), aes(x = Gear, y = Weight_g, fill = Gear)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(outlier.shape = 16, outlier.size = 2) +
  geom_text(
    data = MergedData %>% mutate(Weight_g = if_else(is.na(Weight_g), AverageWeight_g, Weight_g)) %>% select(Gear, Weight_g) %>% na.omit() %>% group_by(Gear) %>% summarize(Weight_g = mean(Weight_g, na.rm = TRUE)),
    aes(label = paste("Mean =", round(Weight_g, 1)), y = Weight_g, family = "serif"),
    vjust = c(3.2, 3.5, 2), color = "black"
  ) +
  scale_fill_manual(values = GearColors()[2:4]) +
  labs(y = "Weight (g)") +
  theme_classic(base_size = 9, base_family = "serif") +
  theme(axis.title.x = element_blank(), legend.position = "none")
# ggsave("output/plots/Weight_Gear_Boxplot.png", width = 6, height = 6)


# Create two-panel figure of violin plot & scatterplot (length x weight by gear)
# --- 1. Streamline by preparing data and summary stats first ---

# Create a clean data frame for plotting, removing NA values
plot_data <- MergedData %>% 
  filter(!is.na(StandardLength_mm))

# Create a summary data frame for labels and stats
summary_data <- plot_data %>%
  group_by(Gear) %>%
  summarize(
    n = n(),
    Mean = round(mean(StandardLength_mm, na.rm = TRUE), 0)
  ) %>%
  mutate(
    # Create the x-axis labels with sample sizes
    x_label = paste0(Gear, "\n(n = ", n, ")"),
    # Create the mean labels for geom_text
    mean_label = paste("Mean =", Mean)
  )

# --- 2. Build the Violin Plot using the prepared data ---

Violin_Plot_SL <- ggplot(data = plot_data, aes(x = Gear, y = StandardLength_mm, fill = Gear)) +
  # Use geom_jitter for outliers to prevent overplotting, using the same filtered data
  geom_jitter(shape = 16, alpha = 0.75, width = 0.1, data = . %>% filter(StandardLength_mm > 250)) +
  geom_violin(alpha = 0.75) +
  labs(y = "Standard Length (mm)", x = NULL) +
  # Use the pre-calculated labels from the summary data frame
  scale_x_discrete(labels = setNames(summary_data$x_label, summary_data$Gear)) +
  scale_fill_manual(values = GearColors()) +
  # Set a lower y-limit to make space for the mean text
  scale_y_continuous(limits = c(-25, NA)) +
  # Add the mean labels using the summary data frame
  geom_text(
    data = summary_data,
    aes(label = mean_label),
    color = "black", 
    size = 3, # Use a standard relative size
    y = -30     # Position the text below the plot
  ) +
  theme_classic(base_size = 9, base_family = "serif") +
  theme(
    legend.position = "none", 
    axis.text = element_text(size = 9, color = "black")
  )

# Display the plot
plot(Violin_Plot_SL)

# --- 3. Save the plot in the correct format for publication ---
ggsave(plot = Violin_Plot_SL, 
       filename = "output/plots/Length_Gear_Violin_SL.eps", 
       device = cairo_ps,
       width = 2.75, 
       height = 3.66, 
       units = "in")
ggsave(plot = Violin_Plot_SL, 
       filename = "output/plots/Length_Gear_Violin_SL.png", 
       width = 2.75, 
       height = 3.66, 
       units = "in",
       dpi = 600)

# # Scatterplot length x weight by gear
#  # ScatterPlot_LW <-
#    ggplot(MergedData %>% filter(!is.na(Weight_g)), aes(x = StandardLength_mm, y = Weight_g, color = Gear)) +
#     scale_color_manual(values = GearColors()) +
#    scale_y_continuous(limits = c(-20, NA)) +
#     geom_point(size = 2, alpha = 0.75) +
#     labs(x = "Standard Length (mm)", y = "Weight (g)") +
#     theme_classic(base_size = 9, base_family = "serif") +
#     theme(legend.position = "none") #Customize axis text to match violin plot,
#   # nrow = 1, labels = Letters(), label_x = 0, align = "hv", label_fontface = "plain", label_size = 12, label_fontfamily = "serif")
# 
# # ggsave("output/plots/LengthWeight_2Panel.png", width = 6, height = 4)
# # ggsave("output/plots/LengthWeight_Scatter.png", width = 6, height = 4)

# Catch by Site
Tables$`Site Catch` <- MergedData %>%
  filter(!is.na(Species)) %>%
  group_by(Site, Gear) %>%
  summarise(Count = sum(Count), .groups = "keep") %>%
  ungroup() %>%
  pivot_wider(names_from = Gear, values_from = Count) %>%
  mutate("Total Catch" = rowSums(select(., -Site)))

# CPUE by Site
Tables$`Site CPUE` <- SampleData %>%
  group_by(Site, Gear) %>%
  summarise(CPUE = mean(CPUE, na.rm = TRUE), .groups = 'keep') %>%
  pivot_wider(names_from = Gear, values_from = CPUE) %>%
  mutate("All Gears" = rowMeans(across(where(is.numeric)), na.rm = TRUE)) %>%
  bind_rows(
    SampleData %>%
      group_by(Gear) %>%
      summarise(CPUE = mean(CPUE, na.rm = TRUE), .groups = 'drop') %>%
      pivot_wider(names_from = Gear, values_from = CPUE) %>%
      mutate(Site = "All Sites", "All Gears" = rowMeans(across(where(is.numeric)), na.rm = TRUE))
  )


# Summary data by gear

# --- 1. Create a table of TOTAL stats (Catch, Richness, Diversity) ---
# This table calculates the TRUE study-wide totals for all 7 combinations
# by pooling all valid samples, regardless of site-visit matching.

# --- First, create all the correct (un-paired) incidence matrices for this table ---
Matrix_Cast <- dcast(MergedData %>% filter(Gear == "Cast Net" & Effort != 0), Species ~ SampleID, value.var = "Count", fun.aggregate = length, fill = 0) %>% mutate(across(where(is.numeric), ~ replace(., . != 0, 1))) %>% filter(!is.na(Species)) %>% column_to_rownames(var = "Species")
Matrix_Cent <- dcast(MergedData %>% filter(Gear == "Centipede Net" & Effort != 0 & SampleID != 'Tiv2_Cent'), Species ~ SampleID, value.var = "Count", fun.aggregate = length, fill = 0) %>% mutate(across(where(is.numeric), ~ replace(., . != 0, 1))) %>% filter(!is.na(Species)) %>% column_to_rownames(var = "Species")
Matrix_Seine <- dcast(MergedData %>% filter(Gear == "Seine" & Effort != 0), Species ~ SampleID, value.var = "Count", fun.aggregate = length, fill = 0) %>% mutate(across(where(is.numeric), ~ replace(., . != 0, 1))) %>% filter(!is.na(Species)) %>% column_to_rownames(var = "Species")
Matrix_Cast_Cent <- dcast(MergedData %>% filter((Gear == "Cast Net" | Gear == "Centipede Net") & Effort != 0 & SampleID != 'Tiv2_Cent'), Species ~ SampleID, value.var = "Count", fun.aggregate = length, fill = 0) %>% mutate(across(where(is.numeric), ~ replace(., . != 0, 1))) %>% filter(!is.na(Species)) %>% column_to_rownames(var = "Species")
Matrix_Cast_Seine <- dcast(MergedData %>% filter((Gear == "Cast Net" | Gear == "Seine") & Effort != 0), Species ~ SampleID, value.var = "Count", fun.aggregate = length, fill = 0) %>% mutate(across(where(is.numeric), ~ replace(., . != 0, 1))) %>% filter(!is.na(Species)) %>% column_to_rownames(var = "Species")
Matrix_Cent_Seine <- dcast(MergedData %>% filter((Gear == "Centipede Net" | Gear == "Seine") & Effort != 0 & SampleID != 'Tiv2_Cent'), Species ~ SampleID, value.var = "Count", fun.aggregate = length, fill = 0) %>% mutate(across(where(is.numeric), ~ replace(., . != 0, 1))) %>% filter(!is.na(Species)) %>% column_to_rownames(var = "Species")
# Matrix_All_Gears is your existing 'StudyWide_Incidence_Matrix' (which is correct)

# --- Now, build the 7-row table of total stats ---
total_stats_final <- bind_rows(
  # Cast Net
  GearSpeciesGrid %>% select(`Cast Net`) %>% summarize("Total Catch (n)" = sum(`Cast Net`, na.rm = TRUE)) %>%
    mutate(Gear = "Cast Net", "Species Richness" = nrow(Matrix_Cast),
           "Shannon Diversity" = exp(diversity(t(Matrix_Cast), groups = "none", index = "shannon")),
           "Simpson Diversity" = diversity(t(Matrix_Cast), groups = "none", index = "invsimpson")),
  # Centipede Net
  GearSpeciesGrid %>% select(`Centipede Net`) %>% summarize("Total Catch (n)" = sum(`Centipede Net`, na.rm = TRUE)) %>%
    mutate(Gear = "Centipede Net", "Species Richness" = nrow(Matrix_Cent),
           "Shannon Diversity" = exp(diversity(t(Matrix_Cent), groups = "none", index = "shannon")),
           "Simpson Diversity" = diversity(t(Matrix_Cent), groups = "none", index = "invsimpson")),
  # Seine
  GearSpeciesGrid %>% select(Seine) %>% summarize("Total Catch (n)" = sum(Seine, na.rm = TRUE)) %>%
    mutate(Gear = "Seine", "Species Richness" = nrow(Matrix_Seine),
           "Shannon Diversity" = exp(diversity(t(Matrix_Seine), groups = "none", index = "shannon")),
           "Simpson Diversity" = diversity(t(Matrix_Seine), groups = "none", index = "invsimpson")),
  # Cast Net & Centipede Net
  GearSpeciesGrid %>% select(`Cast Net`, `Centipede Net`) %>% summarize("Total Catch (n)" = sum(., na.rm = TRUE)) %>%
    mutate(Gear = "Cast Net & Centipede Net", "Species Richness" = nrow(Matrix_Cast_Cent),
           "Shannon Diversity" = exp(diversity(t(Matrix_Cast_Cent), groups = "none", index = "shannon")),
           "Simpson Diversity" = diversity(t(Matrix_Cast_Cent), groups = "none", index = "invsimpson")),
  # Cast Net & Seine
  GearSpeciesGrid %>% select(`Cast Net`, Seine) %>% summarize("Total Catch (n)" = sum(., na.rm = TRUE)) %>%
    mutate(Gear = "Cast Net & Seine", "Species Richness" = nrow(Matrix_Cast_Seine),
           "Shannon Diversity" = exp(diversity(t(Matrix_Cast_Seine), groups = "none", index = "shannon")),
           "Simpson Diversity" = diversity(t(Matrix_Cast_Seine), groups = "none", index = "invsimpson")),
  # Centipede Net & Seine
  GearSpeciesGrid %>% select(`Centipede Net`, Seine) %>% summarize("Total Catch (n)" = sum(., na.rm = TRUE)) %>%
    mutate(Gear = "Centipede Net & Seine", "Species Richness" = nrow(Matrix_Cent_Seine),
           "Shannon Diversity" = exp(diversity(t(Matrix_Cent_Seine), groups = "none", index = "shannon")),
           "Simpson Diversity" = diversity(t(Matrix_Cent_Seine), groups = "none", index = "invsimpson")),
  # All Gears (using the true study-wide objects)
  tibble(
    Gear = "All Gears",
    "Total Catch (n)" = sum(GearSpeciesGrid[,c("Cast Net", "Centipede Net", "Seine")], na.rm = TRUE),
    "Species Richness" = length(unique(na.omit(MergedData$Species))),
    "Shannon Diversity" = exp(diversity(t(StudyWide_Incidence_Matrix), groups = "none", index = "shannon")),
    "Simpson Diversity" = diversity(t(StudyWide_Incidence_Matrix), groups = "none", index = "invsimpson")
  )
)

# --- 2. Create the final table by merging MEAN and TOTAL stats ---
Tables$`Gear Study Summary` <- merge(
  bind_rows(
    SampleData %>%
      summarize(
        "Total Samples" = length(Richness), "Mean Effort" = mean(Effort),
        "Mean Catch (n)" = mean(Abundance, na.rm = T), "Mean CPUE" = mean(CPUE, na.rm = T),
        "Mean Biomass (g)" = mean(Biomass, na.rm = T), "Mean Richness" = mean(Richness, na.rm = T),
        "Mean Shannon Diversity" = mean(Shannon, na.rm = T), "Mean Simpson Diversity" = mean(Simpson, na.rm = T),
        .by = Gear
      ),
    SampleData %>%
      summarize(
        Gear = "All Gears", "Total Samples" = length(Richness), "Mean Effort" = NA,
        "Mean Catch (n)" = mean(Abundance, na.rm = T), "Mean CPUE" = NA,
        "Mean Biomass (g)" = mean(Biomass, na.rm = T), "Mean Richness" = mean(Richness, na.rm = T),
        "Mean Shannon Diversity" = mean(Shannon, na.rm = T), "Mean Simpson Diversity" = mean(Simpson, na.rm = T)
      )
  ),
  total_stats_final,
  by = "Gear", all = TRUE, sort = FALSE
) %>%
  mutate("Effort Unit" = c("Throws", "Net Group Hours", "Hauls", NA, NA, NA, NA), .after = "Mean Effort") %>% 
  arrange(match(Gear, c("Cast Net", "Centipede Net", "Seine", "Cast Net & Centipede Net", "Cast Net & Seine", "Centipede Net & Seine", "All Gears"))) %>%
  relocate(Gear, `Total Samples`, `Total Catch (n)`, `Species Richness`, `Shannon Diversity`, `Simpson Diversity`, `Mean Catch (n)`, `Mean Richness`, `Mean Shannon Diversity`, `Mean Simpson Diversity`, `Mean Biomass (g)`, `Mean CPUE`)


# Dominant species with abundances:
Tables$`Gear Species Abundance` <- cbind(
  SpeciesList %>%
    arrange(desc(Abundance)) %>%
    select("Species (All Gears)" = Species, "All Gears" = Abundance) %>%
    slice_head(n = 5),
  GearSpeciesGrid %>%
    select("Species (Cast Net)" = Species, `Cast Net`) %>%
    arrange(desc(`Cast Net`)) %>%
    slice_head(n = 5),
  GearSpeciesGrid %>%
    select("Species (Centipede Net)" = Species, `Centipede Net`) %>%
    arrange(desc(`Centipede Net`)) %>%
    slice_head(n = 5),
  GearSpeciesGrid %>%
    select("Species (Seine)" = Species, Seine) %>%
    arrange(desc(Seine)) %>%
    slice_head(n = 5)
)


# Unique species by Gear
# Count
GearComboGrid <- GearSpeciesGrid %>%
  mutate(GearCombo = case_when(
    `Cast Net` > 0 & `Centipede Net` > 0 & `Seine` > 0 ~ "All Gears",
    `Cast Net` > 0 & `Centipede Net` > 0 & `Seine` == 0 ~ "Cast & Centipede",
    `Cast Net` > 0 & `Centipede Net` == 0 & `Seine` > 0 ~ "Cast & Seine",
    `Cast Net` == 0 & `Centipede Net` > 0 & `Seine` > 0 ~ "Centipede & Seine",
    `Cast Net` > 0 & `Centipede Net` == 0 & `Seine` == 0 ~ "Cast Only",
    `Cast Net` == 0 & `Centipede Net` > 0 & `Seine` == 0 ~ "Centipede Only",
    `Cast Net` == 0 & `Centipede Net` == 0 & `Seine` > 0 ~ "Seine Only",
    .default = "??"
  ))

Tables$`Gear Combo Richness` <- GearComboGrid %>%
  summarize(Species = n(), .by = "GearCombo") %>%
  arrange(desc(Species))

# ID
Tables$`Gear Exclusives` <- bind_cols(
  GearComboGrid %>%
    filter(GearCombo == "Cast Only") %>%
    select(Species, "Cast Net") %>%
    rename("Cast Net" = Species, "(Cast)" = `Cast Net`) %>%
    bind_rows(tibble(.rows = (19 - length(.[, 1])))),
  GearComboGrid %>%
    filter(GearCombo == "Centipede Only") %>%
    select(Species, "Centipede Net") %>%
    rename("Centipede Net " = Species, "(Cent)" = `Centipede Net`) %>%
    bind_rows(tibble(.rows = (19 - length(.[, 1])))),
  GearComboGrid %>%
    filter(GearCombo == "Seine Only") %>%
    select(Species, "Seine") %>%
    rename("Seine" = Species, "(Seine)" = Seine) %>%
    bind_rows(tibble(.rows = (19 - length(.[, 1]))))
)

# Venn Diagram
venn_fit <- euler(GearSpeciesCheck[, c(2, 4, 3)] == "✓")
Gear_Venn <- plot(venn_fit,
                  labels = list(col = "black", fontfamily = "serif SemiBold", fontsize = 12),
                  quantities = list(col = "black", fontfamily = "serif", fontsize = 12),
                  edges = list(col = "black", alpha = 1),
                  fills = list(fill = GearColors()[c(2, 4, 3)], alpha = 0.9),
                  lty = c(2, 3, 1),
                  lwd = 3
                  )

Gear_Venn

# Saving as vector graphics
# cairo_ps(filename = "output/plots/Venn_Diagram_Final.eps", 
#          width = 5, 
#          height = 5, 
#          family = "serif")
# 
# plot(Gear_Venn)
# 
# dev.off()


# Create Excel file with a sheet for each table in Tables-------------
# Set minimum column width to avoid very narrow cells
options("openxlsx2.minWidth" = 15)
# Create workbook object
wb <- wb_workbook()

# Loop through each dataframe in Tables
for (i in 1:length(Tables)) {
  sheet_name <- names(Tables)[i]

  wb <- wb_add_worksheet(wb, sheet = sheet_name) %>%
    wb_add_data_table(wb,
      sheet = sheet_name, x = Tables[[i]],
      table_style = "TableStyleLight1",
      table_name = gsub(" ", "_", sheet_name),
      dims = wb_dims(from_col = 1),
      na.strings = ""
    ) %>%
    wb_set_col_widths(cols = 1:ncol(Tables[[1]]), widths = "auto") %>%
    wb_add_cell_style(dims = "A1:K100", horizontal = "center", vertical = "center", num_fmt_id = 4)
}

# Save the workbook (And add specific formatting tweaks for some tables)
wb <- wb %>%
  wb_add_font("Species List", dims = "A2:A100", italic = TRUE) %>%
  wb_add_cell_style("Species List", dims = "C2:E100", num_fmt_id = 1) %>%
  wb_save(file = "output/tables/Tables.xlsx", overwrite = TRUE)
