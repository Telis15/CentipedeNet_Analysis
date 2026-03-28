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
library(parallel)

# Load pre-wrangled dataframes and remove what I don't need right now
load("data/WrangledData.RData")
load("output/models/Rarefaction.RData")

# Source the global ggplot2 theme and aesthetic scales
source("scripts/0_PlotTheme.R")

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
    "Shannon Diversity" = exp(diversity(rowSums(StudyWide_Incidence_Matrix), index = "shannon")),
    "Simpson Diversity" = diversity(rowSums(StudyWide_Incidence_Matrix), index = "invsimpson"),
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
# # Species Accumulation Curves & Rarefaction with iNEXT --------------------
# #** These take 15-30 minutes to run**
# 
# # 1. Define the sets to analyze
# inext_sets <- list(
#   Gear_Out_Raw      = list(x = Incidence_Matrices, q = c(0,1,2)),
#   Gear_Out_Single   = list(x = Incidence_Matrices[1:3], q = 0),
#   Gear_Out_Combos   = list(x = Incidence_Matrices[4:7], q = 0),
#   Cast              = list(x = Incidence_Matrices['Cast Net'], q = 0, knots = length(Incidence_Matrices$`Cast Net`)),
#   Cent              = list(x = Incidence_Matrices['Centipede Net'], q = 0, knots = length(Incidence_Matrices$`Centipede Net`)),
#   Seine             = list(x = Incidence_Matrices['Seine'], q = 0, knots = length(Incidence_Matrices$Seine)),
#   Cast_Cent         = list(x = Incidence_Matrices['Cast Net & Centipede Net'], q = 0, knots = length(Incidence_Matrices$`Cast Net & Centipede Net`)),
#   Cast_Seine        = list(x = Incidence_Matrices['Cast Net & Seine'], q = 0, knots = length(Incidence_Matrices$`Cast Net & Seine`)),
#   Cent_Seine        = list(x = Incidence_Matrices['Centipede Net & Seine'], q = 0, knots = length(Incidence_Matrices$`Centipede Net & Seine`)),
#   All_Gears         = list(x = Incidence_Matrices['All Gears'], q = 0, knots = length(Incidence_Matrices$`All Gears`))
# )
# 
# # 2. Run in Parallel
# clu <- makeCluster(detectCores() - 1)
# clusterExport(clu, varlist = c("inext_sets", "Incidence_Matrices"))
# clusterEvalQ(clu, library(iNEXT))
# 
# inext_results <- parLapply(clu, inext_sets, function(params) {
#   do.call(iNEXT, c(params, list(datatype = "incidence_raw", nboot = 1000)))
# })
# 
# stopCluster(clu)
# 
# # 3. Extract results back to individual objects for plotting compatibility
# Gear_Out_Raw        <- inext_results$Gear_Out_Raw
# Gear_Out_Single_Raw <- inext_results$Gear_Out_Single
# Gear_Out_Combos_Raw <- inext_results$Gear_Out_Combos
# Cast                <- inext_results$Cast
# Cent                <- inext_results$Cent
# Seine               <- inext_results$Seine
# Cast_Cent           <- inext_results$Cast_Cent
# Cast_Seine          <- inext_results$Cast_Seine
# Cent_Seine          <- inext_results$Cent_Seine
# All_Gears           <- inext_results$All_Gears
# 
# # 4. Coverage estimates (separate call as it's quick once richness is known)
# SampleEst <- estimateD(Incidence_Matrices, q = 0, datatype = "incidence_raw",
#                         base = "coverage", level = c(0.75, .95), nboot = 1000) %>%
#   rename("Samples (n)" = t) %>%
#   mutate(SC = paste(100 * round(SC, digits = 2), "% Coverage")) %>%
#   select(-Order.q) %>%
#   arrange(SC) %>%
#   add_row(.after = 7)
# 
# # 5. Save the results
# save(
#   Gear_Out_Raw, Gear_Out_Single_Raw, Gear_Out_Combos_Raw,
#   Cast, Cent, Seine,
#   Cast_Cent, Cast_Seine, Cent_Seine, All_Gears,
#   SampleEst,
#   file = "output/models/Rarefaction.RData"
# )

# Sample Coverage Plots for Solo Gears & Combo Gears ----------------------
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
# 
# # Display the plot
# plot(Plot_Name)
# 
# # #If it's a grid you're lookin' for:
# Grid_Name <- plot_grid(Plot_Name1, Plot_Name2, nrow = 2, labels = Letters(), label_x = 0.11, align = "hv", label_fontface = "plain", label_size = 12, label_fontfamily = "serif")
# 
# # ggsave(plot = Plot_Name, "output/plots/Plot_Name.png", width = 6, height = 4, units = "in")


## Sample Coverage Plots for Solo Gears & Combo Gears ----------------------

# --- Streamlining: Create a common theme for both plots ---
common_theme <- function() {
  theme(
      panel.background = element_rect(fill = "transparent"),
      legend.background = element_rect(fill = "transparent"),
      legend.title = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(1, 0),
      legend.justification = c(1, 0),
      legend.key.spacing.y = unit(0.05, "in"),
      legend.key.height = unit(0.15, "in"),
      legend.key.width = unit(0.35, "in")
    )
}

# --- Singles_Coverage_Plot (Native ggplot2 Approach) ---

# 1. Extract and format data
df_singles <- fortify(Gear_Out_Single_Raw, type = 2) %>%
  mutate(Assemblage = factor(Assemblage, levels = c("Cast Net", "Centipede Net", "Seine")))

# 2. Build mapping vectors based on factor levels
single_names <- levels(df_singles$Assemblage)
single_colors <- GearColorsAll()[single_names]
single_shapes <- unlist(GearShapes[single_names])
single_lines  <- GearLines()[single_names]

# 3. Build the plot
Singles_Coverage_Plot <- ggplot(df_singles, aes(x = x, y = y, color = Assemblage, fill = Assemblage)) +
  # Confidence Ribbon
  geom_ribbon(aes(ymin = y.lwr, ymax = y.upr), alpha = 0.2, color = NA, show.legend = FALSE) +

  # Rarefaction Lines (Thick and Bold)
  geom_line(data = subset(df_singles, Method == "Rarefaction"), 
            aes(linetype = Assemblage), linewidth = 1, alpha = 1) +
  
  # Extrapolation Lines (Thin and slightly faded)
  geom_line(data = subset(df_singles, Method == "Extrapolation"), 
            aes(linetype = Assemblage), linewidth = 0.8, alpha = 0.75) +
  # Observed Points
  geom_point(data = subset(df_singles, Method == "Observed"), 
             aes(shape = Assemblage), size = 3) +
             
  # Scales
  scale_color_manual(values = single_colors) +
  scale_fill_manual(values = single_colors) +
  scale_shape_manual(values = single_shapes) +
  scale_linetype_manual(values = single_lines) +
  
  # 4. Use Guides to Force Linetypes into the Legend
  guides(
    color = guide_legend(
      keywidth = unit(0.6, "in"),
      override.aes = list(
        linetype = single_lines,
        shape = single_shapes,
        fill = single_colors,
        linewidth = 0.9,
        size = 2.5 
      )
    ),
    linetype = "none",
    shape = "none",
    fill = "none"
  ) +

  # Labels and Formatting
  labs(x = "Sampling Units", y = "Sample Coverage") +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01), labels = scales::percent) +
  common_theme()


# --- Combos_Coverage_Plot (Native ggplot2 Approach) ---

# 1. Extract and format data
df_combos <- fortify(Gear_Out_Combos_Raw, type = 2) %>%
  mutate(Assemblage = factor(Assemblage, levels = c("All Gears", "Cast Net & Centipede Net", "Cast Net & Seine", "Centipede Net & Seine")))

# 2. Build certain vectors for manual mapping to ensure fixed legend order
combo_names <- levels(df_combos$Assemblage)
combo_colors <- GearColorsAll()[combo_names]
combo_shapes <- unlist(GearShapes[combo_names])
combo_lines  <- GearLines()[combo_names]

# 3. Build the plot
Combos_Coverage_Plot <- ggplot(df_combos, aes(x = x, y = y, color = Assemblage, fill = Assemblage)) +
  # Confidence Ribbon
  # Using drastically lower alpha (0.08) to prevent muddying where ribbons overlap
  geom_ribbon(aes(ymin = y.lwr, ymax = y.upr), alpha = 0.08, color = NA, show.legend = FALSE) +
  
  # Rarefaction Lines (Thick and Bold)
  geom_line(data = subset(df_combos, Method == "Rarefaction"), 
            aes(linetype = Assemblage), linewidth = 1, alpha = 1) +
  
  # Extrapolation Lines (Thin and slightly faded)
  geom_line(data = subset(df_combos, Method == "Extrapolation"), 
            aes(linetype = Assemblage), linewidth = 0.8, alpha = 0.7) +
            
  # Observed Points
  geom_point(data = subset(df_combos, Method == "Observed"), 
             aes(shape = Assemblage), size = 3) +
             
  # Scales (Named vectors ensure correct mapping)
  scale_color_manual(values = combo_colors) +
  scale_fill_manual(values = combo_colors) +
  scale_shape_manual(values = combo_shapes) +
  scale_linetype_manual(values = combo_lines) +
  
  # 4. Use Guides to Force Linetypes into the Legend
  guides(
    color = guide_legend(
      keywidth = unit(0.6, "in"),
      override.aes = list(
        linetype = combo_lines,
        shape = combo_shapes,
        fill = combo_colors, # Ensure shapes are filled
        linewidth = 0.9,
        size = 2.5 
      )
    ),
    linetype = "none",
    shape = "none",
    fill = "none"
  ) +
  # Labels and Formatting
  labs(x = "Sampling Units", y = NULL) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01), labels = scales::percent, limits = c(.20, NA)) +
  common_theme()



# --- Combine Plots ---
Coverage_Grid <-
  plot_grid(Singles_Coverage_Plot, Combos_Coverage_Plot, nrow = 1, labels = c("Individual Gears", "Gear Combinations"), label_x = 0.05, align = "hv", rel_widths = c(1.03, 1), label_fontface = "plain", label_size = 9, label_fontfamily = "serif")

# Display the final plot
plot(Coverage_Grid)

# # --- Save the Final Plot for Publication ---
# ggsave(plot = Coverage_Grid,
#        filename = "output/plots/Coverage_2Panel_Final_Fig6.eps",
#        device = cairo_ps,
#        width = 5.62,
#        height = 3.75,
#        units = "in")
# ggsave(plot = Coverage_Grid,
#        filename = "output/plots/Coverage_2Panel_Final_Fig6.png",
#        width = 5.62,
#        height = 3.75,
#        units = "in",
#        dpi = 600)


# All q Plot | Diversity Accumulation Curve --------------------------------------------------------------

# 1. Extract and format data
df_all_q <- fortify(Gear_Out_Raw, type = 1) %>%
  mutate(Order.q = case_when(
    Order.q == 0 ~ "Species Richness",
    Order.q == 1 ~ "Shannon Diversity",
    Order.q == 2 ~ "Simpson Diversity"
  ),
  Order.q = factor(Order.q, levels = c("Species Richness", "Shannon Diversity", "Simpson Diversity")),
  # Force matching levels to ensure legend order (Performance Trend: All -> Combos -> Singles)
  # Added a blank spacer " " so All Gears can sit in its own column in a 2-row layout
  Assemblage = factor(Assemblage, levels = c(
    "All Gears", "Cast Net & Centipede Net", "Cast Net & Seine", "Centipede Net & Seine",
    " ", "Cast Net", "Centipede Net", "Seine"
  )))

# 2. Build mapping vectors based on factor levels
all_q_names <- levels(df_all_q$Assemblage)
all_q_colors <- c(GearColorsAll(), " " = "transparent")[all_q_names]
all_q_shapes <- c(unlist(GearShapes), " " = NA)[all_q_names]
all_q_lines  <- c(GearLines(), " " = "blank")[all_q_names]

# 3. Build the plot
All_q_Plot <- ggplot(df_all_q, aes(x = x, y = y, color = Assemblage, fill = Assemblage)) +
  # Confidence Ribbon
  # Using drastically lower alpha (0.08) to prevent muddying in the crowded diversity plot
  geom_ribbon(aes(ymin = y.lwr, ymax = y.upr), alpha = 0.08, color = NA, show.legend = FALSE) +
  
  # Rarefaction Lines (Thick and Bold)
  geom_line(data = subset(df_all_q, Method == "Rarefaction"), 
            aes(linetype = Assemblage), linewidth = 1, alpha = 1) +
  
  # Extrapolation Lines (Thin and slightly faded)
  geom_line(data = subset(df_all_q, Method == "Extrapolation"), 
            aes(linetype = Assemblage), linewidth = 0.8, alpha = 0.75) +
            
  # Observed Points
  geom_point(data = subset(df_all_q, Method == "Observed"), 
             aes(shape = Assemblage), size = 3) +
             
  # Scales
  scale_color_manual(values = all_q_colors, breaks = all_q_names, drop = FALSE) +
  scale_fill_manual(values = all_q_colors, breaks = all_q_names, drop = FALSE) +
  scale_shape_manual(values = all_q_shapes, breaks = all_q_names, drop = FALSE) +
  scale_linetype_manual(values = all_q_lines, breaks = all_q_names, drop = FALSE) +
  
  # 4. Use Guides to Force Linetypes into the Legend
  guides(
    color = guide_legend(
      nrow = 2,           # Two rows for a wide layout
      byrow = TRUE,       # Left-to-right, then top-to-bottom
      keywidth = unit(0.5, "in"),
      override.aes = list(
        linetype = all_q_lines,
        shape = all_q_shapes,
        fill = all_q_colors,
        linewidth = c(0.9, 0.9, 0.9, 0.9, 0, 0.9, 0.9, 0.9), # Hide line for spacer
        size = c(2.5, 2.5, 2.5, 2.5, 0, 2.5, 2.5, 2.5)       # Hide point for spacer
      )
    ),
    linetype = "none",
    shape = "none",
    fill = "none"
  ) +
  # Labels and Formatting
  labs(x = "Sampling Units", y = "Effective Species") +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  facet_wrap(~Order.q, scales = "free_y") +
  theme(
    panel.background = element_rect(fill = "transparent"),
    legend.background = element_rect(fill = "transparent"), 
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.justification = "center",
    strip.text = element_text(vjust = 0.5),
    legend.key.height = unit(0.15, "in"),
    legend.key.width = unit(0.5, "in")
  )

# Display the plot
plot(All_q_Plot)

# # --- Save the Final Plot for Publication --- Figure 9
# ggsave(plot = All_q_Plot,
#        filename = "output/plots/All_q_Plot_Final.eps",
#        device = cairo_ps,
#        width = 5.62,
#        height = 3.75,
#        units = "in")
# ggsave(plot = All_q_Plot,
#        filename = "output/plots/All_q_Plot_Final.png",
#        width = 5.62,
#        height = 3.75,
#        units = "in",
#        dpi = 600)



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
           "Shannon Diversity" = exp(diversity(rowSums(Matrix_Cast), index = "shannon")),
           "Simpson Diversity" = diversity(rowSums(Matrix_Cast), index = "invsimpson")),
  # Centipede Net
  GearSpeciesGrid %>% select(`Centipede Net`) %>% summarize("Total Catch (n)" = sum(`Centipede Net`, na.rm = TRUE)) %>%
    mutate(Gear = "Centipede Net", "Species Richness" = nrow(Matrix_Cent),
           "Shannon Diversity" = exp(diversity(rowSums(Matrix_Cent), index = "shannon")),
           "Simpson Diversity" = diversity(rowSums(Matrix_Cent), index = "invsimpson")),
  # Seine
  GearSpeciesGrid %>% select(Seine) %>% summarize("Total Catch (n)" = sum(Seine, na.rm = TRUE)) %>%
    mutate(Gear = "Seine", "Species Richness" = nrow(Matrix_Seine),
           "Shannon Diversity" = exp(diversity(rowSums(Matrix_Seine), index = "shannon")),
           "Simpson Diversity" = diversity(rowSums(Matrix_Seine), index = "invsimpson")),
  # Cast Net & Centipede Net
  GearSpeciesGrid %>% select(`Cast Net`, `Centipede Net`) %>% summarize("Total Catch (n)" = sum(., na.rm = TRUE)) %>%
    mutate(Gear = "Cast Net & Centipede Net", "Species Richness" = nrow(Matrix_Cast_Cent),
           "Shannon Diversity" = exp(diversity(rowSums(Matrix_Cast_Cent), index = "shannon")),
           "Simpson Diversity" = diversity(rowSums(Matrix_Cast_Cent), index = "invsimpson")),
  # Cast Net & Seine
  GearSpeciesGrid %>% select(`Cast Net`, Seine) %>% summarize("Total Catch (n)" = sum(., na.rm = TRUE)) %>%
    mutate(Gear = "Cast Net & Seine", "Species Richness" = nrow(Matrix_Cast_Seine),
           "Shannon Diversity" = exp(diversity(rowSums(Matrix_Cast_Seine), index = "shannon")),
           "Simpson Diversity" = diversity(rowSums(Matrix_Cast_Seine), index = "invsimpson")),
  # Centipede Net & Seine
  GearSpeciesGrid %>% select(`Centipede Net`, Seine) %>% summarize("Total Catch (n)" = sum(., na.rm = TRUE)) %>%
    mutate(Gear = "Centipede Net & Seine", "Species Richness" = nrow(Matrix_Cent_Seine),
           "Shannon Diversity" = exp(diversity(rowSums(Matrix_Cent_Seine), index = "shannon")),
           "Simpson Diversity" = diversity(rowSums(Matrix_Cent_Seine), index = "invsimpson")),
  # All Gears (using the true study-wide objects)
  tibble(
    Gear = "All Gears",
    "Total Catch (n)" = sum(GearSpeciesGrid[,c("Cast Net", "Centipede Net", "Seine")], na.rm = TRUE),
    "Species Richness" = length(unique(na.omit(MergedData$Species))),
    "Shannon Diversity" = exp(diversity(rowSums(StudyWide_Incidence_Matrix), index = "shannon")),
    "Simpson Diversity" = diversity(rowSums(StudyWide_Incidence_Matrix), index = "invsimpson")
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
        Gear = "All Gears", "Total Samples" = NA_real_, "Mean Effort" = NA_real_,
        "Mean Catch (n)" = mean(Abundance, na.rm = T), "Mean CPUE" = NA_real_,
        "Mean Biomass (g)" = mean(Biomass, na.rm = T), "Mean Richness" = mean(Richness, na.rm = T),
        "Mean Shannon Diversity" = mean(Shannon, na.rm = T), "Mean Simpson Diversity" = mean(Simpson, na.rm = T)
      )
  ),
  total_stats_final,
  by = "Gear"
) %>%
  mutate(`Total Samples` = case_when(
    Gear == "Cast Net & Centipede Net" ~ SampleData %>% filter(Gear %in% c("Cast Net", "Centipede Net")) %>% group_by(Site, Visit) %>% filter(n() == 2) %>% ungroup() %>% distinct(Site, Visit) %>% nrow(),
    Gear == "Cast Net & Seine" ~ SampleData %>% filter(Gear %in% c("Cast Net", "Seine")) %>% group_by(Site, Visit) %>% filter(n() == 2) %>% ungroup() %>% distinct(Site, Visit) %>% nrow(),
    Gear == "Centipede Net & Seine" ~ SampleData %>% filter(Gear %in% c("Centipede Net", "Seine")) %>% group_by(Site, Visit) %>% filter(n() == 2) %>% ungroup() %>% distinct(Site, Visit) %>% nrow(),
    Gear == "All Gears" ~ SampleData %>% filter(Gear %in% c("Cast Net", "Centipede Net", "Seine")) %>% group_by(Site, Visit) %>% filter(n() == 3) %>% ungroup() %>% distinct(Site, Visit) %>% nrow(),
    TRUE ~ `Total Samples`
  )) %>%
  mutate("Effort Unit" = c("Throws", "Net Group Hours", "Hauls", NA, NA, NA, NA), .after = "Mean Effort") %>% 
  arrange(match(Gear, c("Cast Net", "Centipede Net", "Seine", "Cast Net & Centipede Net", "Cast Net & Seine", "Centipede Net & Seine", "All Gears"))) %>%
  relocate(Gear, `Total Samples`, `Total Catch (n)`, `Species Richness`, `Shannon Diversity`, `Simpson Diversity`, `Mean Catch (n)`, `Mean Richness`, `Mean Shannon Diversity`, `Mean Simpson Diversity`, `Mean Biomass (g)`, `Mean CPUE`)


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
  theme(
    legend.position = "none", 
    axis.text = element_text(size = 9, color = "black")
  )

# Display the plot
plot(Violin_Plot_SL)

# # --- 3. Save the plot in the correct format for publication ---
# ggsave(plot = Violin_Plot_SL,
#        filename = "output/plots/Length_Gear_Violin_SL.eps",
#        device = cairo_ps,
#        width = 2.75,
#        height = 3.66,
#        units = "in")
# ggsave(plot = Violin_Plot_SL,
#        filename = "output/plots/Length_Gear_Violin_SL.png",
#        width = 2.75,
#        height = 3.66,
#        units = "in",
#        dpi = 600)


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
                  labels = list(labels = c("Cast Net", "Seine", "Centipede Net"), 
                                col = "black", fontfamily = "serif", fontface = "semibold", fontsize = 9),
                  quantities = list(col = "black", fontfamily = "serif", fontsize = 9),
                  edges = list(col = "black", alpha = 1),
                  fills = list(fill = GearColors()[c(2, 4, 3)], alpha = 0.7),
                  lty = c("22", "1343", "73"), # Cast Net, Seine, Centipede Net
                  lwd = 1.5
                  )

Gear_Venn

# # Saving as vector graphics
# cairo_ps(filename = "output/plots/Venn_Diagram_Final.eps",
#          width = 5,
#          height = 5,
#          family = "serif")
# 
# plot(Gear_Venn)
# 
# dev.off()
# ggsave(plot = Gear_Venn, "output/plots/Venn_Diagram_Final.png", 
#        width = 5, height = 5, units = "in", dpi = 600)

# Create Excel file with a sheet for each table in Tables-------------
# Set minimum column width to avoid very narrow cells
options("openxlsx2.minWidth" = 15)
# Create workbook object
wb <- wb_workbook()

# Loop through each dataframe in Tables
for (i in 1:length(Tables)) {
  sheet_name <- names(Tables)[i]

  wb <- wb %>%
    wb_add_worksheet(sheet = sheet_name) %>%
    wb_add_data_table(
      sheet = sheet_name, x = Tables[[i]],
      table_style = "TableStyleLight1",
      table_name = gsub(" ", "_", sheet_name),
      dims = wb_dims(from_col = 1),
      na.strings = ""
    ) %>%
    wb_set_col_widths(sheet = sheet_name, cols = 1:ncol(Tables[[i]]), widths = "auto") %>%
    wb_add_cell_style(sheet = sheet_name, dims = wb_dims(x = Tables[[i]]), horizontal = "center", vertical = "center", num_fmt_id = 4)
}

# Save the workbook (And add specific formatting tweaks for some tables)
wb <- wb %>%
  wb_add_font("Species List", dims = "A2:A100", italic = TRUE) %>%
  wb_add_cell_style("Species List", dims = "C2:E100", num_fmt_id = 1) %>%
  wb_save(file = "output/tables/Tables.xlsx", overwrite = TRUE)

 