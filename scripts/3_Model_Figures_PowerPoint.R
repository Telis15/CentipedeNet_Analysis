source('scripts/00_AI_Export_Utils.R')
# 3_Model_Figures_PowerPoint.R
# Streamlined script for generating Marginal Effects plots for PowerPoint.
# Loads pre-run models from output/models/DredgeListGLMM.RData.
# All output targets a 16:9 landscape aspect ratio (13.33 x 7.5 in).

# 1. Libraries and Data Loading -------------------------------------------
library("glmmTMB")
library("tidyverse")
library("cowplot")
library("modelbased")
library("insight")
library("MuMIn")
library("emmeans")
library("scales")

# Load pre-saved models and data
load("output/models/DredgeListGLMM.RData")

# Ensure ScaleDataEffort is available (may be missing from older .RData)
if (!exists("ScaleDataEffort")) {
  message("ScaleDataEffort missing; recalculating from geardat...")
  ScaleDataEffort <- geardat |>
    group_by(Gear) |>
    summarize(
      Mean = mean(log(Effort), na.rm = TRUE),
      SD = sd(log(Effort), na.rm = TRUE)
    )
}

# Source the global ggplot2 theme and aesthetic scales
source("scripts/0_PlotTheme.R")

# Apply the presentation theme
theme_set(presentation_theme)


# 2. Helper Functions -----------------------------------------------------

create_plot_data <- function(model, global_model, predictor, grouping_var = NULL) {
  by_vars <- if (is.null(grouping_var)) predictor else c(predictor, grouping_var)
  grid_args <- list(global_model, by = by_vars, preserve_range = FALSE)
  if (predictor != "Steepness") grid_args$Steepness <- "Medium"
  if ("logEffort" %in% names(coef(model)) && predictor != "logEffort") {
    grid_args$logEffort <- 0
  }
  newdata <- do.call(insight::get_datagrid, grid_args)
  newdata$Site <- NA
  
  preds <- predict(model, newdata = newdata, se.fit = TRUE, type = "link")
  plot_data <- data.frame(
    newdata,
    Predicted = exp(preds$fit),
    CI_low    = exp(preds$fit - 1.96 * preds$se.fit),
    CI_high   = exp(preds$fit + 1.96 * preds$se.fit)
  )
  
  if (is.null(grouping_var)) {
    plot_data <- plot_data |>
      group_by(!!sym(predictor)) |>
      summarize(across(c(Predicted, CI_low, CI_high), \(x) mean(x, na.rm = TRUE)))
  }
  return(plot_data)
}


# 3. Shared Legend --------------------------------------------------------

shared_legend_horiz <- get_legend(
  ggplot(geardat, aes(x = Gear, y = Abundance, color = Gear, shape = Gear, linetype = Gear, fill = Gear)) +
    geom_ribbon(aes(ymin = 0, ymax = 0), alpha = 0.2) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    scale_color_manual(values = GearColors()) +
    scale_fill_manual(values = GearColors()) +
    scale_shape_manual(values = GearShapes) +
    scale_linetype_manual(values = GearLines()) +
    guides(color = guide_legend()) +
    theme(legend.position = "top", legend.title = element_blank(),
          legend.key.width = rel(2), legend.direction = "horizontal",
          legend.text = element_text(size = 16))
)


# =========================================================================
# 4. ABUNDANCE Marginal Effects
# =========================================================================
# Abundance uses nbinom2 → log-scale y-axis is appropriate.

# --- Daylight x Gear (interaction) ---
p_abund_daylight <- estimate_relation(
  model = AvgMods$AbundInt,
  by = c("DaylightPercent = [fivenum]", "Gear"),
  fixed = list(logEffort = 0, Steepness = "Medium"),
  preserve_range = FALSE
) |>
  ggplot(aes(x = DaylightPercent, y = Predicted, color = Gear, fill = Gear, linetype = Gear)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.15, linetype = 0) +
  labs(title = "Daylight x Gear", y = "Abundance", x = "Percent Daylight") +
  scale_x_continuous(
    breaks = (c(0, 25, 50, 75, 100) - ScaleData$DaylightPercent$Mean) / ScaleData$DaylightPercent$SD,
    labels = c("0%", "25%", "50%", "75%", "100%")
  ) +
  gear_scales

# --- Effort x Gear (overlaid on relative axis) ---
newdata_effort_abund <- insight::get_datagrid(
  AvgMods$AbundInt, by = c("logEffort = [fivenum]", "Gear"),
  preserve_range = FALSE, Steepness = "Medium"
)
preds_effort_abund <- predict(AvgMods$AbundInt, newdata = newdata_effort_abund, se.fit = TRUE, type = "link")
plot_data_effort_abund <- data.frame(
  newdata_effort_abund,
  Predicted = exp(preds_effort_abund$fit),
  CI_low    = exp(preds_effort_abund$fit - 1.96 * preds_effort_abund$se.fit),
  CI_high   = exp(preds_effort_abund$fit + 1.96 * preds_effort_abund$se.fit)
)

p_abund_effort <- ggplot(plot_data_effort_abund, aes(x = logEffort, y = Predicted, color = Gear, fill = Gear, linetype = Gear)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.1, linetype = 0) +
  labs(title = "Effort x Gear", x = "Relative Effort (SD)", y = "Abundance") +
  scale_y_log10(labels = label_number(accuracy = 1)) +
  gear_scales

# --- Occlusion (main effect) ---
plot_data_abund_occ <- create_plot_data(AvgMods$AbundInt, AbundIntFit, "Occlusion")
p_abund_occ <- ggplot(plot_data_abund_occ, aes(x = Occlusion, y = Predicted)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.2, fill = "grey50") +
  labs(title = "Occlusion", x = "Percent Occlusion", y = "Abundance") +
  scale_x_continuous(
    breaks = (c(17, 22, 27) - ScaleData$Occlusion$Mean) / ScaleData$Occlusion$SD,
    labels = c("17%", "22%", "27%")
  )

# --- Season (main effect) ---
plot_data_abund_season <- create_plot_data(AvgMods$AbundInt, AbundIntFit, "Season")
p_abund_season <- ggplot(plot_data_abund_season, aes(x = Season, y = Predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +
  labs(title = "Season", x = "Season", y = "Abundance")

# --- Assemble: interactions on top, legend, main effects below ---
abund_top_row <- plot_grid(p_abund_daylight, p_abund_effort, ncol = 2, align = "hv")
abund_bot_row <- plot_grid(p_abund_occ, p_abund_season, ncol = 2, align = "hv")
final_plot_abund <- plot_grid(
  abund_top_row, shared_legend_horiz, abund_bot_row,
  ncol = 1, rel_heights = c(1, 0.08, 1)
) + theme(plot.background = element_rect(fill = "white", color = NA))


# =========================================================================
# 5. RICHNESS Marginal Effects
# =========================================================================
# Richness uses Poisson → natural scale y-axis is fine (counts are small).

# --- Effort x Gear (overlaid) ---
newdata_effort_rich <- insight::get_datagrid(
  AvgMods$RichInt, by = c("logEffort = [fivenum]", "Gear"),
  preserve_range = FALSE, Steepness = "Medium"
)
newdata_effort_rich$Site <- NA
preds_effort_rich <- predict(AvgMods$RichInt, newdata = newdata_effort_rich, se.fit = TRUE, type = "link")
plot_data_effort_rich <- data.frame(
  newdata_effort_rich,
  Predicted = exp(preds_effort_rich$fit),
  CI_low    = exp(preds_effort_rich$fit - 1.96 * preds_effort_rich$se.fit),
  CI_high   = exp(preds_effort_rich$fit + 1.96 * preds_effort_rich$se.fit)
)

p_rich_effort <- ggplot(plot_data_effort_rich, aes(x = logEffort, y = Predicted, color = Gear, fill = Gear, linetype = Gear)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.1, linetype = 0) +
  labs(title = "Effort x Gear", x = "Relative Effort (SD)", y = "Richness") +
  gear_scales

# --- Season ---
plot_data_rich_season <- create_plot_data(AvgMods$RichInt, RichIntFit, "Season")
p_rich_season <- ggplot(plot_data_rich_season, aes(x = Season, y = Predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +
  labs(title = "Season", x = "Season", y = "Richness") +
  scale_y_continuous(breaks = c(3, 5, 7, 9), minor_breaks = NULL) +
  coord_cartesian(ylim = c(2, 10))

# --- Daylight ---
plot_data_rich_daylight <- create_plot_data(AvgMods$RichInt, RichIntFit, "DaylightPercent")
p_rich_daylight <- ggplot(plot_data_rich_daylight, aes(x = DaylightPercent, y = Predicted)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.2, fill = "grey50") +
  labs(title = "Daylight", x = "Percent Daylight", y = "Richness") +
  scale_x_continuous(
    breaks = (c(0, 25, 50, 75, 100) - ScaleData$DaylightPercent$Mean) / ScaleData$DaylightPercent$SD,
    labels = c("0%", "25%", "50%", "75%", "100%")
  ) +
  scale_y_continuous(breaks = c(3, 5, 7, 9), minor_breaks = NULL) +
  coord_cartesian(ylim = c(2, 10))

# --- Steepness ---
plot_data_rich_steepness <- create_plot_data(AvgMods$RichInt, RichIntFit, "Steepness")
p_rich_steepness <- ggplot(plot_data_rich_steepness, aes(x = Steepness, y = Predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +
  labs(title = "Steepness", x = "Steepness", y = "Richness") +
  scale_x_discrete(limits = c("Low", "Medium", "High")) +
  scale_y_continuous(breaks = c(3, 5, 7, 9), minor_breaks = NULL) +
  coord_cartesian(ylim = c(2, 10))

# --- Occlusion ---
plot_data_rich_occ <- create_plot_data(AvgMods$RichInt, RichIntFit, "Occlusion")
p_rich_occ <- ggplot(plot_data_rich_occ, aes(x = Occlusion, y = Predicted)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.2, fill = "grey50") +
  labs(title = "Occlusion", x = "Percent Occlusion", y = "Richness") +
  scale_x_continuous(
    breaks = (c(17, 22, 27) - ScaleData$Occlusion$Mean) / ScaleData$Occlusion$SD,
    labels = c("17%", "22%", "27%")
  ) +
  scale_y_continuous(breaks = c(3, 5, 7, 9), minor_breaks = NULL) +
  coord_cartesian(ylim = c(2, 10))

# --- Substrate ---
plot_data_rich_mud <- create_plot_data(AvgMods$RichInt, RichIntFit, "MudDominant")
p_rich_mud <- ggplot(plot_data_rich_mud, aes(x = MudDominant, y = Predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +
  labs(title = "Substrate", x = "Substrate", y = "Richness") +
  scale_x_discrete(labels = c("False" = "Other", "True" = "Mud"), limits = rev) +
  scale_y_continuous(breaks = c(3, 5, 7, 9), minor_breaks = NULL) +
  coord_cartesian(ylim = c(2, 10))

# --- Secchi Depth ---
plot_data_rich_secchi <- create_plot_data(AvgMods$RichInt, RichIntFit, "SecchiDepth")
p_rich_secchi <- ggplot(plot_data_rich_secchi, aes(x = SecchiDepth, y = Predicted)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.2, fill = "grey50") +
  labs(title = "Secchi Depth", x = "Secchi Depth (cm)", y = "Richness") +
  scale_x_continuous(
    breaks = (c(20, 75, 125, 175) - ScaleData$SecchiDepth$Mean) / ScaleData$SecchiDepth$SD,
    labels = c("20", "75", "125", "175")
  ) +
  scale_y_continuous(breaks = c(3, 5, 7, 9), minor_breaks = NULL) +
  coord_cartesian(ylim = c(2, 10))

# --- Assemble: 4-panel top row, legend, 3-panel bottom row (full width) ---
rich_top_row <- plot_grid(p_rich_season, p_rich_steepness, p_rich_mud, p_rich_daylight, ncol = 4, align = "hv", axis = "tblr")
rich_bot_row <- plot_grid(p_rich_effort, p_rich_secchi, p_rich_occ, ncol = 3, align = "hv", axis = "tblr")
final_plot_rich <- plot_grid(
  rich_top_row, shared_legend_horiz, rich_bot_row,
  ncol = 1, rel_heights = c(1, 0.08, 1)
) + theme(plot.background = element_rect(fill = "white", color = NA))


# =========================================================================
# 6. SHANNON DIVERSITY Marginal Effects
# =========================================================================
# Shannon uses Tweedie → natural scale y-axis.

# --- Gear ---
plot_data_div_gear <- create_plot_data(AvgMods$DivInt, DivIntFit, "Gear")
p_div_gear <- ggplot(plot_data_div_gear, aes(x = Gear, y = Predicted, color = Gear)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +
  labs(title = "Gear", x = NULL, y = "Shannon") +
  gear_scales +
  scale_x_discrete(labels = c("Cast Net" = "Cast\nNet", "Centipede Net" = "Centipede\nNet", "Seine" = "Seine"))

# --- Daylight ---
plot_data_div_daylight <- create_plot_data(AvgMods$DivInt, DivIntFit, "DaylightPercent")
p_div_daylight <- ggplot(plot_data_div_daylight, aes(x = DaylightPercent, y = Predicted)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.2, fill = "grey50") +
  labs(title = "Daylight", x = "Percent Daylight", y = "Shannon") +
  scale_x_continuous(
    breaks = (c(0, 25, 50, 75, 100) - ScaleData$DaylightPercent$Mean) / ScaleData$DaylightPercent$SD,
    labels = c("0%", "25%", "50%", "75%", "100%")
  )

# --- Effort ---
plot_data_div_effort <- create_plot_data(AvgMods$DivInt, DivIntFit, "logEffort")
p_div_effort <- ggplot(plot_data_div_effort, aes(x = logEffort, y = Predicted)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.2, fill = "grey50") +
  labs(title = "Effort", x = "Relative Effort (SD)", y = "Shannon")

# --- Occlusion ---
plot_data_div_occ <- create_plot_data(AvgMods$DivInt, DivIntFit, "Occlusion")
p_div_occ <- ggplot(plot_data_div_occ, aes(x = Occlusion, y = Predicted)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.2, fill = "grey50") +
  labs(title = "Occlusion", x = "Percent Occlusion", y = "Shannon") +
  scale_x_continuous(
    breaks = (c(17, 22, 27) - ScaleData$Occlusion$Mean) / ScaleData$Occlusion$SD,
    labels = c("17%", "22%", "27%")
  )

# --- Steepness ---
plot_data_div_steepness <- create_plot_data(AvgMods$DivInt, DivIntFit, "Steepness")
p_div_steepness <- ggplot(plot_data_div_steepness, aes(x = Steepness, y = Predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +
  labs(title = "Steepness", x = "Steepness", y = "Shannon") +
  scale_x_discrete(limits = c("Low", "Medium", "High"))

# --- Substrate ---
plot_data_div_mud <- create_plot_data(AvgMods$DivInt, DivIntFit, "MudDominant")
p_div_mud <- ggplot(plot_data_div_mud, aes(x = MudDominant, y = Predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +
  labs(title = "Substrate", x = "Substrate", y = "Shannon") +
  scale_x_discrete(labels = c("False" = "Other", "True" = "Mud"), limits = rev)

# --- Assemble: 3 cols x 2 rows, landscape ---
final_plot_div <- plot_grid(
  p_div_daylight, p_div_effort, p_div_occ,
  p_div_gear, p_div_steepness, p_div_mud,
  ncol = 3, align = "hv"
)


# =========================================================================
# 7. SIMPSON DIVERSITY Marginal Effects
# =========================================================================
# Simpson uses Tweedie → natural scale y-axis.

# --- Gear ---
plot_data_sim_gear <- create_plot_data(AvgMods$SimInt, SimIntFit, "Gear")
p_sim_gear <- ggplot(plot_data_sim_gear, aes(x = Gear, y = Predicted, color = Gear)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +
  labs(title = "Gear", x = NULL, y = "Simpson") +
  gear_scales +
  scale_x_discrete(labels = c("Cast Net" = "Cast\nNet", "Centipede Net" = "Centipede\nNet", "Seine" = "Seine"))

# --- Daylight ---
plot_data_sim_daylight <- create_plot_data(AvgMods$SimInt, SimIntFit, "DaylightPercent")
p_sim_daylight <- ggplot(plot_data_sim_daylight, aes(x = DaylightPercent, y = Predicted)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.2, fill = "grey50") +
  labs(title = "Daylight", x = "Percent Daylight", y = "Simpson") +
  scale_x_continuous(
    breaks = (c(0, 25, 50, 75, 100) - ScaleData$DaylightPercent$Mean) / ScaleData$DaylightPercent$SD,
    labels = c("0%", "25%", "50%", "75%", "100%")
  )

# --- Assemble: 1x2 landscape ---
final_plot_sim <- plot_grid(p_sim_daylight, p_sim_gear, ncol = 2, align = "hv", rel_widths = c(1.2, 1))


# =========================================================================
# 8. SUMMARY 2x2 Gear Grid
# =========================================================================

# Gear-only plots for Abundance & Richness (for the summary grid)
plot_data_abund_gear <- create_plot_data(AvgMods$AbundInt, AbundIntFit, "Gear")
p_abund_gear <- ggplot(plot_data_abund_gear, aes(x = Gear, y = Predicted, color = Gear)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +
  labs(title = "Abundance", x = NULL, y = "Abundance") +
  scale_y_log10(labels = label_number(accuracy = 1)) +
  gear_scales +
  scale_x_discrete(labels = c("Cast Net" = "Cast\nNet", "Centipede Net" = "Centipede\nNet", "Seine" = "Seine"))

plot_data_rich_gear <- create_plot_data(AvgMods$RichInt, RichIntFit, "Gear")
p_rich_gear <- ggplot(plot_data_rich_gear, aes(x = Gear, y = Predicted, color = Gear)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +
  labs(title = "Richness", x = NULL, y = "Richness") +
  gear_scales +
  scale_x_discrete(labels = c("Cast Net" = "Cast\nNet", "Centipede Net" = "Centipede\nNet", "Seine" = "Seine"))

# Relabel for summary context
p_div_gear_summary <- p_div_gear + labs(title = "Shannon Diversity")
p_sim_gear_summary <- p_sim_gear + labs(title = "Simpson Diversity")

row1 <- plot_grid(p_abund_gear, p_rich_gear, ncol = 2, align = "hv")
row2 <- plot_grid(p_div_gear_summary, p_sim_gear_summary, ncol = 2, align = "hv")

final_metric_gear_2x2 <- plot_grid(
  row1, row2,
  ncol = 1
)


# =========================================================================
# 9. Save Output (all landscape 16:9)
# =========================================================================
if (!dir.exists("output/plots/PowerPoint")) dir.create("output/plots/PowerPoint")

# Standard PPT slide content area ≈ 13.33 x 7.5 in
ggsave(plot = final_plot_abund,       "output/plots/PowerPoint/MarginalEffects_Abundance_PPT.png", width = 13.33, height = 7.5, units = "in", dpi = 300)
ggsave(plot = final_plot_rich,        "output/plots/PowerPoint/MarginalEffects_Richness_PPT.png",  width = 16, height = 8.5,  units = "in", dpi = 300)
ggsave(plot = final_plot_div,         "output/plots/PowerPoint/MarginalEffects_Shannon_PPT.png",   width = 13.33, height = 7.5, units = "in", dpi = 300)
ggsave(plot = final_plot_sim,         "output/plots/PowerPoint/MarginalEffects_Simpson_PPT.png",   width = 13.33, height = 5,   units = "in", dpi = 300)
ggsave(plot = final_metric_gear_2x2,  "output/plots/PowerPoint/Metric_Gear_2x2_PPT.png",          width = 13.33, height = 7.5, units = "in", dpi = 300)

message("Done! Figures saved to output/plots/PowerPoint/")
