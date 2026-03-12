# AFS Global Plot Theme Settings (Times Roman, 9pt text, simple background)
# NOTE: This script must be sourced AFTER `data/WrangledData.RData` is loaded.

# Set the global ggplot2 theme
theme_set(
  theme_classic(
    base_family = "serif", 
    base_size = 9
  ) +
    theme(
      text = element_text(colour = "black")
    )
)

# Set the global discrete color and fill palettes to Okabe-Ito (colorblind safe)
options(
  ggplot2.discrete.colour = palette.colors(palette = "Okabe-Ito"),
  ggplot2.discrete.fill = palette.colors(palette = "Okabe-Ito")
)

# --- Visual Helper Functions ---

# Standard Gear Colors (Okabe-Ito)
GearColors <- function() {
  return(c("All Gears" = "#000000", "Cast Net" = "#E69F00", "Centipede Net" = "#56B4E9", "Seine" = "#009E73"))
}

# Full Gear Combinations Palette
GearColorsAll <- function() {
  return(c(
    "All Gears" = "#000000",                 # Black
    "Cast Net" = "#E69F00",                  # Orange
    "Cast Net & Centipede Net" = "#CC79A7",  # Pink/Rose
    "Cast Net & Seine" = "#785EF0",          # Violet
    "Centipede Net" = "#56B4E9",             # Sky Blue
    "Centipede Net & Seine" = "#B22222",     # Red
    "Seine" = "#009E73"                     # Green
  ))
}

GearShapes <- list(
  "All Gears" = 18,               # Filled diamond
  "Cast Net" = 16,                # Filled circle
  "Centipede Net" = 15,           # Filled square
  "Seine" = 17,                   # Filled triangle
  "Cast Net & Centipede Net" = 22, # Open square with border 
  "Cast Net & Seine" = 23,        # Open diamond with border 
  "Centipede Net & Seine" = 24    # Open triangle with border
)

GearLines <- function() {
  return(c(
    "Cast Net" = "22",          # "dashed"
    "Centipede Net" = "73",     # "longdash"
    "Seine" = "1343",           # "dotdash"
    "Cast Net & Centipede Net" = "1446",
    "Cast Net & Seine" = "491549",
    "Centipede Net & Seine" = "188888",
    "All Gears" = "solid"
  ))
}

SiteColors <- function() {
  palette.colors(6, "Tableau")[1:6]
}

# Multi-panel labeling helper (No periods/brackets per AFS)
Letters <- function() {
  return(c("A", "B", "C", "D", "E", "F", "G", "H", "I"))
}

# --- Shared Plotting Logic ---

# A list of common aesthetic scales for Gear mapping
gear_scales <- list(
  scale_color_manual(values = GearColors()),
  scale_fill_manual(values = GearColors()),
  scale_linetype_manual(values = GearLines()),
  scale_shape_manual(values = GearShapes)
)

# A reusable theme mapping for standard plots that hides redundant legends and manages margins
plot_theme <- 
  theme(legend.position = "none", 
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"))
