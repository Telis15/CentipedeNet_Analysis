source('scripts/00_AI_Export_Utils.R')
# This script handles ONLY the processing of occlusion images using countcolors.

library(tidyverse)
# library(scales)
# library(countcolors)

# Color matrix from eyedropper tool color sampling of yellow background in a variety of images to account for different lighting & shadows.
# colorMatrix <- matrix(ncol = 3, c(
#     0.533333333333333, 0.376470588235294, 0, 
#     ...
#     .949,.898,.18), byrow = T)
# colnames(colorMatrix) <- c("R", "G", "B")
# 
# 
# #Testing with newly sampled colors
# ColorCount <- countColorsInDirectory(".", color.range = "spherical", center = colorMatrix, radius = rep(.1, length(colorMatrix[,1])), target.color = "magenta", return.indicator = T, save.indicator = "./Masked")

#setwd("..")
# save.image("data/OcclusionData.RData")

# OcclDF <- map_dfr(ColorCount, ~tibble(Openness = .x[2])) |>
#     mutate(Image = names(ColorCount)) |>
#     mutate(Occlusion = 100*(1 - as.numeric(Openness))) |>
#     mutate(Net = str_sub_all(Image, end = -2)) |>
#     separate(Net, into = c("Site", "Group"), sep = "_") |>
#     mutate(Gear = "Centipede Net") |>
#     mutate(.by = c(Site, Gear, Group), OcclusionGroupAvg = mean(Occlusion), OcclusionGroupSTDev = sd(Occlusion)) |>
#     mutate(.by = Site, OcclusionSiteAvg = mean(Occlusion), OcclusionSiteSTDev = sd(Occlusion)) |>
#     select(Site, Gear, Group, OcclusionGroupAvg, OcclusionGroupSTDev, OcclusionSiteAvg, OcclusionSiteSTDev) |>
#     unique()
# # rm(ColorCount, colorMatrix)
# 
# save(OcclDF, file = "data/Occlusion.RData")

load("data/Occlusion.RData")

# Export Occlusion data frame to AI_Ready_Data
write.csv(OcclDF, "output/AI_Ready_Data/OcclDF.csv", row.names = FALSE)
