# This script handles ONLY the processing of occlusion images using countcolors.

 
library(tidyverse)
# library(scales)
library(countcolors)
# 
# Color matrix from eyedropper tool color sampling of yellow background in a variety of images to account for different lighting & shadows.
colorMatrix <- matrix(ncol = 3, c(
    0.533333333333333, 0.376470588235294, 0, 
    0.56078431372549, 0.345098039215686, 0, 
    0.635294117647059, 0.454901960784314, 0, 
    0.682352941176471, 0.462745098039216, 0, 
    0.709803921568627, 0.482352941176471, 0, 
    0.725490196078431, 0.509803921568627, 0, 
    0.811764705882353, 0.607843137254902, 0, 
    0.823529411764706, 0.635294117647059, 0, 
    0.870588235294118, 0.701960784313725, 0, 
    0.905882352941176, 0.733333333333333, 0, 
    0.913725490196078, 0.788235294117647, 0.431372549019608, 
    0.941176470588235, 0.729411764705882, 0, 
    0.941176470588235, 0.8, 0.00392156862745098, 
    0.964705882352941, 0.823529411764706, 0, 
    0.996078431372549, 0.87843137254902, 0.372549019607843, 
    1, 0.670588235294118, 0, 1, 
    0.768627450980392, 0.00392156862745098, 1, 
    0.803921568627451, 0, 1, 
    0.83921568627451, 0, 1, 
    0.866666666666667, 0.0784313725490196, 
    0.949, 0.886, 0.02, 
    .949,.914,.42,
    .949,.663,.133,
    .949,.898,.18), byrow = T)
# colnames(colorMatrix) <- c("R", "G", "B")
# 
# 
# #Testing with newly sampled colors
# ColorCount <- countColorsInDirectory(".", color.range = "spherical", center = colorMatrix, radius = rep(.1, length(colorMatrix[,1])), target.color = "magenta", return.indicator = T, save.indicator = "./Masked")

#setwd("..")
# save.image("data/OcclusionData.RData")

OcclDF <- map_dfr(ColorCount, ~tibble(Openness = .x[2])) %>%
    mutate(Image = names(ColorCount)) %>%
    mutate(Occlusion = 100*(1 - as.numeric(Openness))) %>%
    mutate(Net = str_sub_all(Image, end = -2)) %>%
    separate(Net, into = c("Site", "Group"), sep = "_") %>%
    mutate(Gear = "Centipede Net") %>%
    mutate(.by = c(Site, Gear, Group), OcclusionGroupAvg = mean(Occlusion), OcclusionGroupSTDev = sd(Occlusion)) %>%
    mutate(.by = Site, OcclusionSiteAvg = mean(Occlusion), OcclusionSiteSTDev = sd(Occlusion)) %>%
    select(Site, Gear, Group, OcclusionGroupAvg, OcclusionGroupSTDev, OcclusionSiteAvg, OcclusionSiteSTDev) %>%
    unique()
# rm(ColorCount, colorMatrix)

save(OcclDF, file = "data/Occlusion.RData")
