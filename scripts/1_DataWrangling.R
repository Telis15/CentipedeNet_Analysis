
#**Verified to be working on 3/6/2026, with changes made to utilize RProject and containerization to via GitHub**

# Packages ----------------------------------------------------------------
library(suncalc)
library(rfishbase)
library(vegan)
library(tidyverse)
library(openxlsx2)
# load("./WrangledData.RData")
# NOTE: Plotting helper functions (GearColors, GearShapes, Letters, etc.) 
# have been moved to scripts/0_PlotTheme.R to centralize styling.


# Read in Site Base data
SiteData <- wb_to_df("data/Sampling Data.xlsx",
  sheet = "AllSamples_Info") %>%
  select(-c(Secchi1_cm, Secchi2_cm, mVpH, mVORP, `DO_%`, `Conductivity_mS/cm`, `AbsConductivity_mS/cm`, Resistance_MΩcm, TDS_ppt, σt, AtmosphericPressure_PSI))

# Build true DateTimes natively before interpolation
SiteData <- SiteData %>%
  mutate(
    Date = as_date(Arrival_DateTime), 
    StartTime = if_else(!is.na(StartTime), ymd_hms(paste(Date, StartTime), tz = "America/Costa_Rica", quiet = TRUE), as.POSIXct(NA, tz = "America/Costa_Rica")),
    EndTime   = if_else(!is.na(EndTime), ymd_hms(paste(Date, EndTime), tz = "America/Costa_Rica", quiet = TRUE), as.POSIXct(NA, tz = "America/Costa_Rica"))
  ) %>%
  # Handle cross-midnight by adding 1 day to EndTime if it's earlier than StartTime
  mutate(
    EndTime = if_else(!is.na(StartTime) & !is.na(EndTime) & hour(EndTime) < hour(StartTime), EndTime + days(1), EndTime)
  ) %>%
  relocate(Date, .before = Arrival_DateTime)

# Impute missing Start/End values, then generate Sample IDs
SiteData <- SiteData %>%
  mutate(
    StartTime = if_else(!is.na(StartTime), StartTime, min(StartTime, na.rm = TRUE)),
    EndTime   = if_else(!is.na(EndTime), EndTime, max(EndTime, na.rm = TRUE)), 
    .by = c(Site, Visit)
  ) %>%
  mutate(SampleID = paste(substr(Site, 1, 3), substr(Visit, 1, 1), "_", substr(Gear, 1, 4), sep = ""), .before = Site) %>% 
  filter(Effort > 0)


## Clearing Excel-averaged values from Seine & Cast Net to prep for later summarizing by Site
SiteData[which(SiteData$Gear != "Centipede Net"), match("SecchiDepth_cm", names(SiteData)):match("Temperature_C", names(SiteData))] <- NA

# Add site-wide descriptors Dominant Substrate (descriptive), Mud Dominance (0,1) & Steepness (0,1,2)
SiteData <- left_join(SiteData, wb_to_df("data/Sampling Data.xlsx",
  sheet = "Site_SubstrateSteepness"))

# Calculate Sunrise/Sunset and add to SiteData
SiteData <- left_join(SiteData, getSunlightTimes(
  data = (SiteData %>%
    select(Date, Latitude, Longitude) %>%
      distinct() %>% 
    rename(date = "Date", lat = "Latitude", lon = "Longitude")),
  tz = "America/Costa_Rica", keep = c("sunrise", "sunset")
), by = join_by(Date == date, Latitude == lat, Longitude == lon), relationship = "many-to-many") %>% 
  rename(Sunrise = sunrise, Sunset = sunset)

# Calculate daylight hours & percent of sample during daylight for all samples
SiteData <- SiteData %>%
  mutate(DaylightHrs = parse_number(as.character(case_when(
    StartTime > Sunrise & EndTime < Sunset ~ difftime(EndTime, StartTime, units = "hours"),
    StartTime < Sunset & EndTime > Sunset ~ difftime(Sunset, StartTime, units = "hours"),
    StartTime < Sunrise & EndTime > Sunrise ~ difftime(EndTime, Sunrise, units = "hours"),
    .default = difftime(EndTime, EndTime, units = "hours")
  ))), .after = EffortUnit) %>% 
  select(-(Sunrise:Sunset)) %>%
  mutate(DaylightPercent = DaylightHrs / as.numeric(difftime(EndTime, StartTime, units = "hours")) * 100, .after = DaylightHrs)

# Read-in Slope/Aspect & Merge with SiteData. Calculated in ArcGIS Pro from "Costa Rica DEM" by j_nelson, sourced from 30 meter NASA SRTM elevation models.
SiteData <- left_join(SiteData, wb_to_df("data/Sites_SlopeAspect.xlsx")) %>%
  relocate(c(SampleID, Site, Visit, Gear, Group)) %>%
  arrange(SampleID)


# Occlusion Image Analysis is completely handled by scripts/Occlusion Processing.R
# Merge Site Data with Occlusion Data
load("data/Occlusion.RData")
OcclDF$Site <- str_replace_all(OcclDF$Site, "PuntaMorales", "Punta Morales")
OcclDF$Group <- as.numeric(OcclDF$Group)

# Join Occlusion data with SiteData
SiteData <- left_join(SiteData, OcclDF[, c(1:5)]) %>%
  relocate(c(SampleID, Site, Visit, Gear, Group))
rm(OcclDF)

# Read in fish data & remove two unidentified small fish. Likely small enough to escape, potentially too small to ID even if preserved. Assumed to be non-novel species already represented in dataset.
IndivData <- wb_to_df("data/Sampling Data.xlsx", sheet = "AllSamples_Fish", cols = "A:I") %>%
  filter(Species != "Anchoa sp. *" & Species != "Juvenil perciformes*") %>%
  select(-c(Collected, Notes))

# Batch weights & counts when n>30 for a species in a single gear
BatchData <- wb_to_df("data/Sampling Data.xlsx", sheet = "AllSamples_Batch") %>%
  select(-TotalWeight_g)

# Combine individual and batch fish data, expanding Batch to 1 row per fish
FishData <- bind_rows(IndivData, as.data.frame(lapply(BatchData, rep, BatchData$Number))) %>%
  select(-Number) %>%
  mutate(SampleID = paste(substr(Site, 1, 3), substr(Visit, 1, 1), "_", substr(Gear, 1, 4), sep = ""), .before = Site)

# Merging Fish data with Site data
MergedData <- full_join(FishData, 
                        SiteData %>%
                          group_by(Site, Visit, Gear, Group) %>% 
                          uncount(
                            ifelse(Gear == 'Centipede Net' | Group == "0" | is.na(Group), 1, max(Effort, na.rm = TRUE)), 
                            .remove = FALSE
                          ) %>%
                          group_by(SampleID) %>%
                          mutate(Group = if_else(Gear == 'Centipede Net', Group, Group + row_number() - 1)) %>% 
                          ungroup(), 
                        by = join_by(SampleID, Site, Visit, Gear, Group)) %>%
  mutate(Count = ifelse(is.na(Species), 0, 1), .after = Species) %>%
  relocate(c(SampleID, Site, Visit, Gear, Group)) %>%
  relocate(Notes, .after = last_col())

rm(IndivData, FishData, BatchData)

## Create list of Species augmented with taxonomic & life history data from FishBase 
load("data/SpeciesList.RData")

# SpeciesList <- MergedData %>%
#   filter(!is.na(Species)) %>%
#   summarize(Abundance = sum(Count), .by = Species) %>%
#   mutate(SpeciesRank = rank(-Abundance, ties.method = "last"), .before = Species) %>%
#   left_join((fb_tbl("species") %>%
#   left_join(fb_tbl("families") %>%
#               select(FamCode, Family), by = "FamCode") %>%
#   mutate(Species = paste(Genus, Species)) %>%
#   filter(Species %in% MergedData$Species) %>%
#   collect()), by = "Species") %>%
#   relocate(SpeciesRank, Species, Abundance, Genus, Family, FBname) %>%
#   mutate(FBname = case_when(
#     Species == "Bathygobius andrei" ~ "Estuarine frillfin",
#     Species == "Diapterus brevirostris" ~ "Peruvian mojarra",
#     Species == "Gerres simillimus" ~ "Yellowfin mojarra",
#     Species == "Sphoeroides rosenblatti" ~ "Oval puffer",
#     .default = FBname))
# 
# # Check for missing data
# filter(SpeciesList, is.na(FBname) | is.na(Genus))

# save(SpeciesList, file = "data/SpeciesList.RData")

FamilyList <- SpeciesList %>%
  select(Family, Abundance) %>%
  summarize(Abundance = sum(Abundance), .by = Family) %>%
  mutate(FamilyRank = rank(-Abundance, ties.method = "first"), .before = Family) %>%
  arrange(FamilyRank)


# Calculate Sample averages + standard deviations for Occlusion, Slope, Aspect, and Site+Visit averages + standard deviations for all other metrics
# Fill values from centipede nets to all gears
SampleData <- SiteData %>%
  select(-c(Group)) %>%
  filter(Effort > 0) %>%
  mutate(Occlusion = mean(OcclusionGroupAvg, na.rm = T), Occlusion_SD = sd(OcclusionGroupAvg, na.rm = T), .by = Site) %>%
  select(-c(OcclusionGroupAvg, OcclusionGroupSTDev)) %>%
  reframe(across(where(is.character), first),
    Effort = if_else(Gear == "Centipede Net", sum(Effort), max(Effort)),
    StartTime = min(StartTime),
    EndTime = max(EndTime),
    Latitude = median(Latitude),
    Longitude = median(Longitude),
    Occlusion = first(Occlusion),
    Occlusion_SD = first(Occlusion_SD),
    MudDominant = first(MudDominant),
    Steepness = first(Steepness),
    across(
      where(is.numeric) & !c(Effort, Latitude, Longitude, Slope, Aspect, Occlusion, Occlusion_SD, MudDominant, Steepness),
      list(Avg = mean, SD = sd)
    ),
    .by = c(Site, Visit, Gear)
  ) %>%
  select(-c(pH_SD, Temperature_C_SD, DO_PPM_SD, DaylightHrs_SD, DaylightPercent_SD)) %>%
  rename_with(~ str_remove_all(.x, "(_Avg|_C|_cm|_PPM)")) %>% # Tidy up names
  relocate(SampleID, Site, Visit, Gear, Season, Effort, EffortUnit) %>%
  group_by(Site, Visit) %>%
  fill(SecchiDepth:Temperature, .direction = "downup") %>% # Fill Cent values across other gears at same visit
  ungroup() %>%
  distinct() %>%
  arrange(SampleID)

# Create Community matrix (Filtering invalid centipede net sample tiv2)
Sample_Community_Matrix <- MergedData %>%
  filter(Effort > 0 & SampleID != 'Tiv2_Cent') %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = Species,
    values_from = Count,
    values_fn = length,
    values_fill = 0
  ) %>%
  column_to_rownames(var = "SampleID") %>%
  t() %>%
  as.data.frame() %>%
  filter(rownames(.) != "NA")

# # Add Shannon & Simpson diversity to Sample data (Don't forget to set to NA when Abundance = 0 before using these data!)
SampleData <- left_join(SampleData, tibble(SampleID = colnames(Sample_Community_Matrix), Shannon = exp(diversity(t(Sample_Community_Matrix), index = "shannon")), Simpson = diversity(t(Sample_Community_Matrix), index = "invsimpson")), by = join_by(SampleID))

# Average Effort by Gear
summarise(SampleData, MeanEffort = mean(Effort, na.rm = TRUE), Unit = first(EffortUnit), .by = Gear) %>%
  mutate(GroupEffort = case_when(Gear == "Centipede Net" ~ MeanEffort / 3))

# Add Fish Sample diversity()summary Values & Standardized CPUE (Standard Sample = 10 Cast Net Throws, 4.64 hr soak of 1 centipede net group, 1.5 seine hauls)
SampleData <- SampleData %>%
  left_join(
    MergedData %>%
      filter(Effort > 0) %>%
      reframe(
        Richness = length(unique(na.omit(Species))),
        Abundance = if (Richness > 0) sum(Count, na.rm = TRUE) else 0,
        Biomass = if (Richness > 0) sum(Weight_g, AverageWeight_g, na.rm = TRUE) else 0,
        AvgLength = if (Richness > 0) mean(StandardLength_mm, na.rm = TRUE) else NaN,
        MinLength = if (Richness > 0) min(StandardLength_mm, na.rm = TRUE) else NaN,
        MaxLength = if (Richness > 0) max(StandardLength_mm, na.rm = TRUE) else NaN,
        AvgWeight = if (Richness > 0) mean(c(AverageWeight_g, Weight_g), na.rm = TRUE) else NaN,
        MinWeight = if (Richness > 0) min(c(AverageWeight_g, Weight_g), na.rm = TRUE) else NaN,
        MaxWeight = if (Richness > 0) max(AverageWeight_g, Weight_g, na.rm = TRUE) else NaN,
        .by = SampleID
      ) %>%
      distinct(),
    by = "SampleID"
  ) %>%
  mutate(
    StandardEffort = case_when(
      Gear == "Cast Net" ~ Effort / 10,
      Gear == "Centipede Net" ~ (Effort / 3) / 4.64,
      Gear == "Seine" ~ Effort / 1.5
    ),
    CPUE = Abundance / StandardEffort,
    logCPUE = log(CPUE + 1),
    BPUE = Biomass / StandardEffort,
    Shannon = if_else(Richness == 0, 0, Shannon),
    Simpson = if_else(Richness == 0, 0, Simpson)
  ) %>%
  mutate(across(everything(), ~ replace(.x, is.infinite(.x) | is.nan(.x), NA))) %>%
  relocate(SampleID, Abundance, CPUE, logCPUE, Richness, Shannon, Simpson, BPUE, Effort, EffortUnit, StandardEffort) %>%
  arrange(SampleID)

# Species abundances by gear
GearSpeciesGrid <- MergedData %>%
  filter(!is.na(Species)) %>%
  select(Species, Gear) %>%
  mutate("Cast Net" = sum(if_else(Gear == "Cast Net", 1, 0)), "Centipede Net" = sum(if_else(Gear == "Centipede Net", 1, 0)), "Seine" = sum(if_else(Gear == "Seine", 1, 0)), .by = Species) %>%
  select(-Gear) %>%
  distinct() %>%
  arrange(Species)

# Checkbox species by gear
GearSpeciesCheck <- MergedData %>%
  filter(!is.na(Species)) %>%
  select(Species, Gear) %>%
  mutate("Cast Net" = if_else(sum(if_else(Gear == "Cast Net", 1, 0)) > 0, "✓", ""), "Centipede Net" = if_else(sum(if_else(Gear == "Centipede Net", 1, 0)) > 0, "✓", ""), "Seine" = if_else(sum(if_else(Gear == "Seine", 1, 0)) > 0, "✓", ""), .by = Species) %>%
  select(-Gear) %>%
  arrange(Species) %>%
  distinct()

# --- NEW: Identify complete site visits BEFORE building matrices ---
# This ensures we only pool data from visits where all required gears were present.
# We use SiteData (which has all valid samples, Effort > 0), exclude invalid Tiv2_Cent
valid_samples <- SiteData %>%
  filter(SampleID != 'Tiv2_Cent') %>%
  select(SampleID, Gear) %>%
  mutate(SiteVisit = substr(SampleID, 1, 4)) %>%
  distinct(SiteVisit, Gear)

# --- 2. Find which SiteVisits have which valid gears ---
visits_with_cast <- valid_samples %>% filter(Gear == "Cast Net") %>% pull(SiteVisit)
visits_with_cent <- valid_samples %>% filter(Gear == "Centipede Net") %>% pull(SiteVisit) # "Tiv2" is now correctly excluded
visits_with_seine <- valid_samples %>% filter(Gear == "Seine") %>% pull(SiteVisit)

# --- 3. Find the intersections (visits that have ALL required gears) ---
visits_cast_cent <- intersect(visits_with_cast, visits_with_cent)
visits_cast_seine <- intersect(visits_with_cast, visits_with_seine) # This will correctly include "Tiv2"
visits_cent_seine <- intersect(visits_with_cent, visits_with_seine)
visits_all_gears <- intersect(visits_cast_seine, visits_with_cent) # This will correctly exclude "Tiv2"


# --- MODIFIED: Build Incidence_Matrices using the complete visit lists ---
Incidence_Matrices <- list(
  "Cast Net" = MergedData %>%
    filter(Gear == "Cast Net" & Effort != 0) %>%
    pivot_wider(id_cols = Species, names_from = SampleID, values_from = Count, values_fn = length, values_fill = 0) %>%
    mutate(across(where(is.numeric), ~ replace(., . != 0, 1))) %>%
    filter(!is.na(Species)) %>%
    column_to_rownames(var = "Species"),
  
  # This one is correct (filters the single invalid sample)
  "Centipede Net" = MergedData %>%
    filter(Gear == "Centipede Net" & Effort != 0 & SampleID != 'Tiv2_Cent') %>%
    pivot_wider(id_cols = Species, names_from = SampleID, values_from = Count, values_fn = length, values_fill = 0) %>%
    mutate(across(where(is.numeric), ~ replace(., . != 0, 1))) %>%
    filter(!is.na(Species)) %>%
    column_to_rownames(var = "Species"),
  
  "Seine" = MergedData %>%
    filter(Gear == "Seine" & Effort != 0) %>%
    pivot_wider(id_cols = Species, names_from = SampleID, values_from = Count, values_fn = length, values_fill = 0) %>%
    mutate(across(where(is.numeric), ~ replace(., . != 0, 1))) %>%
    filter(!is.na(Species)) %>%
    column_to_rownames(var = "Species"),
  
  # FIX: Filter MergedData to only include site visits that have BOTH gears
  "Cast Net & Centipede Net" = MergedData %>%
    filter(substr(SampleID, 1, 4) %in% visits_cast_cent & (Gear == "Cast Net" | Gear == "Centipede Net")) %>%
    mutate(SiteVisit = substr(SampleID, 1, 4)) %>%
    pivot_wider(id_cols = Species, names_from = SiteVisit, values_from = Count, values_fn = length, values_fill = 0) %>%
    mutate(across(where(is.numeric), ~ replace(., . != 0, 1))) %>%
    filter(!is.na(Species)) %>%
    column_to_rownames(var = "Species"),
  
  # FIX: Filter MergedData to only include site visits that have BOTH gears
  "Cast Net & Seine" = MergedData %>%
    filter(substr(SampleID, 1, 4) %in% visits_cast_seine & (Gear == "Cast Net" | Gear == "Seine")) %>%
    mutate(SiteVisit = substr(SampleID, 1, 4)) %>%
    pivot_wider(id_cols = Species, names_from = SiteVisit, values_from = Count, values_fn = length, values_fill = 0) %>%
    mutate(across(where(is.numeric), ~ replace(., . != 0, 1))) %>%
    filter(!is.na(Species)) %>%
    column_to_rownames(var = "Species"),
  
  # FIX: Filter MergedData to only include site visits that have BOTH gears
  "Centipede Net & Seine" = MergedData %>%
    filter(substr(SampleID, 1, 4) %in% visits_cent_seine & (Gear == "Centipede Net" | Gear == "Seine")) %>%
    mutate(SiteVisit = substr(SampleID, 1, 4)) %>%
    pivot_wider(id_cols = Species, names_from = SiteVisit, values_from = Count, values_fn = length, values_fill = 0) %>%
    mutate(across(where(is.numeric), ~ replace(., . != 0, 1))) %>%
    filter(!is.na(Species)) %>%
    column_to_rownames(var = "Species"),
  
  # FIX: Filter MergedData to only include site visits that have ALL THREE gears
  "All Gears" = MergedData %>%
    filter(substr(SampleID, 1, 4) %in% visits_all_gears) %>%
    mutate(SiteVisit = substr(SampleID, 1, 4)) %>%
    pivot_wider(id_cols = Species, names_from = SiteVisit, values_from = Count, values_fn = length, values_fill = 0) %>%
    mutate(across(where(is.numeric), ~ replace(., . != 0, 1))) %>%
    filter(!is.na(Species)) %>%
    column_to_rownames(var = "Species")
)

# Save core wrangled dataframes for analysis and modeling
save(MergedData, SampleData, SpeciesList, FamilyList, Sample_Community_Matrix, Incidence_Matrices, file = "data/WrangledData.RData")