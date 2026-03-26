
# Libraries and reading/loading data --------------------------------------
library("MuMIn")
library("car")
library("emmeans")
library("scales")
library("parallel")
library("multcomp")
library("multcompView")
library("DHARMa")
library("performance")
library("openxlsx2")
library("cowplot")
library("modelbased")
library("glmmTMB")
library("tidyverse")


load("data/WrangledData.RData")
load("output/models/DredgeListGLMM.RData")

rm(FamilyList, GearSpeciesCheck, GearSpeciesGrid, MergedData, Sample_Community_Matrix, SiteData, SpeciesList, Incidence_Matrices)

# Source the global ggplot2 theme and aesthetic scales
source("scripts/0_PlotTheme.R")
# Data Wrangling ----------------------------------------------------------
#**SCALING EFFORT AGAIN. Switched to GLMMs with Effort x Gear Interaction**
# Removing semi-failed centipede net sample where nets were never fully inundated:
SampleData <- SampleData %>% 
  filter(!(Gear == 'Centipede Net' & Effort < 3))

CorMat <- cor(SampleData %>%
                select(-c(SampleID, logCPUE, Shannon, Richness, Simpson, Effort, EffortUnit, StartTime, EndTime, Season, Site, Gear, DominantSubstrate)) %>%
                na.omit())

# Create dataset with scaled & centered numeric predictors.
dat <- SampleData %>%
  mutate(logEffort = log(Effort)) %>%
  mutate(logEffort = scale(logEffort)[,1], .by = Gear) %>%
  mutate(across(c(Season, Site, Gear), as.factor)) %>%
  mutate(Steepness = factor(Steepness, levels = c(0,1,2), labels = c("Low", "Medium", "High"))) %>%
  mutate(MudDominant = factor(MudDominant, levels = c(0,1), labels = c("False", "True"))) %>%
  mutate(Occlusion = scale(Occlusion)[,1], DaylightPercent = scale(DaylightPercent)[,1], SecchiDepth = scale(SecchiDepth)[,1], Temperature = scale(Temperature)[,1])

# Scaled Data Summary (Means & SDs to report)
ScaleData <- SampleData %>%
  select(Occlusion:Temperature) %>%
  summarise(across(everything(),
                   .fns = list(Mean = ~ mean(., na.rm = TRUE), SD = ~ sd(., na.rm = TRUE))
  )) %>%
  pivot_longer(everything(),
               names_to = c("Predictor", ".value"),
               names_pattern = "(.*)_(.*)"
  ) %>%
  group_by(Predictor) %>%
  summarise(Mean = Mean[1], SD = SD[1]) %>%
  split(.$Predictor) %>%
  purrr::map(~ as.list(.x)[-1])

ScaleDataEffort <- SampleData %>%
  mutate(logEffort = log(Effort)) %>%
  group_by(Gear) %>%
  summarize(
    Mean = mean(logEffort, na.rm = TRUE),
    SD = sd(logEffort, na.rm = TRUE)
  )

# Effort Mean & SD by Gear
SummaryData_Effort <- SampleData %>%
  select(Effort, Gear) %>%
  group_by(Gear) %>%
    summarize(across(where(is.numeric),
                   .fns = list(Mean = ~ mean(., na.rm = TRUE), SD = ~ sd(., na.rm = TRUE))
  )) %>%
  ungroup() %>%
  pivot_longer(-Gear,
               names_to = c("Predictor", ".value"),
               names_pattern = "(.*)_(.*)"
  ) %>%
  group_by(Predictor, Gear) %>%
  summarise(Mean = Mean[1], SD = SD[1]) %>%
  pivot_wider(names_from = Gear, values_from = c(Mean, SD)) %>%
  split(.$Predictor) %>%
  purrr::map(~ as.list(.x)[-1])

# NA-Free Dataset for Modelling
geardat <- dat %>%
  select(Effort, logEffort, Abundance, Richness, Shannon, Simpson, MudDominant, Steepness, Season, Gear, DaylightPercent, SecchiDepth, Occlusion, Temperature, Site) %>%
  na.omit()


# Finding appropriate distributions for each response variable -------------
# (These models take a minute to run and are saved in RData files. 
# They are preserved here in an inactive block for reproducibility and reviewer access.)
# run_distribution_tests <- FALSE
# 
# if (run_distribution_tests) {
  # Abundance: Selecting best distribution for Abundance ~ Gear
  GearAbund_Gaussian <- glmmTMB(Abundance ~ Gear*logEffort + (1 | Site),
                            data = geardat,
                            family = gaussian)
  
  GearAbund_Poisson <- glmmTMB(Abundance ~ Gear*logEffort + (1 | Site),
                           data = geardat,
                           family = poisson)
  
  GearAbund_NBinom <- glmmTMB(Abundance ~ Gear*logEffort + (1 | Site),
                          data = geardat,
                          family = nbinom1())
  
  GearAbund_NBinom2 <- glmmTMB(Abundance ~ Gear*logEffort + (1 | Site),
                           data = geardat,
                           family = nbinom2())
  
  
  # Compare the models using AICc
  model.sel(GearAbund_Gaussian, GearAbund_Poisson, GearAbund_NBinom, GearAbund_NBinom2)
  # nbinom2 is the best fitting distribution by a long shot
  
  # Diagnostic plots:
  simulateResiduals(fittedModel = GearAbund_Poisson, plot = TRUE)
  simulateResiduals(fittedModel = GearAbund_NBinom, plot = TRUE)
  simulateResiduals(fittedModel = GearAbund_NBinom2, plot = TRUE) # QQ plot looks pretty good. Predicted vs Residual shows consistent overestimation for the highest values. Could be because of nonlinear impacts of other predictors/interactions that aren't present in this base model?
  testOutliers(GearAbund_Poisson)
  testDispersion(GearAbund_NBinom)
  testDispersion(GearAbund_NBinom2)
  Anova(GearAbund_NBinom2)
  cld(emmeans(GearAbund_NBinom2, specs = ~ Gear, type = "response"), Letters = letters)
  #**# There is a significant Gear:logEffort interaction**
  #* The base model for Abundance should be glmmTMB(Abundance ~ Gear*logEffort + (1 | Site), data = geardat, family = nbinom2())
  
  
  # Richness: Selecting best distribution for Richness ~ Gear
  GearRich_Gaussian <- glmmTMB(Richness ~ Gear*logEffort + (1 | Site),
                                data = geardat,
                                family = gaussian)
  
  GearRich_Poisson <- glmmTMB(Richness ~ Gear*logEffort + (1 | Site),
                               data = geardat,
                               family = poisson)

  GearRich_NBinom <- glmmTMB(Richness ~ Gear*logEffort + (1 | Site),
                              data = geardat,
                              family = nbinom1())
  
  GearRich_NBinom2 <- glmmTMB(Richness ~ Gear*logEffort + (1 | Site),
                               data = geardat,
                               family = nbinom2()) #Did not converge
  
  
  # Compare the models using AICc
  model.sel(GearRich_Gaussian, GearRich_Poisson, GearRich_NBinom, GearRich_NBinom2)
  #** Poisson has lowest AICc, followed by NBinom1 (delta = 2.45)*
  
  # Diagnostic plots:
  simulateResiduals(fittedModel = GearRich_Gaussian, plot = TRUE)
  simulateResiduals(fittedModel = GearRich_Poisson, plot = TRUE) 
  simulateResiduals(fittedModel = GearRich_NBinom, plot = TRUE)
  simulateResiduals(fittedModel = GearRich_NBinom2, plot = TRUE)
  
  testOutliers(GearRich_Poisson) 
  testDispersion(GearRich_Poisson)
  Anova(GearRich_Poisson) 
  summary(GearRich_Poisson)
  #* The base model for Richness should be glmmTMB(Richness ~ Gear*logEffort + (1 | Site), data = geardat, family = poisson)
  
  
  
  # Shannon: Selecting best distribution for Shannon ~ Gear
  GearShannon_Gaussian <- glmmTMB(Shannon ~ Gear + (1 | Site),
                               data = geardat,
                               family = gaussian)
  
  GearShannon_Tweedie <- glmmTMB(Shannon ~ Gear + (1 | Site),
                              data = geardat,
                              family = tweedie(link = "log"))
  
  
  model.sel(GearShannon_Gaussian, GearShannon_Tweedie)
  
  simulateResiduals(fittedModel = GearShannon_Gaussian, plot = TRUE)
  simulateResiduals(fittedModel = GearShannon_Tweedie, plot = TRUE)
  testOutliers(GearShannon_Gaussian)
  testOutliers(GearShannon_Tweedie)
  testDispersion(GearShannon_Gaussian)
  testDispersion(GearShannon_Tweedie)
  #* The base model for Shannon should be glmmTMB(Shannon ~ Gear + (1 | Site), data = geardat, family = tweedie(link = "log"))
  
  # Simpson: Selecting the best distribution for Simpson ~ Gear
  GearSimpson_Gaussian <- glmmTMB(Simpson ~ Gear + (1 | Site),
                                  data = geardat,
                                  family = gaussian)
  
  GearSimpson_Tweedie <- glmmTMB(Simpson ~ Gear + (1 | Site),
                                 data = geardat,
                                 family = tweedie())
  
  model.sel(GearSimpson_Gaussian, GearSimpson_Tweedie)
  simulateResiduals(fittedModel = GearSimpson_Gaussian, plot = TRUE)
  simulateResiduals(fittedModel = GearSimpson_Tweedie, plot = TRUE)
  testOutliers(GearSimpson_Gaussian)
  testOutliers(GearSimpson_Tweedie)
  testDispersion(GearSimpson_Gaussian)
  testDispersion(GearSimpson_Tweedie)
  #* The base model for Simpson should be glmmTMB(Simpson ~ Gear + (1 | Site), data = geardat, family = tweedie(link = "log"))
# }


# Fit the Simple GLMMs (Gear as the only predictor) ---
# Abundance Model
fit_abund <- glmmTMB(Abundance ~ Gear + (1 | Site), data = geardat, family = nbinom2)
# Richness Model
fit_rich <- glmmTMB(Richness ~ Gear + (1 | Site), data = geardat, family = poisson)
# Shannon Diversity Model
fit_div <- glmmTMB(Shannon ~ Gear + (1 | Site), data = geardat, family = tweedie)
# Simpson Diversity Model
fit_sim <- glmmTMB(Simpson ~ Gear + (1 | Site), data = geardat, family = tweedie)
# 
#Gear differences in response variables, quick checks ----------------------------------
  # Does Effort differ among gears?
  GearEffort <- lm(Effort ~ factor(Gear), data = geardat)
anova(GearEffort)
summary(GearEffort)
# Yes. Effort differs significantly among gears, so a simple offset would be inappropriate. A Gear:Effort interaction will be included.

# Does Gear significantly impact the response variables?
Anova(fit_abund) # Yes
Anova(fit_rich) # Yes
Anova(fit_div) # Yes
summary(fit_sim) # Yes
# Also, because these indices are theoretically (and mathematically) descriptions of the balance of species, expecting them to scale with increased effort may be unfounded.(Gotelli et al., 2024)

# Tukey plots -------------------------------------------------------------
# https://schmidtpaul.github.io/dsfair_quarto/ch/summaryarticles/compactletterdisplay.html


# --- 2. Perform Post-Hoc Tests (emmeans + cld) ---
cld_abund <- cld(emmeans(fit_abund, specs = ~ Gear, type = "response"), Letters = letters)
cld_rich  <- cld(emmeans(fit_rich, specs = ~ Gear, type = "response"), Letters = letters)
cld_div   <- cld(emmeans(fit_div, specs = ~ Gear, type = "response"), Letters = letters)
cld_sim   <- cld(emmeans(fit_sim, specs = ~ Gear, type = "response"), Letters = letters)

# Function to get pairwise p-values
get_pairwise_pvals <- function(model, name) {
  print(paste("---", name, "Pairwise Comparisons ---"))
  pairs_out <- pairs(emmeans(model, specs = ~ Gear, type = "response"))
  print(pairs_out)
}

# Run for each model
get_pairwise_pvals(fit_abund, "Abundance")
get_pairwise_pvals(fit_rich, "Richness")
get_pairwise_pvals(fit_div, "Shannon")
get_pairwise_pvals(fit_sim, "Simpson")

# --- 2. Create the Individual Plots with Nudged Layout ---

# Abundance Plot (with pseudo-log y-axis)
p_abund <- ggplot(geardat, aes(x = as.numeric(as.factor(Gear)), y = Abundance)) +
  geom_jitter(aes(x = as.numeric(as.factor(Gear)) - 0.2), shape = 16, color = "gray50", alpha = 0.75, width = 0.05, height = 0) +
  geom_boxplot(aes(group = Gear, fill = Gear), alpha = 0.75, width = 0.1, outlier.shape = NA) +
  geom_point(data = cld_abund, aes(x = as.numeric(as.factor(Gear)), y = response), size = 1, color = "red", position = position_nudge(x = 0.2)) +
  geom_errorbar(data = cld_abund, aes(x = as.numeric(as.factor(Gear)), y = NULL, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1, color = "red", position = position_nudge(x = 0.2)) +
  geom_text(data = cld_abund, aes(x = as.numeric(as.factor(Gear)), y = asymp.UCL, label = str_trim(.group)), 
            # size = 4, 
            position = position_nudge(x = 0.3), hjust = 0) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10), breaks = c(0, 2, 10, 50, 150, 500), labels = label_number(accuracy = 1)) +
  scale_x_continuous(breaks = 1:3, labels = sort(unique(as.character(geardat$Gear))), expand = expansion(mult = c(0.05, 0.10))) +
  scale_fill_manual(values = GearColors()) +
  labs(y = "Abundance (log scale)", x = NULL, title = "Abundance", subtitle = "") +
  theme(legend.position = "none", plot.title.position = "plot", plot.title = element_text(vjust = -1))

# Richness Plot
p_rich <- ggplot(geardat, aes(x = as.numeric(as.factor(Gear)), y = Richness)) +
  geom_jitter(aes(x = as.numeric(as.factor(Gear)) - 0.2), shape = 16, color = "gray50", alpha = 0.75, width = 0.05, height = 0) +
  geom_boxplot(aes(group = Gear, fill = Gear), alpha = 0.75, width = 0.1, outlier.shape = NA) +
  geom_point(data = cld_rich, aes(x = as.numeric(as.factor(Gear)), y = rate), size = 1, color = "red", position = position_nudge(x = 0.2)) +
  geom_errorbar(data = cld_rich, aes(x = as.numeric(as.factor(Gear)), y = NULL, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1, color = "red", position = position_nudge(x = 0.2)) +
  geom_text(data = cld_rich, aes(x = as.numeric(as.factor(Gear)), y = asymp.UCL, label = str_trim(.group)), 
            # size = 4, 
            position = position_nudge(x = 0.3), hjust = 0) +
  scale_y_continuous(breaks = c(0,3,6,9,12)) +
  scale_x_continuous(breaks = 1:3, labels = sort(unique(as.character(geardat$Gear))), expand = expansion(mult = c(0.05, 0.10))) +
  scale_fill_manual(values = GearColors()) +
  labs(y = "Species Richness", x = NULL, title ="Species Richness", subtitle = "") +
  theme(legend.position = "none", plot.title.position = "plot", plot.title = element_text(vjust = -1))

# Shannon Diversity Plot
p_div <- ggplot(geardat, aes(x = as.numeric(as.factor(Gear)), y = Shannon)) +
  geom_jitter(aes(x = as.numeric(as.factor(Gear)) - 0.2), shape = 16, color = "gray50", alpha = 0.75, width = 0.05, height = 0) +
  geom_boxplot(aes(group = Gear, fill = Gear), alpha = 0.75, width = 0.1, outlier.shape = NA) +
  geom_point(data = cld_div, aes(x = as.numeric(as.factor(Gear)), y = response), size = 1, color = "red", position = position_nudge(x = 0.2)) +
  geom_errorbar(data = cld_div, aes(x = as.numeric(as.factor(Gear)), y = NULL, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1, color = "red", position = position_nudge(x = 0.2)) +
  geom_text(data = cld_div, aes(x = as.numeric(as.factor(Gear)), y = asymp.UCL, label = str_trim(.group)), 
            # size = 4, 
            position = position_nudge(x = 0.3), hjust = 0) +
  scale_x_continuous(breaks = 1:3, labels = sort(unique(as.character(geardat$Gear))), expand = expansion(mult = c(0.05, 0.10))) +
  scale_fill_manual(values = GearColors()) +
  labs(y = "Shannon Diversity", x = NULL, title = "Shannon Diversity", subtitle = "") +
  theme(legend.position = "none", plot.title.position = "plot", plot.title = element_text(vjust = -1))

# Simpson Diversity Plot
p_sim <- ggplot(geardat, aes(x = as.numeric(as.factor(Gear)), y = Simpson)) +
  geom_jitter(aes(x = as.numeric(as.factor(Gear)) - 0.2), shape = 16, color = "gray50", alpha = 0.75, width = 0.05, height = 0) +
  geom_boxplot(aes(group = Gear, fill = Gear), alpha = 0.75, width = 0.1, outlier.shape = NA) +
  geom_point(data = cld_sim, aes(x = as.numeric(as.factor(Gear)), y = response), size = 1, color = "red", position = position_nudge(x = 0.2)) +
  geom_errorbar(data = cld_sim, aes(x = as.numeric(as.factor(Gear)), y = NULL, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1, color = "red", position = position_nudge(x = 0.2)) +
  geom_text(data = cld_sim, aes(x = as.numeric(as.factor(Gear)), y = asymp.UCL, label = str_trim(.group)), 
            # size = 4, 
            position = position_nudge(x = 0.3), hjust = 0) +
  scale_y_continuous(breaks = c(1,3,5,7)) +
  scale_x_continuous(breaks = 1:3, labels = sort(unique(as.character(geardat$Gear))), expand = expansion(mult = c(0.05, 0.10))) +
  scale_fill_manual(values = GearColors()) +
  labs(y = "Simpson Diversity", x = NULL, title = "Simpson Diversity", subtitle = "") +
  theme(legend.position = "none", plot.title.position = "plot", plot.title = element_text(vjust = -1))

# --- 3. Assemble the Final 2x2 Grid ---
final_tukey_plot <- plot_grid(
  p_abund, p_rich, p_div, p_sim,
  ncol = 2,
  align = "hv" 
  )  +
  plot_theme

final_tukey_plot
# ggsave(plot = final_tukey_plot, "output/plots/GearTukeyGrid.eps", device = cairo_ps,  width = 5.62, height = 5.62, units = "in")
ggsave(plot = final_tukey_plot, "output/plots/GearTukeyGrid.png",  width = 5.62, height = 5.62, units = "in", dpi = 300)


# Model Selection------------------------------------------------
# Full models with Interactions, updated for GLMMs

# Global model for Abundance
AbundIntFit <- glmmTMB(
  Abundance ~ Gear * logEffort + Gear * (Occlusion + SecchiDepth + DaylightPercent) +
    Steepness + MudDominant + Season + Temperature + (1 | Site),
  data = geardat,
  family = nbinom2, na.action = "na.fail"
)

# Global model for Richness
RichIntFit <- glmmTMB(
  Richness ~ Gear * logEffort + Gear * (Occlusion + SecchiDepth + DaylightPercent) +
    Steepness + MudDominant + Season + Temperature + (1 | Site),
  data = geardat,
  family = poisson(link = "log"), na.action = "na.fail"
)

# Global model for Shannon Diversity (Semi-continuous Data with Zeros)
DivIntFit <- glmmTMB(
  Shannon ~ Gear * logEffort + Gear * (Occlusion + SecchiDepth + DaylightPercent) +
    Steepness + MudDominant + Season + Temperature + (1 | Site),
  data = geardat,
  family = tweedie(link = "log"), na.action = "na.fail"
)

# Global model for Simpson Diversity (Semi-continuous Data with Zeros)
SimIntFit <- glmmTMB(
  Simpson ~ Gear * logEffort + Gear * (Occlusion + SecchiDepth + DaylightPercent) +
    Steepness + MudDominant + Season + Temperature + (1 | Site),
  data = geardat,
  family = tweedie(), na.action = "na.fail"
)

# Check VIFs on global models
# --- 1. Variance Inflation Factors (Collinearity) ---
# We use performance::check_collinearity() because it natively handles the complex covariance structure of glmmTMB models, unlike car::vif().
check_collinearity(AbundIntFit)
check_collinearity(RichIntFit)
check_collinearity(DivIntFit)
check_collinearity(SimIntFit)

# --- 2. Residual Diagnostics (DHARMa) ---
# Run on global models to validate distributional assumptions before dredging.
# Saves plots to output/plots/Diagnostics/ for documentation.
if(!dir.exists("output/plots/Diagnostics")) dir.create("output/plots/Diagnostics")

run_dharma <- function(model, name) {
  print(paste("--- DHARMa Diagnostics:", name, "---"))
  sim_res <- simulateResiduals(fittedModel = model, plot = FALSE)
  print(testResiduals(sim_res))
  print(testDispersion(sim_res))
  png(filename = paste0("output/plots/Diagnostics/DHARMa_", name, ".png"), 
      width = 800, height = 500)
  plot(sim_res)
  dev.off()
}

run_dharma(AbundIntFit, "Abundance_Global")
run_dharma(RichIntFit,  "Richness_Global")
run_dharma(DivIntFit,   "Shannon_Global")
run_dharma(SimIntFit,   "Simpson_Global")

# --- Create Null Models for R-squared Calculation ---

# Null model for Shannon Diversity
DivNullFit <- glmmTMB(
  Shannon ~ 1 + (1 | Site),
  data = geardat,
  family = tweedie(link = "log"),
  na.action = "na.fail"
)

# Null model for Simpson Diversity
SimNullFit <- glmmTMB(
  Simpson ~ 1 + (1 | Site),
  data = geardat,
  family = tweedie(link = "log"),
  na.action = "na.fail"
)

# High VIFs for interaction terms & involved main effects. Check "global" models without interactions to determine if multicollinearity is present:
check_additive_vif <- function(data, family, response_var) {
  f <- as.formula(paste(response_var, "~ Gear + logEffort + MudDominant + Steepness + Season + DaylightPercent + SecchiDepth + Occlusion + Temperature + (1 | Site)"))
    fit <- glmmTMB(f, data = data, family = family, na.action = "na.fail")
  print(paste("--- VIF for", response_var, "(Additive Model) ---"))
  print(performance::check_collinearity(fit))
}

# Run the check for each response
check_additive_vif(geardat, nbinom2, "Abundance")
check_additive_vif(geardat, poisson, "Richness")
check_additive_vif(geardat, tweedie(link="log"), "Shannon")
check_additive_vif(geardat, tweedie(link="log"), "Simpson")

# # Dredging------------------------------------------------------------
# #** These will take some time (~30-90 mins. depending on computer specs) **
# # Create the cluster object (leave 1 core free for OS stability)
# clu <- makeCluster(detectCores() - 1)
# 
# # Export your data frame to each worker node
# clusterExport(clu, varlist = c("geardat", "DivNullFit", "SimNullFit"))
# 
# # Load the required modeling package on all workers
# clusterEvalQ(clu, library(glmmTMB))
# clusterEvalQ(clu, library(MuMIn))
# 
# # --- 2. Run All dredge() Calls in Parallel ---
# 
# # A. Dredge the Abundance model
# AbundIntTest <- dredge(
#   AbundIntFit,
#   fixed = ~cond(Gear) + cond(logEffort) + cond(Gear:logEffort),
#   extra = "r.squaredGLMM",
#   cluster = clu
# )
# 
# # B. Dredge the Richness model
# RichIntTest <- dredge(
#   RichIntFit,
#   fixed = ~cond(Gear) + cond(logEffort) + cond(Gear:logEffort),
#   extra = "r.squaredGLMM",
#   cluster = clu
# )
# 
# # C. Dredge the Shannon Diversity model
# DivIntTest <- dredge(
#   DivIntFit,
#   fixed = ~cond(Gear),
#   extra = list(R2.LR = function(x) r.squaredLR(x, null = DivNullFit)),
#   cluster = clu
# )
# 
# # D. Dredge the Simpson Diversity model
# SimIntTestAll <- dredge(
#   SimIntFit,
#   fixed = ~cond(Gear),
#   extra = list(R2.LR = function(x) r.squaredLR(x, null = SimNullFit)),
#   cluster = clu
# )
# # Filter out any models that failed to converge
# SimIntTest <- subset(SimIntTestAll, !is.na(AICc), recalc.delta = TRUE)
# 
# 
# # --- 3. Stop the Cluster ---
# 
# # IMPORTANT: Always stop the cluster to release system resources
# stopCluster(clu)
# 
# # #Play some tunes when done running.
# # shell.exec("https://youtu.be/yebNIHKAC4A?si=J7uIqIlqzWXGIdKg")
# 
# 
# # AIC Table Building ------------------------------------------------------
# # A new, more robust function to create a clean AICc table from a dredge object
# BuildTable <- function(dredge_object, response_name) {
#   dredge_subset <- subset(dredge_object, !nested(.))
#   model_calls <- attr(dredge_subset, "model.calls")
#   
#   predictor_strings <- purrr::map_chr(model_calls, function(call) {
#     terms <- attr(terms(call$formula), "term.labels")
#     terms <- grep("\\|", terms, value = TRUE, invert = TRUE)
#     sort_key <- dplyr::case_when(
#       terms == "logEffort" ~ 1,
#       stringr::str_detect(terms, "logEffort") & stringr::str_detect(terms, ":") ~ 2,
#       stringr::str_detect(terms, "Gear") ~ 3,
#       TRUE ~ 4
#     )
#     terms <- terms[order(sort_key, terms)]
#     if (length(terms) == 0) {
#       "(Intercept only)"
#     } else {
#       paste(terms, collapse = " + ")
#     }
#   })
#   
#   final_table_raw <- as.data.frame(dredge_subset) %>%
#     rownames_to_column(var = "Model") %>%
#     mutate(
#       Predictors = predictor_strings,
#       Response = response_name,
#       weight_tot = cumsum(weight)
#     )
#   
#   # --- New, robust R-squared logic ---
#   has_glmm_r2 <- "r.squaredGLMM1" %in% colnames(final_table_raw)
#   has_lr_r2 <- "R2.LR" %in% colnames(final_table_raw)
#   
#   if (has_glmm_r2) {
#     final_table_processed <- final_table_raw %>%
#       mutate(
#         R2_Marginal = r.squaredGLMM1,    # R2m (Fixed effects)
#         R2_Conditional = r.squaredGLMM2 # R2c (Fixed + Random)
#       )
#   } else if (has_lr_r2) {
#     final_table_processed <- final_table_raw %>%
#       mutate(
#         R2_Marginal = R2.LR # Pseudo-R2 (Fixed effects)
#       )
#   } else {
#     final_table_processed <- final_table_raw # No R2 columns found
#   }
#   
#   # Select the final, clean set of columns
#   final_table_clean <- final_table_processed %>%
#     select(
#       Response, Model, Predictors, df, logLik, AICc, delta, weight, 
#       any_of(c("R2_Marginal", "R2_Conditional")), # Select them if they exist
#       weight_tot
#     )
#   
#   return(final_table_clean)
# }
# 
# # Build AIC Tables
# AbundIntTable <- BuildTable(AbundIntTest, "Abundance")
# RichIntTable <- BuildTable(RichIntTest, "Richness")
# DivIntTable <- BuildTable(DivIntTest, "Shannon")
# SimIntTable <- BuildTable(SimIntTest, "Simpson")
# 
# # # Build Coefficient tables for top model or averaged model if multiple delta <= 2.0
# nrow(AbundIntTable %>% filter(delta <= 2)) # Single top model
# nrow(RichIntTable %>% filter(delta <= 2)) # 7 Similarly parsimonious models
# nrow(DivIntTable %>% filter(delta <= 2)) # 4 Similarly parsimonious models
# nrow(SimIntTable %>% filter(delta <= 2)) # Single top model
# 
# TopAvgCoefs <- bind_rows(
#   # Abundance Top Model Coefficients (K=68)
#   coefTable(AbundIntTest[1,]) %>%
#     as.data.frame() %>%
#     rownames_to_column(var = "Predictor") %>%
#     mutate(Response = "Abundance", .before = 1) %>%
#     rename_with(
#       ~ .x %>%
#         str_remove("X141.") %>%
#         str_replace("\\..", ".")) %>%
#     mutate(Predictor = Predictor %>%
#              str_remove_all("cond|disp") %>%
#              str_remove_all("[()]") %>%
#              str_replace_all(":", "*") %>%
#              str_replace_all("\\s", "") %>%
#              str_replace("^Int$", "Intercept")) %>%
#     select(-df) %>%
#     mutate(Model = "Top",Response = "Abundance", .before = 1)  %>%
#     rename("Std. Error" = Std.Error),
#   
#   # Richness Avg Model Coefficients
#   coefTable(model.avg(RichIntTest, subset = delta <= 2 & !nested(.))) %>%
#     as.data.frame() %>%
#     rownames_to_column(var = "Predictor") %>%
#     mutate(Model = "Averaged", Response = "Richness", .before = 1) %>%
#     mutate(Predictor = Predictor %>%
#              str_remove_all("cond|disp") %>%
#              str_remove_all("[()]") %>%
#              str_replace_all(":", "*") %>%
#              str_replace_all("\\s", "") %>%
#              str_replace("^Int$", "Intercept")) %>%
#     select(-df),
# 
#   # Shannon Diversity Avg Model Coefficients
#   coefTable(model.avg(DivIntTest, subset = delta <= 2 & !nested(.))) %>%
#     as.data.frame() %>%
#     rownames_to_column(var = "Predictor") %>%
#     mutate(Model = "Averaged", Response = "Shannon", .before = 1) %>%
#     mutate(Predictor = Predictor %>%
#              str_remove_all("cond|disp") %>%
#              str_remove_all("[()]") %>%
#              str_replace_all(":", "*") %>%
#              str_replace_all("\\s", "") %>%
#              str_replace("^Int$", "Intercept")) %>%
#     select(-df),
# 
#   # Simpson Diversity Top Model Coefficients
#   coefTable(SimIntTest[1,]) %>%
#     as.data.frame() %>%
#     rownames_to_column(var = "Predictor") %>%
#     mutate(Response = "Simpson", .before = 1) %>%
#     rename_with(
#       ~ .x %>%
#         str_remove("X1.") %>%
#         str_replace("\\..", ".")) %>%
#     mutate(Predictor = Predictor %>%
#              str_remove_all("cond|disp") %>%
#              str_remove_all("[()]") %>%
#              str_replace_all(":", "*") %>%
#              str_replace_all("\\s", "") %>%
#              str_replace("^Int$", "Intercept")) %>%
#     select(-df) %>%
#     mutate(Model = "Top",Response = "Simpson", .before = 1)  %>%
#     rename("Std. Error" = Std.Error)
#   
# ) %>%
#   rename(SE = "Std. Error") %>%
#   mutate(LCL = Estimate - (1.96 * SE),
#          UCL = Estimate + (1.96 * SE),
#          Sig = if_else((LCL > 0 & UCL > 0) | (LCL < 0 & UCL < 0), "*", "")) %>% 
#   # --- Final Polishing Mutate (as in your working code) ---
#   mutate(
#     Predictor = Predictor %>%
#       str_replace_all("\\*", " * ") %>%
#       str_replace_all("logEffort", "Effort") %>% 
#       str_replace_all("GearCentipedeNet", "Centipede Net") %>%
#       str_replace_all("GearSeine", "Seine") %>%
#       str_replace_all("DaylightPercent", "Daylight Percent") %>%
#       str_replace_all("MudDominantTrue", "Mud Dominant") %>%
#       str_replace_all("SeasonRainy", "Rainy Season") %>%
#       str_replace_all("SecchiDepth", "Secchi Depth") %>%
#       str_replace_all("SteepnessMedium", "Steepness (Medium)") %>%
#       str_replace_all("SteepnessHigh", "Steepness (High)")
#   ) %>%
#   
#   # --- New, Final Sorting Logic ---
#   mutate(
#     sort_key = case_when(
#       Predictor == "Intercept" ~ 1,
#       Predictor == "logEffort" ~ 2,
#       str_detect(Predictor, "\\*") ~ 3, # All interactions
#       TRUE ~ 4 # All other main effects
#     )
#   ) %>%
#   arrange(Response, sort_key, Predictor) %>% # Sorts by Response, then your new rules
#   select(-sort_key) # Remove the helper column


# Saving to Excel
# Set Number Format
# options("openxlsx2.numFmt" = "#,##0.00")
# 
# # Single xlsx file with 1) Stacked AIC tables per gear, 2) Stacked coefficient tables for top models, 3) Individual gear/response tables
# write_xlsx(
#   list(
#     "ModelList" = bind_rows(AbundIntTable, RichIntTable, DivIntTable, SimIntTable),
#     "TopAvgCoefficients" = TopAvgCoefs),
#   file = "output/tables/AICTables_GLMM-ReducedMod.xlsx", as_table = TRUE, tableStyle = "TableStyleLight1", na.strings = ""
# )
# 
# # Creating list of top/averaged models
# # Abundance Top Model
# AbundIntTop <- get.models(AbundIntTest, subset = 1)[[1]]
# # Richness Model Averaging
# RichIntTop <- get.models(RichIntTest, subset = delta <= 2 & !nested(.))
# RichIntAvg <- model.avg(RichIntTop, fit = TRUE)
# # Shannon Diversity Model Averaging
# DivIntTop <- get.models(DivIntTest, subset = delta <= 2 & !nested(.)) %>%
#   purrr::keep(~ !is.na(AIC(.)))
# DivIntAvg <- model.avg(DivIntTop, fit = TRUE)
# # Simpson Diversity Top Model
# SimIntTop <- get.models(SimIntTest, subset = 1)[[1]]
# 
# AvgMods <- list(
#   AbundInt = AbundIntTop,
#   RichInt  = RichIntAvg,
#   DivInt   = DivIntAvg,
#   SimInt   = SimIntTop
# )
# save(geardat, ScaleData, AbundIntFit, RichIntFit, DivIntFit, SimIntFit, AbundIntTest, RichIntTest, DivIntTest, SimIntTest, SimIntTestAll, TopAvgCoefs, AbundIntTable, RichIntTable, DivIntTable, SimIntTable, AvgMods, file = "output/models/DredgeListGLMM.RData")

# Marginal Effects ----------------------------------- 
# Dummy plot to extract shared legends
shared_legend_horiz <- get_legend(
  ggplot(data = geardat, aes(x = Gear, y = Abundance, color = Gear, shape = Gear, linetype = Gear, fill = Gear)) +
    geom_ribbon(aes(ymin = 0, ymax = 0), alpha = 0.2) + geom_line(linewidth = 1) + geom_point(size = 2) +
    scale_color_manual(values = GearColors()) + scale_fill_manual(values = GearColors()) +
    scale_shape_manual(values = GearShapes) + scale_linetype_manual(values = GearLines()) +
    guides(color = guide_legend()) +
    theme(legend.position = "top", legend.title = element_blank(), legend.key.width = rel(1), legend.direction = "horizontal")
)
shared_legend_vert <- get_legend(
  ggplot(data = geardat, aes(x = Gear, y = Abundance, color = Gear, shape = Gear, linetype = Gear, fill = Gear)) +
    geom_ribbon(aes(ymin = 0, ymax = 0), alpha = 0.2) + geom_line(linewidth = 1) + geom_point(size = 2) +
    scale_color_manual(values = GearColors()) + scale_fill_manual(values = GearColors()) +
    scale_shape_manual(values = GearShapes) + scale_linetype_manual(values = GearLines()) +
    guides(color = guide_legend(nrow = 3)) +
    theme(legend.position = "top", legend.title = element_blank(), legend.key.width = rel(1), legend.direction = "vertical")
)


# 1. Helper Functions for Plotting

# --- A. Helper to create a back-transformed effort axis ---
create_effort_axis_scale <- function(gear_name, scale_df, raw_breaks) {
  gear_params <- filter(scale_df, Gear == gear_name)
  mu <- gear_params$Mean
  sigma <- gear_params$SD
  scaled_breaks <- (log(raw_breaks) - mu) / sigma
  formatted_labels <- ifelse(raw_breaks < 5, sprintf("%.1f", raw_breaks), sprintf("%.0f", raw_breaks))
  scale_x_continuous(breaks = scaled_breaks, labels = formatted_labels)
}

# --- B. Main helper function for Richness & Diversity plots ---
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
  plot_data <- data.frame(newdata, Predicted = exp(preds$fit), CI_low = exp(preds$fit - 1.96 * preds$se.fit), CI_high = exp(preds$fit + 1.96 * preds$se.fit))
  
  if (is.null(grouping_var)) {
    plot_data <- plot_data %>% group_by(!!sym(predictor)) %>% summarize(across(c(Predicted, CI_low, CI_high), \(x) mean(x, na.rm = TRUE)))
  }
  return(plot_data)
}


# 2. Abundance Figure

# --- Create Component Plots ---
p_abund_daylight <- estimate_relation(model = AvgMods$AbundInt, by = c("DaylightPercent = [fivenum]", "Gear"), fixed = list(logEffort = 0, Steepness = "Medium"), preserve_range = FALSE) %>%
  ggplot(aes(x = DaylightPercent, y = Predicted, color = Gear, fill = Gear, linetype = Gear)) +
  geom_line(linewidth = 1) + geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.2, linetype = 0) +
  labs(title = "Daylight x Gear", subtitle = "", y = "Abundance", x = "Percent Daylight") +
  scale_x_continuous(breaks = (c(0, 25, 50, 75, 100) - ScaleData$DaylightPercent$Mean) / ScaleData$DaylightPercent$SD, labels = c("0%", "25%", "50%", "75%", "100%")) +
  gear_scales + 
  theme(legend.position.inside = c(0.2,0.8), legend.title = element_blank(), legend.position = "inside", plot.title.position = "plot", plot.title = element_text(vjust = -1))

plot_data_abund_occ <- create_plot_data(AvgMods$AbundInt, AbundIntFit, "Occlusion"); p_abund_occ <- ggplot(plot_data_abund_occ, aes(x = Occlusion, y = Predicted)) + geom_line(linewidth = 1, color = "black") + geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.2, fill = "grey50") + labs(title = "Occlusion", subtitle = "", x = "Percent Occlusion", y = "Abundance") + scale_x_continuous(breaks = (c(17, 22, 27) - ScaleData$Occlusion$Mean) / ScaleData$Occlusion$SD, labels = c("17%", "22%", "27%")) + theme(plot.title.position = "plot", plot.title = element_text(vjust = -1))

plot_data_abund_season <- create_plot_data(AvgMods$AbundInt, AbundIntFit, "Season"); p_abund_season <- ggplot(plot_data_abund_season, aes(x = Season, y = Predicted)) + geom_point(size = 2, color = "black") + geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2, color = "black") + labs(title = "Season", subtitle = "", x = "Season", y = "Abundance") + theme(plot.title.position = "plot", plot.title = element_text(vjust = -1))

newdata_effort_abund <- insight::get_datagrid(AvgMods$AbundInt, by = c("logEffort = [fivenum]", "Gear"), preserve_range = FALSE, Steepness = "Medium")
preds_effort_abund <- predict(AvgMods$AbundInt, newdata = newdata_effort_abund, se.fit = TRUE, type = "link")
plot_data_effort_abund <- data.frame(newdata_effort_abund, Predicted = exp(preds_effort_abund$fit), CI_low = exp(preds_effort_abund$fit - 1.96 * preds_effort_abund$se.fit), CI_high = exp(preds_effort_abund$fit + 1.96 * preds_effort_abund$se.fit))
p_abund_effort_list <- c("Cast Net", "Centipede Net", "Seine") %>%
  map(~ {
    gear_name <- .x; main_title <- if (gear_name == "Cast Net") "Effort x Gear" else ""; x_lab <- switch(gear_name, "Cast Net" = "Throws", "Centipede Net" = "Net-Group-Hours", "Seine" = "Hauls"); raw_breaks <- switch(gear_name, "Cast Net" = exp(seq(log(10), log(20), length.out = 5)), "Centipede Net" = exp(seq(log(8), log(19), length.out = 5)), "Seine" = exp(seq(log(1.0), log(2.0), length.out = 5))); x_scale <- create_effort_axis_scale(gear_name, ScaleDataEffort, raw_breaks)
    
    # Calculate scaled limits for the x-axis
    gear_params <- filter(ScaleDataEffort, Gear == gear_name)
    raw_limits <- range(raw_breaks)
    scaled_limits <- (log(raw_limits) - gear_params$Mean) / gear_params$SD
    
    filter(plot_data_effort_abund, Gear == gear_name) %>% ggplot(aes(x = logEffort, y = Predicted, color = Gear, fill = Gear)) +
      geom_line(linewidth = 1) + geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.2, linetype = 0) +
      labs(title = main_title, subtitle = gear_name, x = x_lab, y = "Abundance\n(log scale)") + x_scale +
      scale_y_log10(breaks = c(1, 5, 25, 125, 625), labels = scales::label_number(accuracy = 1)) +
      coord_cartesian(xlim = scaled_limits, ylim = c(1, 625), expand = expansion(mult = 0)) + # Apply corrected x-axis limits
      scale_color_manual(values = GearColors()[gear_name]) + scale_fill_manual(values = GearColors()[gear_name]) +
      theme(plot.subtitle = element_text(hjust = 1), plot.title = element_text(vjust = -1), legend.position = "none", plot.title.position = "plot")
  })

p_abund_effort <- plot_grid(plotlist = p_abund_effort_list, ncol = 1, align = "v")

# --- Assemble Abundance Grid ---
left_col_abund <- plot_grid(p_abund_daylight, p_abund_occ, p_abund_season, ncol = 1, align = "v")
final_plot_abund <- plot_grid(left_col_abund, p_abund_effort, ncol = 2) + plot_theme


# 3. Richness Figure

# --- Interaction Plots ---
newdata_effort_rich <- insight::get_datagrid(AvgMods$RichInt, by = c("logEffort = [fivenum]", "Gear"), preserve_range = FALSE, Steepness = "Medium")
newdata_effort_rich$Site <- NA
preds_effort_rich <- predict(AvgMods$RichInt, newdata = newdata_effort_rich, se.fit = TRUE, type = "link")
plot_data_effort_rich <- data.frame(newdata_effort_rich, Predicted = exp(preds_effort_rich$fit), CI_low = exp(preds_effort_rich$fit - 1.96 * preds_effort_rich$se.fit), CI_high = exp(preds_effort_rich$fit + 1.96 * preds_effort_rich$se.fit))
p_rich_effort_list <- c("Cast Net", "Centipede Net", "Seine") %>%
  map(~ {
    gear_name <- .x; main_title <- if (gear_name == "Cast Net") "Effort x Gear" else ""; x_lab <- switch(gear_name, "Cast Net" = "Throws", "Centipede Net" = "Net-Group-Hours", "Seine" = "Hauls"); raw_breaks <- switch(gear_name, "Cast Net" = exp(seq(log(10), log(20), length.out = 5)), "Centipede Net" = exp(seq(log(8), log(19), length.out = 5)), "Seine" = exp(seq(log(1.0), log(2.0), length.out = 5))); x_scale <- create_effort_axis_scale(gear_name, ScaleDataEffort, raw_breaks)
    
    gear_params <- filter(ScaleDataEffort, Gear == gear_name)
    raw_limits <- range(raw_breaks)
    scaled_limits <- (log(raw_limits) - gear_params$Mean) / gear_params$SD
    
    ggplot(filter(plot_data_effort_rich, Gear == gear_name), aes(x = logEffort, y = Predicted, color = Gear, fill = Gear)) +
      geom_line(linewidth = 1) + geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.2, linetype = 0) +
      labs(title = main_title, subtitle = gear_name, x = x_lab, y = "Richness") + x_scale +
      coord_cartesian(xlim = scaled_limits, ylim = c(0, 30), expand = expansion(mult = 0)) +
      scale_color_manual(values = GearColors()[gear_name]) + scale_fill_manual(values = GearColors()[gear_name]) +
      theme(plot.subtitle = element_text(hjust = 1), legend.position = "none", plot.title.position = "plot", plot.title = element_text(vjust = -1))
  })

# --- Main Effect Plots ---
plot_data_rich_mud <- create_plot_data(AvgMods$RichInt, RichIntFit, "MudDominant"); p_rich_mud <- ggplot(plot_data_rich_mud, aes(x = MudDominant, y = Predicted)) + geom_point(size = 2, color = "black") + geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2, color = "black") + labs(title = "Substrate", subtitle = "", x = "Substrate", y = "Richness") + scale_x_discrete(labels = c("False" = "Other", "True" = "Mud Dominant"), limits = rev)+ 
  scale_y_continuous(breaks = c(3,5,7,9), minor_breaks = NULL) + coord_cartesian(ylim = c(2,10)) +
  theme(plot.title.position = "plot", plot.title = element_text(vjust = -1))

plot_data_season <- create_plot_data(AvgMods$RichInt, RichIntFit, "Season"); p_rich_season <- ggplot(plot_data_season, aes(x = Season, y = Predicted)) + geom_point(size = 2, color = "black") + geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2, color = "black") + labs(title = "Season", subtitle = "", x = "Season", y = "Richness") + 
  scale_y_continuous(breaks = c(3,5,7,9), minor_breaks = NULL) + coord_cartesian(ylim = c(2,10)) +
  theme(plot.title.position = "plot", plot.title = element_text(vjust = -1))

plot_data_daylight <- create_plot_data(AvgMods$RichInt, RichIntFit, "DaylightPercent"); p_rich_daylight <- ggplot(plot_data_daylight, aes(x = DaylightPercent, y = Predicted)) + geom_line(linewidth = 1, color = "black") + geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.2, fill = "grey50") + labs(title = "Daylight", subtitle = "", x = "Percent Daylight", y = "Richness") + scale_x_continuous(breaks = (c(0, 25, 50, 75, 100) - ScaleData$DaylightPercent$Mean) / ScaleData$DaylightPercent$SD, labels = c("0%", "25%", "50%", "75%", "100%")) + 
  scale_y_continuous(breaks = c(3,5,7,9), minor_breaks = NULL) + coord_cartesian(ylim = c(2,10)) +
  theme(plot.title.position = "plot", plot.title = element_text(vjust = -1))

plot_data_secchi <- create_plot_data(AvgMods$RichInt, RichIntFit, "SecchiDepth"); p_rich_secchi <- ggplot(plot_data_secchi, aes(x = SecchiDepth, y = Predicted)) + geom_line(linewidth = 1, color = "black") + geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.2, fill = "grey50") + labs(title = "Secchi Depth", subtitle = "", x = "Secchi Depth (cm)", y = "Richness") + scale_x_continuous(breaks = (c(20, 75, 125, 175) - ScaleData$SecchiDepth$Mean) / ScaleData$SecchiDepth$SD, labels = c("20", "75", "125", "175")) + 
  scale_y_continuous(breaks = c(3,5,7,9), minor_breaks = NULL) + coord_cartesian(ylim = c(2,10)) +
  theme(plot.title.position = "plot", plot.title = element_text(vjust = -1))

plot_data_occ <- create_plot_data(AvgMods$RichInt, RichIntFit, "Occlusion"); p_rich_occ <- ggplot(plot_data_occ, aes(x = Occlusion, y = Predicted)) + geom_line(linewidth = 1, color = "black") + geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.2, fill = "grey50") + labs(title = "Occlusion", subtitle = "", x = "Percent Occlusion", y = "Richness") + scale_x_continuous(breaks = (c(17, 22, 27) - ScaleData$Occlusion$Mean) / ScaleData$Occlusion$SD, labels = c("17%", "22%", "27%")) + 
  scale_y_continuous(breaks = c(3,5,7,9), minor_breaks = NULL) + coord_cartesian(ylim = c(2,10)) +
  theme(plot.title.position = "plot", plot.title = element_text(vjust = -1))

plot_data_steepness <- create_plot_data(AvgMods$RichInt, RichIntFit, "Steepness"); p_rich_steepness <- ggplot(plot_data_steepness, aes(x = Steepness, y = Predicted)) + geom_point(size = 2, color = "black") + geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2, color = "black") + labs(title = "Steepness", subtitle = "", x = "Steepness", y = "Richness") + scale_x_discrete(limits = c("Low", "Medium", "High")) + 
  scale_y_continuous(breaks = c(3,5,7,9), minor_breaks = NULL) + coord_cartesian(ylim = c(2,10)) +
  theme(plot.title.position = "plot", plot.title = element_text(vjust = -1))

# --- Assemble Richness Grid ---
final_plot_rich <- plot_grid(p_rich_effort_list[[1]], p_rich_season, p_rich_daylight, p_rich_effort_list[[2]], p_rich_steepness, p_rich_occ, p_rich_effort_list[[3]], p_rich_mud, p_rich_secchi, ncol = 3, align = "hv", axis = "tblr", rel_widths = c(1,0.85,1)) + plot_theme


# 4. Diversity Figures

# --- Shannon Diversity ---
plot_data_div_gear <- create_plot_data(AvgMods$DivInt, DivIntFit, "Gear")
p_div_gear <- ggplot(plot_data_div_gear, aes(x = Gear, y = Predicted, color = Gear)) + geom_point(size = 3) + geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2, linewidth = 1) + labs(title = "Gear", subtitle = "", x = "Gear", y = "Shannon Diversity") +  gear_scales + theme(legend.position = "none", plot.title.position = "plot", plot.title = element_text(vjust = -1))

plot_data_div_effort <- create_plot_data(AvgMods$DivInt, DivIntFit, "logEffort")

p_div_effort <- ggplot(plot_data_div_effort, aes(x = logEffort, y = Predicted)) + 
  geom_line(linewidth = 1, color = "black") + 
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.2, fill = "grey50") + 
  labs(
    title = "Effort", 
    subtitle = "", 
    x = "Relative Effort (SD)",
    y = "Shannon Diversity"
  ) + 
  theme(plot.title.position = "plot", plot.title = element_text(vjust = -1))

plot_data_div_mud <- create_plot_data(AvgMods$DivInt, DivIntFit, "MudDominant")
p_div_mud <- ggplot(plot_data_div_mud, aes(x = MudDominant, y = Predicted)) + geom_point(size = 2, color = "black") + geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2, color = "black") + labs(title = "Substrate", subtitle = "", x = "Substrate", y = "Shannon Diversity") + scale_x_discrete(labels = c("False" = "Other", "True" = "Mud Dominant"), limits = rev) + theme(plot.title.position = "plot", plot.title = element_text(vjust = -1))

plot_data_div_daylight <- create_plot_data(AvgMods$DivInt, DivIntFit, "DaylightPercent")
p_div_daylight <- ggplot(plot_data_div_daylight, aes(x = DaylightPercent, y = Predicted)) + geom_line(linewidth = 1, color = "black") + geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.2, fill = "grey50") + labs(title = "Daylight", subtitle = "", x = "Percent Daylight", y = "Shannon Diversity") + scale_x_continuous(breaks = (c(0, 25, 50, 75, 100) - ScaleData$DaylightPercent$Mean) / ScaleData$DaylightPercent$SD, labels = c("0%", "25%", "50%", "75%", "100%"))+theme(plot.title.position = "plot", plot.title = element_text(vjust = -1))

plot_data_div_occ <- create_plot_data(AvgMods$DivInt, DivIntFit, "Occlusion")
p_div_occ <- ggplot(plot_data_div_occ, aes(x = Occlusion, y = Predicted)) + geom_line(linewidth = 1, color = "black") + geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.2, fill = "grey50") + labs(title = "Occlusion", subtitle = "", x = "Percent Occlusion", y = "Shannon Diversity") + scale_x_continuous(breaks = (c(17, 22, 27) - ScaleData$Occlusion$Mean) / ScaleData$Occlusion$SD, labels = c("17%", "22%", "27%")) + theme(plot.title.position = "plot", plot.title = element_text(vjust = -1))

plot_data_div_steepness <- create_plot_data(AvgMods$DivInt, DivIntFit, "Steepness")
p_div_steepness <- ggplot(plot_data_div_steepness, aes(x = Steepness, y = Predicted)) + geom_point(size = 2, color = "black") + geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2, color = "black") + labs(title = "Steepness", subtitle = "", x = "Steepness", y = "Shannon Diversity") + scale_x_discrete(limits = c("Low", "Medium", "High")) + theme(plot.title.position = "plot", plot.title = element_text(vjust = -1))


final_plot_div <- plot_grid(p_div_daylight, p_div_effort, p_div_gear, p_div_occ, p_div_steepness, p_div_mud,ncol = 3, align = "hv") + plot_theme

# --- Simpson Diversity ---
plot_data_sim_gear <- create_plot_data(AvgMods$SimInt, SimIntFit, "Gear"); p_sim_gear <- ggplot(plot_data_sim_gear, aes(x = Gear, y = Predicted, color = Gear)) + geom_point(size = 3) + geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2, linewidth = 1) + labs(title = "Gear", subtitle = "", x = "Gear", y = "Simpson Diversity") + gear_scales + theme(legend.position = "none", plot.title.position = "plot", plot.title = element_text(vjust = -1))

plot_data_sim_daylight <- create_plot_data(AvgMods$SimInt, SimIntFit, "DaylightPercent"); p_sim_daylight <- ggplot(plot_data_sim_daylight, aes(x = DaylightPercent, y = Predicted)) + geom_line(linewidth = 1, color = "black") + geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.2, fill = "grey50") + labs(title = "Daylight", subtitle = "", x = "Percent Daylight", y = "Simpson Diversity") + scale_x_continuous(breaks = (c(0, 25, 50, 75, 100) - ScaleData$DaylightPercent$Mean) / ScaleData$DaylightPercent$SD, labels = c("0%", "25%", "50%", "75%", "100%")) + theme(plot.title.position = "plot", plot.title = element_text(vjust = -1))

final_plot_sim <- plot_grid(p_sim_daylight, p_sim_gear, ncol = 2, align = "hv", axis = "tblr", rel_widths = c(1.25, 1)) + plot_theme


# 5. Display Final Plots
final_plot_abund
final_plot_rich
final_plot_div
final_plot_sim


# ggsave(plot = final_plot_abund, "output/plots/MarginalEffects_Abundance.eps", device = cairo_ps,  width = 5.62, height = 7.25, units = "in")
# ggsave(plot = final_plot_abund, "output/plots/MarginalEffects_Abundance.png",  width = 5.62, height = 7.25, units = "in", dpi = 600)
# 
# ggsave(plot = final_plot_rich, "output/plots/MarginalEffects_Richness.eps", device = cairo_ps,  width = 5.62, height = 7.25, units = "in")
# ggsave(plot = final_plot_rich, "output/plots/MarginalEffects_Richness.png",  width = 5.62, height = 7.25, units = "in", dpi = 600)
# 
# ggsave(plot = final_plot_div, "output/plots/MarginalEffects_Shannon.eps", device = cairo_ps,  width = 5.62, height = 4.21, units = "in")
# ggsave(plot = final_plot_div, "output/plots/MarginalEffects_Shannon.png",  width = 5.62, height = 4.21, units = "in", dpi = 600)
# 
# ggsave(plot = final_plot_sim, "output/plots/MarginalEffects_Simpson.eps", device = cairo_ps,  width = 5.62, height = 3.75, units = "in")
# ggsave(plot = final_plot_sim, "output/plots/MarginalEffects_Simpson.png",  width = 5.62, height = 3.75, units = "in", dpi = 600)


# --- 6. 2x2 Metric ~ Gear Summary ---

p_abund_effort_combined <- ggplot(plot_data_effort_abund, aes(x = logEffort, y = Predicted, color = Gear, fill = Gear, linetype = Gear)) +
  geom_line(linewidth = 1) + 
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.1, linetype = 0) +
  labs(title = "Abundance", x = "Relative Effort (SD)", y = "Abundance\n(log scale)") + 
  scale_y_log10(labels = scales::label_number(accuracy = 1)) +
  gear_scales + 
  theme(legend.position = "none", plot.title.position = "plot", plot.title = element_text(vjust = -1))

p_rich_effort_combined <- ggplot(plot_data_effort_rich, aes(x = logEffort, y = Predicted, color = Gear, fill = Gear, linetype = Gear)) +
  geom_line(linewidth = 1) + 
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.1, linetype = 0) +
  labs(title = "Richness", x = "Relative Effort (SD)", y = "Richness\n(log scale)") + 
  scale_y_log10(labels = scales::label_number(accuracy = 1)) +
  gear_scales +
  theme(legend.position = "none", plot.title.position = "plot", plot.title = element_text(vjust = -1))

p_div_gear_clean <- p_div_gear + labs(title = "Shannon Diversity", x = "Gear", y = "Shannon Diversity") + theme(legend.position = "none")
p_sim_gear_clean <- p_sim_gear + labs(title = "Simpson Diversity", x = "Gear", y = "Simpson Diversity") + theme(legend.position = "none")

# Layout components
row1 <- plot_grid(p_abund_effort_combined, p_rich_effort_combined, ncol = 2, align = "hv", axis = "tblr")
row2 <- plot_grid(p_div_gear_clean, p_sim_gear_clean, ncol = 2, align = "hv", axis = "tblr")

shared_legend_horiz <- get_legend(
  p_abund_effort_combined + 
    theme(legend.position = "bottom", legend.title = element_blank(), legend.key.width = rel(2))
)

# Assemble with legend in the middle
final_metric_gear_2x2_legend <- plot_grid(
  row1, shared_legend_horiz, row2,
  ncol = 1, rel_heights = c(1, 0.15, 1)
)

final_metric_gear_2x2_legend
# ggsave(plot = final_metric_gear_2x2_legend, "output/plots/Metric_Gear_2x2.png", width = 5.62, height = 6.5, units = "in", dpi = 300)
