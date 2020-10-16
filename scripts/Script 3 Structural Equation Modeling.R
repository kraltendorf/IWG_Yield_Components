# Project: IWG_Yield_Components
# Script 3 - Structural Equation Modeling
# Author: Kayla R. Altendorf
# Date: 05/25/2020

# Required Packages::
library("lavaan")
library("lavaanPlot")

# Using the same vectors and data as used in previous scripts (e.g. envs, path, emmeans_all, etc.)

#### Step 1: Scale Data ####
# scale data to have a mean of 0 and sd 1, and make as data frame
dat_stp17 <- as.data.frame(scale(emmeans_all[[1]][,-14:-15])) # removing loc and year identifiers
dat_stp18 <- as.data.frame(scale(emmeans_all[[2]][,-14:-15]))
dat_tli17 <- as.data.frame(scale(emmeans_all[[3]][,-14:-15]))
dat_tli18 <- as.data.frame(scale(emmeans_all[[4]][,-14:-15])) 



#### Step 2: Set Up the Initial Model  ####
# the term to the left of the "~" is the trait explained by the traits on the right
# each term is fit with an acronym that defines the estimate for that path
# naming these makes calculating the effects easier 

# KEY:
#1. STD = stem_diameter
#2. TGW = thousand_grain_weight
#3. FLA = flag_leaf_area
#4. ANT = anthesis_score
#5. SPS = spikelets_per_spike
#6. FPS = florets_per_spikelet
#7. RNC = reproductive_tiller_ct
#8. YPS = yield_per_spike
#9. YPP = yield_per_plant
#10. FSU = floret_site_utilization
#11. SPD = spikelet_density
#12. HIT = height
#13. EMP = emergence_percent

# example: STD_TGW is term defining the effect of stem_diameter on thousand_grain_weight
# example: STD_YPS_Indirect := (STD_FSU * FSU_YPS) + (STD_TGW * TGW_YPS) is the total indirect effect of 
# stem diameter on yield per spike via its effects on floret site utilization, and yield per spike

initial_model <- ' 
thousand_grain_weight ~ STD_TGW*stem_diameter + FLA_TGW*flag_leaf_area + ANT_TGW*anthesis_score + SPD_TGW*spikelet_density + HIT_TGW*height + EMP_TGW*emergence_percent + RTC_TGW*reproductive_tiller_ct
floret_site_utilization ~ STD_FSU*stem_diameter + FLA_FSU*flag_leaf_area + ANT_FSU*anthesis_score + SPD_FSU*spikelet_density + HIT_FSU*height + EMP_FSU*emergence_percent + RTC_FSU*reproductive_tiller_ct
yield_per_spike ~ FSU_YPS*floret_site_utilization + TGW_YPS*thousand_grain_weight + SPS_YPS*spikelets_per_spike + FPS_YPS*florets_per_spikelet 
yield_per_plant ~ YPS_YPP*yield_per_spike + RTC_YPP*reproductive_tiller_ct

# code for calculating direct and indirect effects
# the variables are named by the symbols ":="
# first for yield per spike (YPS)
# FPS 
FPS_YPS_Direct := FPS_YPS

# ANT
ANT_YPS_Indirect := (ANT_FSU * FSU_YPS) + (ANT_TGW * TGW_YPS)

# STD
STD_YPS_Indirect := (STD_FSU * FSU_YPS) + (STD_TGW * TGW_YPS)

# SPD 
SPD_YPS_Indirect := (SPD_FSU * FSU_YPS) + (SPD_TGW * TGW_YPS)

# HIT 
HIT_YPS_Indirect := (HIT_FSU * FSU_YPS) + (HIT_TGW * TGW_YPS)

# FLA 
FLA_YPS_Indirect := (FLA_FSU * FSU_YPS) + (FLA_TGW * TGW_YPS)

# RTC
RTC_YPS_Indirect := (RTC_FSU * FSU_YPS) + (RTC_TGW * TGW_YPS)

# FSU
FSU_YPS_Direct := FSU_YPS

# TGW 
TGW_YPS_Direct := TGW_YPS

# SPS
SPS_YPS_Direct := SPS_YPS

# EMP
EMP_YPS_Indirect := (EMP_TGW * TGW_YPS) + (EMP_FSU * FSU_YPS)

## second for yield per plant
# FPS
FPS_YPP_Indirect := (FPS_YPS * YPS_YPP) 

# HIT
HIT_YPP_Indirect := (HIT_FSU * FSU_YPS * YPS_YPP) + (HIT_TGW * TGW_YPS * YPS_YPP)

# STD
STD_YPP_Indirect := (STD_FSU * FSU_YPS * YPS_YPP) + (HIT_TGW * TGW_YPS * YPS_YPP)

# ANT
ANT_YPP_Indirect := (ANT_FSU * FSU_YPS * YPS_YPP) + (ANT_TGW * TGW_YPS * YPS_YPP)

# EMP
EMP_YPP_Indirect := (EMP_FSU * FSU_YPS * YPS_YPP) + (EMP_TGW * TGW_YPS * YPS_YPP)

# SPD
SPD_YPP_Indirect := (SPD_FSU * FSU_YPS * YPS_YPP) + (SPD_TGW * TGW_YPS * YPS_YPP)

# FLA
FLA_YPP_Indirect := (FLA_FSU * FSU_YPS * YPS_YPP) + (FLA_TGW * TGW_YPS * YPS_YPP)

# FSU
FSU_YPP_Indirect := (FSU_YPS * YPS_YPP)

# TGW
TGW_YPP_Indirect := (TGW_YPS * YPS_YPP)

# SPS
SPS_YPP_Indirect := (SPS_YPS * YPS_YPP)

# RTC
RTC_YPP_Direct := (RTC_YPP)
RTC_YPP_Indirect := (RTC_FSU * FSU_YPS * YPS_YPP) + (RTC_TGW * TGW_YPS * YPS_YPP)

# YPS_YPP
YPS_YPP_Direct := (YPS_YPP)' 



#### Step 3: Run the Model for STP 2017 ####

# run model in lavaan using the sem function
mod <- sem(initial_model, data=dat_stp17, fixed.x = FALSE)  

# visualize the model
lavaanPlot(model = mod, node_options = list(shape = "box"), stand = T, edge_options = list(color = "gray"), coef = T, covs = T, stars = T,  sig=0.05)

# look at model summary to assess fit
summary(mod, fit.measures = TRUE, standardized = T)
# things to look at: RMSEA, CFI, chisq test results

# results from attempt 1
# RMSEA 0.086
# CFI 0.963

# check out the modification indices
modificationIndices(mod, sort. = T)

# according to the modification indices, it is suggested to add floret_site_utilization ~ florets_per_spikelet

#### Step 4: Modify the Model ####
modified_model <- ' 
thousand_grain_weight ~ STD_TGW*stem_diameter + FLA_TGW*flag_leaf_area + ANT_TGW*anthesis_score + SPD_TGW*spikelet_density + HIT_TGW*height + EMP_TGW*emergence_percent + RTC_TGW*reproductive_tiller_ct
floret_site_utilization ~ STD_FSU*stem_diameter + FLA_FSU*flag_leaf_area + ANT_FSU*anthesis_score + SPD_FSU*spikelet_density + HIT_FSU*height + EMP_FSU*emergence_percent + RTC_FSU*reproductive_tiller_ct + FPS_FSU*florets_per_spikelet ### changed here
yield_per_spike ~ FSU_YPS*floret_site_utilization + TGW_YPS*thousand_grain_weight + SPS_YPS*spikelets_per_spike + FPS_YPS*florets_per_spikelet 
yield_per_plant ~ YPS_YPP*yield_per_spike + RTC_YPP*reproductive_tiller_ct

# first for yield per spike (YPS)
# FPS 
FPS_YPS_Direct := FPS_YPS
FPS_YPS_Indirect := FPS_FSU * FSU_YPS ### changed here

# ANT
ANT_YPS_Indirect := (ANT_FSU * FSU_YPS) + (ANT_TGW * TGW_YPS)

# STD
STD_YPS_Indirect := (STD_FSU * FSU_YPS) + (STD_TGW * TGW_YPS)

# SPD 
SPD_YPS_Indirect := (SPD_FSU * FSU_YPS) + (SPD_TGW * TGW_YPS)

# HIT 
HIT_YPS_Indirect := (HIT_FSU * FSU_YPS) + (HIT_TGW * TGW_YPS)

# FLA 
FLA_YPS_Indirect := (FLA_FSU * FSU_YPS) + (FLA_TGW * TGW_YPS)

# RTC
RTC_YPS_Indirect := (RTC_FSU * FSU_YPS) + (RTC_TGW * TGW_YPS)

# FSU
FSU_YPS_Direct := FSU_YPS

# TGW 
TGW_YPS_Direct := TGW_YPS

# SPS
SPS_YPS_Direct := SPS_YPS

# EMP
EMP_YPS_Indirect := (EMP_TGW * TGW_YPS) + (EMP_FSU * FSU_YPS)

## second for yield per plant
# FPS
FPS_YPP_Indirect := (FPS_FSU * FSU_YPS * YPS_YPP) + (FPS_YPS * YPS_YPP) ### changed here

# HIT
HIT_YPP_Indirect := (HIT_FSU * FSU_YPS * YPS_YPP) + (HIT_TGW * TGW_YPS * YPS_YPP)

# STD
STD_YPP_Indirect := (STD_FSU * FSU_YPS * YPS_YPP) + (HIT_TGW * TGW_YPS * YPS_YPP)

# ANT
ANT_YPP_Indirect := (ANT_FSU * FSU_YPS * YPS_YPP) + (ANT_TGW * TGW_YPS * YPS_YPP)

# EMP
EMP_YPP_Indirect := (EMP_FSU * FSU_YPS * YPS_YPP) + (EMP_TGW * TGW_YPS * YPS_YPP)

# SPD
SPD_YPP_Indirect := (SPD_FSU * FSU_YPS * YPS_YPP) + (SPD_TGW * TGW_YPS * YPS_YPP)

# FLA
FLA_YPP_Indirect := (FLA_FSU * FSU_YPS * YPS_YPP) + (FLA_TGW * TGW_YPS * YPS_YPP)

# FSU
FSU_YPP_Indirect := (FSU_YPS * YPS_YPP)

# TGW
TGW_YPP_Indirect := (TGW_YPS * YPS_YPP)

# SPS
SPS_YPP_Indirect := (SPS_YPS * YPS_YPP)

# RTC
RTC_YPP_Direct := (RTC_YPP)
RTC_YPP_Indirect := (RTC_FSU * FSU_YPS * YPS_YPP) + (RTC_TGW * TGW_YPS * YPS_YPP)

# YPS_YPP
YPS_YPP_Direct := (YPS_YPP)'

# run model in lavaan using the sem function
mod <- sem(modified_model, data=dat_stp17, fixed.x = FALSE)  

# visualize the model
lavaanPlot(model = mod, node_options = list(shape = "box"), stand = T, edge_options = list(color = "gray"), coef = T, covs = T, stars = T,  sig=0.05)

# look at model summary to assess fit
summary(mod, fit.measures = TRUE, standardized = T)
# things to look at: RMSEA, CFI, chisq test results

# results from attempt 2
# RMSEA 0.72
# CFI 0.975 
# this is an acceptable fit

# look at r2 value for variance explained by exogenous variables
inspect(mod, 'r2')

# make a dataframe of path coefficients
estimates <- as.data.frame(parameterEstimates(mod))
estimates_output <- estimates %>% 
  filter(op %in% c(":=", "~")) # filter out ~ paths and defined paths := 
# write output
write.csv(estimates_output, paste(path, "/estimates_stp17.csv", sep = ""))

# repare the dataframe to be used in ggplot
estimates1 <- estimates %>% 
  filter(op == ":=") %>%   # keep only the assigned variables
  mutate(est = ifelse(test = (pvalue > 0.05), yes = 0, # and if it's not significant change the estimate to zero 
                      no = est)) %>%  # that way it's not in the figure
  mutate(ci.upper = ifelse(test = (pvalue > 0.05), yes = NA,
                           no = ci.upper)) %>%
  mutate(ci.lower = ifelse(test = (pvalue > 0.05), yes = NA,
                           no = ci.lower)) %>%
  separate(col = label, into = c("trait", "on", "type"), sep = "_") 

# rename traits with their pub names 
# make a key to left merge the names together
trait <- c("ANT", "EMP", "FLA", "FSU", "FPS", "HIT", "RTC", "SPD", "SPS", "STD", "TGW", "YPP", "YPS")
# remove units from pub names for this figure since the paths are essentially unitless
pub_names_no_units <- c("Anthesis Score", "Spike Emergence", "Flag Leaf Area", "Floret Site Utilization", "Florets per Spikelet", 
               "Height", "Reproductive Tillers", "Spikelet Density", "Spikelets per Spike",
               "Stem Diameter", "Thousand Grain Weight", "Yield per Plant", "Yield per Spike")

trait_key <- data.frame(pub_names_no_units, trait)
estimates2 <- left_join(estimates1, trait_key, by = "trait")

# rename dataframe for this env
estimates2$loc <- rep("STP", length(estimates2$trait))
estimates2$year <- rep("2017", length(estimates2$trait))

# write it out for the year
estimates_stp_17 <- estimates2 # will call upon this later



#### Step 3: Run the Model for STP 2018  ####
# using the same initial model from STP 2017, and editing as needed

# run model in lavaan using the sem function
# since the same model fit both stp 17 and stp 18, use the same code and iterate through
mod <- sem(model, data=dat_stp18, fixed.x = FALSE)  

# visualize the model
lavaanPlot(model = mod, node_options = list(shape = "box"), stand = T, edge_options = list(color = "gray"), coef = T, covs = T, stars = T,  sig=0.05)

# look at model summary to assess fit
summary(mod, fit.measures = TRUE, standardized = T)

# check out the modification indices
modificationIndices(mod, sort. = T)

# results from attempt 1
# RMSEA 0.115
# CFI 0.924
# modification indices again suggest adding: floret_site_utilization ~ floret_ct_avg
mod <- sem(modified_model, data=dat_stp18, fixed.x = FALSE)  # now using the same modified model

# visualize the model
lavaanPlot(model = mod, node_options = list(shape = "box"), stand = T, edge_options = list(color = "gray"), coef = T, covs = T, stars = T,  sig=0.05)

# look at model summary to assess fit
summary(mod, fit.measures = TRUE, standardized = T)

# results from attempt 2
# RMSEA 0.059
# CFI 0.98 # this yielded a major improvement in model fit -- meaning it's the same model as STP 2017

# look at r2 value for variance explained by exogenous variables
inspect(mod, 'r2')

# make a dataframe of path coefficients
estimates <- as.data.frame(parameterEstimates(mod))

# export full path results for supplementary tables
estimates_output <- estimates %>% 
  filter(op %in% c(":=", "~")) # filter out ~ paths and defined paths := 

# check year at the end of file path 
write.csv(estimates_output, paste(path, "/estimates_stp18.csv", sep = ""))

# prepare estimates for ggplot
estimates1 <- estimates %>% 
  filter(op == ":=") %>% 
  mutate(est = ifelse(test = (pvalue > 0.05), yes = 0,
                      no = est)) %>% 
  mutate(ci.upper = ifelse(test = (pvalue > 0.05), yes = NA,
                           no = ci.upper)) %>%
  mutate(ci.lower = ifelse(test = (pvalue > 0.05), yes = NA,
                           no = ci.lower)) %>%
  separate(col = label, into = c("trait", "on", "type"), sep = "_") 

# using pre-defined trait_key
estimates2 <- left_join(estimates1, trait_key, by = "trait")

# rename dataframe for this env - edit year
estimates2$loc <- rep("STP", length(estimates2$trait))
estimates2$year <- rep("2018", length(estimates2$trait))

# write it out for the year
estimates_stp_18 <- estimates2



#### Step 4: Run the Model for TLI 2017 ####         
# run using the initial_model
mod <- sem(initial_model, data=dat_tli17, fixed.x = FALSE)  # run model
lavaanPlot(model = mod, node_options = list(shape = "box"), stand = T, edge_options = list(color = "gray"), coef = T, covs = T, stars = T,  sig=0.05)
summary(mod, fit.measures = TRUE, standardized = T)
modificationIndices(mod, sort. = T)

# attempt 1
# RMSEA 0.097
# CFI 0.946 # pretty good model fit, but modification indices suggest add FPS to explain FSU
# while there are other paths higher up on the list, this is the first path that can 
# be logically added (assuming components effect eachother in a sequential order)

# run it again using the modified model
mod <- sem(modified_model, data=dat_tli17, fixed.x = FALSE)  # run model
lavaanPlot(model = mod, node_options = list(shape = "box"), stand = T, edge_options = list(color = "gray"), coef = T, covs = T, stars = T,  sig=0.05)
summary(mod, fit.measures = TRUE, standardized = T)
modificationIndices(mod, sort. = T)

# attempt 2
# RMSEA 0.086
# CFI 0.959 # add height to explain yield per spike

#### Step 5: Modify the Model ####

modified_model2 <- ' 
thousand_grain_weight ~ STD_TGW*stem_diameter + FLA_TGW*flag_leaf_area  + ANT_TGW*anthesis_score + SPD_TGW*spikelet_density + HIT_TGW*height + EMP_TGW*emergence_percent + FSU_TGW*floret_site_utilization + RTC_TGW*reproductive_tiller_ct
floret_site_utilization ~ STD_FSU*stem_diameter + FLA_FSU*flag_leaf_area  + ANT_FSU*anthesis_score + SPD_FSU*spikelet_density + HIT_FSU*height + EMP_FSU*emergence_percent + RTC_FSU*reproductive_tiller_ct + FPS_FSU*florets_per_spikelet
yield_per_spike ~ FSU_YPS*floret_site_utilization + TGW_YPS*thousand_grain_weight + SPS_YPS*spikelets_per_spike + FPS_YPS*florets_per_spikelet + HIT_YPS*height
yield_per_plant ~ YPS_YPP*yield_per_spike + RTC_YPP*reproductive_tiller_ct

## for YPS
# FPS
FPS_YPS_Direct := FPS_YPS
FPS_YPS_Indirect := FPS_FSU * FSU_YPS

# ANT
ANT_YPS_Indirect := (ANT_FSU * FSU_YPS) + (ANT_TGW * TGW_YPS)

# STD
STD_YPS_Indirect := (STD_FSU * FSU_YPS) + (STD_TGW * TGW_YPS)

# SPD 
SPD_YPS_Indirect := (SPD_FSU * FSU_YPS) + (SPD_TGW * TGW_YPS)

# HIT 
HIT_YPS_Indirect := (HIT_FSU * FSU_YPS) + (HIT_TGW * TGW_YPS)
HIT_YPS_Direct := (HIT_YPS)

# FLA 
FLA_YPS_Indirect := (FLA_FSU * FSU_YPS) + (FLA_TGW * TGW_YPS)

# FSU
FSU_YPS_Direct := FSU_YPS
FSU_YPS_Indirect := FSU_TGW * TGW_YPS

# TGW 
TGW_YPS_Direct := TGW_YPS

# SPS
SPS_YPS_Direct := SPS_YPS

# EMP
EMP_YPS_Indirect := (EMP_TGW * TGW_YPS) + (EMP_FSU * FSU_YPS)

# RTC
RTC_YPS_Indirect := (RTC_TGW * TGW_YPS) + (RTC_FSU * FSU_YPS)

## for YPP
# FPS
FPS_YPP_Indirect := (FSU_YPS * YPS_YPP) + (FPS_YPS * YPS_YPP) # edit here FPS_FSU

# HIT
HIT_YPP_Indirect := (HIT_FSU * FSU_YPS * YPS_YPP) + (HIT_TGW * TGW_YPS * YPS_YPP) + (HIT_YPS * YPS_YPP)

# STD
STD_YPP_Indirect := (STD_FSU * FSU_YPS * YPS_YPP) + (HIT_TGW * TGW_YPS * YPS_YPP)

# ANT
ANT_YPP_Indirect := (ANT_FSU * FSU_YPS * YPS_YPP) + (ANT_TGW * TGW_YPS * YPS_YPP)

# EMP
EMP_YPP_Indirect := (EMP_FSU * FSU_YPS * YPS_YPP) + (EMP_TGW * TGW_YPS * YPS_YPP)

# SPD
SPD_YPP_Indirect := (SPD_FSU * FSU_YPS * YPS_YPP) + (SPD_TGW * TGW_YPS * YPS_YPP)

# FLA
FLA_YPP_Indirect := (FLA_FSU * FSU_YPS * YPS_YPP) + (FLA_TGW * TGW_YPS * YPS_YPP)

# FSU
FSU_YPP_Indirect := (FSU_YPS * YPS_YPP) + (FSU_TGW * TGW_YPS * YPS_YPP)

# TGW
TGW_YPP_Indirect := (TGW_YPS * YPS_YPP)

# SPS
SPS_YPP_Indirect := (SPS_YPS * YPS_YPP)

# RTC
RTC_YPP_Direct := (RTC_YPP)
RTC_YPP_Indirect := (RTC_FSU * FSU_YPS * YPS_YPP) + (RTC_TGW * TGW_YPS * YPS_YPP)

# YPS_FPY
YPS_YPP_Direct := (YPS_YPP)'

#### Step 5: Test the Model  ####

mod <- sem(modified_model2, data=dat_tli17, fixed.x = FALSE)  # run model
lavaanPlot(model = mod, node_options = list(shape = "box"), stand = T, edge_options = list(color = "gray"), coef = T, covs = T, stars = T,  sig=0.05)
summary(mod, fit.measures = TRUE, standardized = T)
modificationIndices(mod, sort. = T)

# attempt 3
# RMSEA 0.069
# CFI 0.976
# looks good, proceed

inspect(mod, 'r2')
estimates <- as.data.frame(parameterEstimates(mod))
estimates_output <- estimates %>% 
  filter(op %in% c(":=", "~"))

write.csv(estimates_output, paste(path, "/estimates_tli17.csv", sep = ""))

estimates1 <- estimates %>% 
  filter(op == ":=") %>% 
  mutate(est = ifelse(test = (pvalue > 0.05), yes = 0,
                      no = est)) %>% 
  mutate(ci.upper = ifelse(test = (pvalue > 0.05), yes = NA,
                           no = ci.upper)) %>%
  mutate(ci.lower = ifelse(test = (pvalue > 0.05), yes = NA,
                           no = ci.lower)) %>%
  separate(col = label, into = c("trait", "on", "type"), sep = "_") 

estimates2 <- left_join(estimates1, trait_key, by = "trait")
estimates2$loc <- rep("TLI", length(estimates2$trait))
estimates2$year <- rep("2017", length(estimates2$trait))

estimates_tli_17 <- estimates2


#### Step 6: Run the Model for TLI 2018 ####
# since this environment does not have yield per plant, 
# we remove paths in the initial model that have yield per plant
initial_model <- ' 
thousand_grain_weight ~ STD_TGW*stem_diameter + FLA_TGW*flag_leaf_area + ANT_TGW*anthesis_score + SPD_TGW*spikelet_density + HIT_TGW*height + EMP_TGW*emergence_percent + RTC_TGW*reproductive_tiller_ct
floret_site_utilization ~ STD_FSU*stem_diameter + FLA_FSU*flag_leaf_area + ANT_FSU*anthesis_score + SPD_FSU*spikelet_density + HIT_FSU*height + EMP_FSU*emergence_percent + RTC_FSU*reproductive_tiller_ct
yield_per_spike ~ FSU_YPS*floret_site_utilization + TGW_YPS*thousand_grain_weight + SPS_YPS*spikelets_per_spike + FPS_YPS*florets_per_spikelet 

# FPS
FPS_YPS_Direct := FPS_YPS

# ANT
ANT_YPS_Indirect := (ANT_FSU * FSU_YPS) + (ANT_TGW * TGW_YPS)

# STD
STD_YPS_Indirect := (STD_FSU * FSU_YPS) + (STD_TGW * TGW_YPS)

# SPD 
SPD_YPS_Indirect := (SPD_FSU * FSU_YPS) + (SPD_TGW * TGW_YPS)

# HIT 
HIT_YPS_Indirect := (HIT_FSU * FSU_YPS) + (HIT_TGW * TGW_YPS)

# FLA 
FLA_YPS_Indirect := (FLA_FSU * FSU_YPS) + (FLA_TGW * TGW_YPS)

# RTC
RTC_YPS_Indirect := (RTC_FSU * FSU_YPS) + (RTC_TGW * TGW_YPS)

# FSU
FSU_YPS_Direct := FSU_YPS

# TGW 
TGW_YPS_Direct := TGW_YPS

# SPS
SPS_YPS_Direct := SPS_YPS

# EMP
EMP_YPS_Indirect := (EMP_TGW * TGW_YPS) + (EMP_FSU * FSU_YPS)' 

# run the initial_model on TLI 2018 data
mod <- sem(initial_model, data=dat_tli18, fixed.x = FALSE)  # run model
lavaanPlot(model = mod, node_options = list(shape = "box"), stand = T, edge_options = list(color = "gray"), coef = T, covs = T, stars = T,  sig=0.05)
summary(mod, fit.measures = TRUE, standardized = T)
modificationIndices(mod, sort. = T)
inspect(mod, 'r2')

# results from attempt 1
# RMSEA 0.050
# CFI 0.989
# it's a good fit right away

# make a dataframe for ggplot
estimates <- as.data.frame(parameterEstimates(mod))

# export full path results
estimates_output <- estimates %>% 
  filter(op %in% c(":=", "~")) # filter out ~ paths and defined paths

write.csv(estimates_output, paste(path, "/estimates_tli18.csv", sep = ""))

estimates1 <- estimates %>% 
  filter(op == ":=") %>% 
  mutate(est = ifelse(test = (pvalue > 0.05), yes = 0,
                      no = est)) %>% 
  mutate(ci.upper = ifelse(test = (pvalue > 0.05), yes = NA,
                           no = ci.upper)) %>%
  mutate(ci.lower = ifelse(test = (pvalue > 0.05), yes = NA,
                           no = ci.lower)) %>%
  separate(col = label, into = c("trait", "on", "type"), sep = "_") 

estimates2 <- left_join(estimates1, trait_key, by = "trait")
estimates2$loc <- rep("TLI", length(estimates2$trait))
estimates2$year <- rep("2018", length(estimates2$trait))
estimates_tli_18 <- estimates2

#### Step 7: Prepare the Estimates for the Figure ####
# yield on a per plant basis was not measured at TLI in 2018, so these values need to be
# changed to zero for the sake of the plot, we'll borrow the values from TLI 2017

to_make_zero <- estimates_tli_17 %>% 
  filter(on == "YPP") %>%
  mutate(year = "2018", 
         est = 0, 
         ci.lower = 0, 
         ci.upper = 0)

estimates_tli_18_zero_fpy <- bind_rows(estimates_tli_18,to_make_zero)

# combine years within environments
all_loc_est_stp <- bind_rows(estimates_stp_17, estimates_stp_18)
all_loc_est_tli <- bind_rows(estimates_tli_17, estimates_tli_18_zero_fpy) # using the adjusted TLI18 dataset

#### Step 8: Create the Figure ####
# stacked ggplot to show direct and indirect effects
desired_order <- c("Stem Diameter", "Spikelet Density", "Flag Leaf Area", "Height", 
                      "Spike Emergence", "Anthesis Score", "Spikelets per Spike", "Florets per Spikelet", 
                      "Floret Site Utilization", "Thousand Grain Weight", "Yield per Spike", 
                      "Reproductive Tillers")

# since there are some cases where there are "stacked" plots, we have to adjust the position of the upper and lower
# confidence intervals to the SUM so they are in the appropriate positions
# this happens in two cases: STP 2017 in YPP and TLI 2017 in YPS


# fixing it for STP
all_loc_est_stp[24,11] <- all_loc_est_stp[24,11] + all_loc_est_stp[23,11]
all_loc_est_stp[24,12] <- all_loc_est_stp[24,12] + all_loc_est_stp[23,12]



stp <- ggplot() + 
  geom_bar(aes(y = est, x = pub_names_no_units, fill = type), color = "black", data = all_loc_est_stp, stat = "identity", position = position_stack(reverse =TRUE)) + 
  geom_errorbar(data = all_loc_est_stp, aes(x = pub_names_no_units, ymin = ci.lower, ymax = ci.upper),  width = 0.2)  + 
  facet_grid(cols = vars(on), rows = vars(year))  +
  coord_flip() +
  ggtitle("STP")+ 
  scale_x_discrete(limits = rev(desired_order)) + 
  scale_fill_manual(values=c("#0072B2", "#56B4E9")) + 
  scale_y_continuous(breaks = c(-.50, -0.25, 0, 0.25, 0.5, 0.75, 1), limits =c(-0.5, 1)) + 
  theme_bw() +
  xlab("Trait") + 
  ylab("") + 
  theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust =1 ), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.position = "bottom",  
        strip.text.x = element_text(size = 25, margin = margin(.4,0,.4,0, "cm")), 
        strip.background.x = element_rect(fill="white", color = "black"),
        strip.background.y = element_rect(fill="white"),
        strip.text.y = element_text(size = 25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0)



View(all_loc_est_tli)

# stacked ggplot to show direct and indirect effects
# manually adjust positions of error bars for TLI
all_loc_est_tli[6,11] <- all_loc_est_tli[6,11] + all_loc_est_tli[7,11]
all_loc_est_tli[6,12] <- all_loc_est_tli[6,12] + all_loc_est_tli[7,12]

tli <- ggplot() + 
  geom_bar(aes(y = est, x = pub_names_no_units, fill = type), color = "black", data = all_loc_est_tli, stat = "identity", position = position_stack(reverse =TRUE)) + 
  geom_errorbar(data = all_loc_est_tli, aes(x = pub_names_no_units, ymin = ci.lower, ymax = ci.upper),  width = 0.2, stat =  "identity")  + 
  facet_grid(cols = vars(on), rows = vars(year)) + 
  coord_flip() + 
  theme_bw() +
  ggtitle("TLI")+ 
  scale_x_discrete(limits = rev(desired_order)) + 
  scale_fill_manual(values=c("#C19417", "#f0cc65")) + 
  scale_y_continuous(breaks = c(-.50, -0.25, 0, 0.25, 0.5, 0.75, 1), limits =c(-0.5, 1)) + 
  xlab("") + 
  ylab("") +
  theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust =1 ), 
        axis.text.y = element_text(size = 20, color = "white"), 
        axis.title.x = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.position = "bottom",  
        strip.text.x = element_text(size = 25, margin = margin(.4,0,.4,0, "cm")), 
        strip.text.y = element_text(size = 25),
        strip.background.x = element_rect(colour = "black", fill=NA), 
        strip.background.y = element_rect(colour = "black", fill=NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0)

# using cowplot
plot <- plot_grid(stp, tli, ncol= 2)
# draw a single x axis title
ggdraw(add_sub(plot, "Standardized Path Coefficients", vpadding=grid::unit(0,"lines"),y=9, x=0.5, vjust=4.5, size = 20))
# export 1600 X 1000 tiff 




#### Step 9: Creating Nice Tables of Model Estimates ####
# read in data
stp17 <- read.csv(paste(path, "/estimates_stp17.csv", sep = ""), header = T) %>% mutate(loc = "STP",  year = "2017")
stp18 <- read.csv(paste(path, "/estimates_stp18.csv", sep = ""), header = T) %>% mutate(loc = "STP",  year = "2018")
tli17 <- read.csv(paste(path, "/estimates_tli17.csv", sep = ""), header = T) %>% mutate(loc = "TLI",  year = "2017")
tli18 <- read.csv(paste(path, "/estimates_tli18.csv", sep = ""), header = T) %>% mutate(loc = "TLI",  year = "2018")

all_dat <- rbind(stp17, stp18, tli17, tli18)

# create a key of trait names
key <- data.frame(trait  = c(trait, all_traits),  pub_name = rep(pub_names_no_units, 2))

# first take out different trait tiers
upper_tier <- all_dat %>% 
  filter(op == "~" & ! lhs %in% c("ten_spike_yield_per_spike_adj_thresh", "full_plant_plus_ten_spikes_yield_adj_thresh")) %>%
  dplyr::rename(trait = rhs, 
         yield = lhs) %>%
  mutate(effect = "Direct") %>%
  select(-X, -op, -label) %>%
  select(trait, yield, effect, est, se, z, pvalue, ci.lower, ci.upper, loc, year)


lower_tier <- all_dat %>% filter(op == ":=") %>%  # columns with this include both direct and indirect effects
  select(-X, -rhs, -label, -op) %>% # remove unnecesary columns
  separate(lhs, into = c("trait", "yield", "effect"), sep = "_") # separate out the first col

effects <- rbind(upper_tier, lower_tier)
effects <- left_join(effects, key, by = "trait")

colnames(key)[1] <- c("yield")
effects <- left_join(effects, key, by = "yield")

# format the table
effects <- effects %>% 
  select(pub_name.x, effect, pub_name.y, est, se, z, pvalue, ci.lower, ci.upper, loc, year) %>%
  dplyr::rename(Trait1 = pub_name.x, 
         Trait2 = pub_name.y, 
         Effect = effect) %>%
  mutate(est = round(est, 2), 
         se = round(se, 2), 
         z = round(z, 2), 
         pvalue = round(pvalue, 3), 
         ci.lower = round(ci.lower, 2), 
         ci.upper = round(ci.upper, 2))

# write out each table separately
# set vectors
loc <- c("STP", "TLI")
year <- c("2017", "2018")

for (i in 1:length(loc)) {
  for (j in 1:length(year)) {
    dat <- effects[effects$loc == loc[i],]
    dat1 <- dat[dat$year == year[j],]
    dat2 <- dat1 %>% select(-loc, -year) %>% arrange(Effect)
    write.csv(dat2, paste(path, "/paramter_estimates_table_", tolower(loc[i]), substr(year[j], 3, 4), ".csv", sep = ""), row.names = F)
  }
}
