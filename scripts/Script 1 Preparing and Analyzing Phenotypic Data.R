# Project: IWG_Yield_Components
# Script 1 - Preparing and Analyzing Phenotypic Data
# Author: Kayla R. Altendorf
# Date: 05/21/2020

# Required Pacakges:
library("dplyr") 
library("tidyr")
library("lme4")
library("lmerTest")
library("emmeans")
library("stringr")
library("reshape")
library("plyr")
library("multcompView")

# Declare where you want the output to go 
path <- c("/users/kayla.altendorf/OneDrive - USDA/Publications/IWG Yield Components/Scripts for GitHub/output/")
dir.create(path)


#### Step 1: Load Phenotypic Data #### 
# this is the same format that is available from the intermediate wheatgrass database
# using this dataset requires a bit of fanagalling, but for posterity it's best to have
# one version of the data in use
dat <- read.table("/Users/kayla.altendorf/OneDrive - USDA/Publications/IWG Yield Components/Scripts for GitHub/data/NAM_Data.txt", header = T, sep = "\t")

# select the traits that will be used for this analysis
traits <- c("FLFWD", "FLFLNG", "PTHT", "STMNUM", "STMDIA", "SPKDEN", "HDEMPER", "ZDK", 
            "FLORSPK", "SDSFLOR", "SPLYLD", "SPKYLD", "SDMG", "SPKHD")
dat1 <- dat %>% filter(trait_id %in% traits)

# I prefer to work with full length names, so I'll sub them in here
my_names <- c("flag_leaf_width", "flag_leaf_length", "height", "reproductive_tiller_ct", 
              "stem_diameter", "spikelet_density", "emergence_percent", "anthesis_score", "florets_per_spikelet", 
              "floret_site_utilization", "yield_per_plant", 
              "yield_per_spike", "thousand_grain_weight", "spikelets_per_spike")

trait_names <- data.frame(trait_id = traits, trait_id_full = my_names)
dat2 <- left_join(dat1, trait_names, by = "trait_id")


#### Step 2: Format Phenotypic Data #### 
# get rid of the unnecessary columns
# note: since spikelets per spike was taken three times, the samples need separate names so they can be averaged
dat3 <- dat2 %>% 
  dplyr::rename(year = phenotype_year, # rename cols to match personal preference
                famID = family_name, 
                col = range) %>%
  mutate(loc = substr(experiment_id, 4, 6), # extract location
         trait_id_full = case_when(trait_id_full == "spikelets_per_spike" ~  paste(trait_id_full, sample_number, sep = ""), 
                                   trait_id_full != "spikelets_per_spike" ~ paste(trait_id_full)), # give each subsample of spikelets_per_spike a unique name
         parent = substr(germplasm_id, 6, 6)) %>% # extract parent (C for common, d for donor, p for parent)
  filter(parent != "P") %>% # filter to exclude parents, as this project deals only with progeny
  select(famID, germplasm_id, loc, year, rep, trait_id_full, phenotype_value, plant_id) %>% 
  pivot_wider(names_from = trait_id_full, values_from = phenotype_value) %>% 
  select(-plant_id) %>% # pivot to wide format
  mutate(spikelets_per_spike1 = as.numeric(as.character(spikelets_per_spike1)), # make various traits as numeric
         spikelets_per_spike2 = as.numeric(as.character(spikelets_per_spike2)), 
         spikelets_per_spike3 = as.numeric(as.character(spikelets_per_spike3)),
         flag_leaf_length = as.numeric(as.character(flag_leaf_length)), 
         flag_leaf_width = as.numeric(as.character(flag_leaf_width)), 
         loc = str_replace(loc, "SAL", "TLI")) # replace SAL (Salina) with TLI

# take the average for spikelets and then create a column for merging
# and calculate flag leaf area (length * width)
dat4 <- dat3 %>% mutate(spikelets_per_spike = rowMeans(dat3[,10:12], na.rm = T)) %>% # calculate the average spikelets per spike
  select(-spikelets_per_spike1, -spikelets_per_spike2, -spikelets_per_spike3) %>%
  mutate(merge_col = paste(germplasm_id, loc, year, rep, sep = "_"), 
         flag_leaf_area = (flag_leaf_width*0.1) * flag_leaf_length) %>% # convert width mm to cm 
  select(-famID, -germplasm_id, -loc, -year, -rep, -flag_leaf_width, -flag_leaf_length)


# change feekes scores to coded values for easier analysis
feekes <- data.frame(anthesis_score = c("49", "51", "53", "55", "57", "59", "61", "65", "69", "71"), coded = 1:10)
dat5 <- left_join(dat4, feekes, by = "anthesis_score") %>% 
  select(-anthesis_score) %>%
  dplyr::rename(anthesis_score = coded)

# the database includes *all* entries, including parents, plants that were identified later as selfs, and so on.
# futhermore, the population itself is not balanced (e.g. unequal numbers of individuals within families, 
# entries per location and year), which causes problems with the ANOVA
# to address this, we will load in a backbone dataset, which is balanced using NA values,
# and we'll format the data to match. 

# read in the backbone csv
backbone <- read.csv("/Users/kayla.altendorf/OneDrive - USDA/Publications/IWG Yield Components/Scripts for GitHub/data/backbone.csv")
# plantID3 is the balanced plantID
# example:
backbone %>% group_by(loc, year, rep, famID) %>% tally() # all have 133 entries

# left join with the data
dat6 <- left_join(backbone, dat5, by = "merge_col") %>% select(-merge_col)
dat6 %>% group_by(loc, year, rep, famID) %>% tally() # make sure it's still balanced

# change all selfs to NA
for (i in 1:nrow(dat6)) {
  if (! is.na(dat6$self[i])) { # if self is not NA (e.g. if it is outcross or self)
    dat6[i, 13:25] <- NA # change all phneotype data columns to NA
  }
}

# write out the final dataset
write.csv(dat6, paste(path, "data.csv", sep = "/"), row.names = F)

# read it back in to convert everything to numerics
dat <- read.csv(paste(path, "data.csv", sep = "/"))
dat$year <- as.factor(dat$year)
dat$rep <- as.factor(dat$rep)

#### Step 3: Analysis of Variance for Combined Analysis #### 
traits <- colnames(dat)[13:25] # grab the trait names from column headers

# make data as factor
dat$rep <- as.factor(dat$rep)
dat$year <- as.factor(dat$year)
dat$plantID3 <- as.factor(dat$plantID3)
dat$famID <- as.factor(dat$famID)

# run through each trait
traits

# height
model <- lmer(height ~ famID * loc * year + (1|loc:rep) + (1|loc:famID:plantID3), data = dat)
anova <- anova(model, type="II")
plot(model)
dir.create(paste(path, "height", sep = "/"))
emmeans_loc_year <- CLD(emmeans(model, ~  year * loc), Letters = c(LETTERS))
write.table(anova, paste(path, "/height/height_anova.txt", sep = ""), quote = F, row.names = T, col.names = NA, sep = "\t")
write.table(emmeans_loc_year, paste(path, "/height/height_emmeans_loc_year.txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

# reproductive_tiller_ct
model <- lmer(sqrt(reproductive_tiller_ct) ~ famID * loc * year + (1|loc:rep) + (1|loc:famID:plantID3), data = dat)
anova <- anova(model)
plot(model)
dir.create(paste(path, "reproductive_tiller_ct", sep = "/"))
emmeans_loc_year <- CLD(emmeans(model, ~  year * loc),  type = "response", Letters = c(LETTERS))
write.table(anova, paste(path, "/reproductive_tiller_ct/reproductive_tiller_ct_anova.txt", sep = ""), quote = F, row.names = T, col.names = NA, sep = "\t")
write.table(emmeans_loc_year, paste(path, "/reproductive_tiller_ct/reproductive_tiller_ct_emmeans_loc_year.txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

# florets_per_spikelet
model <- lmer(florets_per_spikelet ~ famID * loc * year + (1|loc:rep) + (1|loc:famID:plantID3), data = dat)
anova <- anova(model, type="II")
plot(model)
dir.create(paste(path, "florets_per_spikelet", sep = "/"))
emmeans_loc_year <- CLD(emmeans(model, ~  year * loc), Letters = c(LETTERS))
write.table(anova, paste(path, "/florets_per_spikelet/florets_per_spikelet_anova.txt", sep = ""), quote = F, row.names = T, col.names = NA, sep = "\t")
write.table(emmeans_loc_year, paste(path, "/florets_per_spikelet/florets_per_spikelet_emmeans_loc_year.txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

# emergence_percent
model <- lmer(emergence_percent ~ famID * loc * year + (1|loc:rep) + (1|loc:famID:plantID3), data = dat)
anova <- anova(model, type="II")
plot(model)
dir.create(paste(path, "emergence_percent", sep = "/"))
emmeans_loc_year <- CLD(emmeans(model, ~  year * loc), Letters = c(LETTERS))
write.table(anova, paste(path, "/emergence_percent/emergence_percent_anova.txt", sep = ""), quote = F, row.names = T, col.names = NA, sep = "\t")
write.table(emmeans_loc_year, paste(path, "/emergence_percent/emergence_percent_emmeans_loc_year.txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

# stem_diameter
model <- lmer(stem_diameter ~ famID * loc * year + (1|loc:rep) + (1|loc:famID:plantID3), data = dat)
anova <- anova(model, type="II")
plot(model)
dir.create(paste(path, "stem_diameter", sep = "/"))
emmeans_loc_year <- CLD(emmeans(model, ~  year * loc), Letters = c(LETTERS))
write.table(anova, paste(path, "/stem_diameter/stem_diameter_anova.txt", sep = ""), quote = F, row.names = T, col.names = NA, sep = "\t")
write.table(emmeans_loc_year, paste(path, "/stem_diameter/stem_diameter_emmeans_loc_year.txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

# yield_per_plant 
dat.yield_per_plant<- dat %>% # making a new data frame to account for missing data
  mutate(loc_rep=paste(loc,"_",rep, sep = ""))
model <- lmer(sqrt(yield_per_plant) ~ famID * loc + (1|plantID3) + (1|loc:plantID3), data = filter(dat.yield_per_plant, year==2017 & loc_rep == c("STP_1", "TLI_2")))
anova <- anova(model, type="II")
plot(model)
dir.create(paste(path, "yield_per_plant", sep = "/"))
emmeans_loc_year <- CLD(emmeans(model, ~  loc), type = "response", Letters = c(LETTERS))
write.table(anova, paste(path, "/yield_per_plant/yield_per_plant_anova.txt", sep = ""), quote = F, row.names = T, col.names = NA, sep = "\t")
write.table(emmeans_loc_year, paste(path, "/yield_per_plant/yield_per_plant_emmeans_loc_year.txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

# floret_site_utilization
model <-lmer(sqrt(floret_site_utilization) ~ famID * loc * year + (1|loc:rep) + (1|loc:famID:plantID3), data = dat)
anova <- anova(model, type="II")
plot(model)
dir.create(paste(path, "floret_site_utilization", sep = "/"))
emmeans_loc_year <- CLD(emmeans(model, ~  year * loc),  type = "response", Letters = c(LETTERS))
write.table(anova, paste(path, "/floret_site_utilization/floret_site_utilization_anova.txt", sep = ""), quote = F, row.names = T, col.names = NA, sep = "\t")
write.table(emmeans_loc_year, paste(path, "/floret_site_utilization/floret_site_utilization_emmeans_loc_year.txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

# spikelet_density
model <-lmer(sqrt(spikelet_density) ~ famID * loc * year + (1|loc:rep) + (1|loc:famID:plantID3), data = dat)
anova <- anova(model, type="II")
plot(model)
dir.create(paste(path, "spikelet_density", sep = "/"))
emmeans_loc_year <- CLD(emmeans(model, ~  year * loc),  type = "response", Letters = c(LETTERS))
write.table(anova, paste(path, "/spikelet_density/spikelet_density_anova.txt", sep = ""), quote = F, row.names = T, col.names = NA, sep = "\t")
write.table(emmeans_loc_year, paste(path, "/spikelet_density/spikelet_density_emmeans_loc_year.txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

# yield_per_spike
model <-lmer(sqrt(yield_per_spike) ~ famID * loc * year + (1|loc:rep) + (1|loc:famID:plantID3), data = dat)
anova <- anova(model, type="II")
plot(model)
dir.create(paste(path, "yield_per_spike", sep = "/"))
emmeans_loc_year <- CLD(emmeans(model, ~  year * loc),  type = "response", Letters = c(LETTERS))
write.table(anova, paste(path, "/yield_per_spike/yield_per_spike_anova.txt", sep = ""), quote = F, row.names = T, col.names = NA, sep = "\t")
write.table(emmeans_loc_year, paste(path, "/yield_per_spike/yield_per_spike_emmeans_loc_year.txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

# thousand_grain_weight
model <- lmer(thousand_grain_weight ~ famID * loc * year + (1|loc:rep) + (1|loc:famID:plantID3), data = dat)
anova <- anova(model, type="II")
plot(model)
dir.create(paste(path, "thousand_grain_weight", sep = "/"))
emmeans_loc_year <- CLD(emmeans(model, ~  year * loc), Letters = c(LETTERS))
write.table(anova, paste(path, "/thousand_grain_weight/thousand_grain_weight_anova.txt", sep = ""), quote = F, row.names = T, col.names = NA, sep = "\t")
write.table(emmeans_loc_year, paste(path, "/thousand_grain_weight/thousand_grain_weight_emmeans_loc_year.txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

# spikelets_per_spike
model <- lmer(spikelets_per_spike ~ famID * loc * year + (1|loc:rep) + (1|loc:famID:plantID3), data = dat)
anova <- anova(model, type="II")
plot(model)
dir.create(paste(path, "spikelets_per_spike", sep = "/"))
emmeans_loc_year <- CLD(emmeans(model, ~  year * loc), Letters = c(LETTERS))
write.table(anova, paste(path, "/spikelets_per_spike/spikelets_per_spike_anova.txt", sep = ""), quote = F, row.names = T, col.names = NA, sep = "\t")
write.table(emmeans_loc_year, paste(path, "/spikelets_per_spike/spikelets_per_spike_emmeans_loc_year.txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

# flag_leaf_area
model <-lmer(log(flag_leaf_area) ~ famID * loc * year + (1|loc:rep) + (1|loc:famID:plantID3), data = dat)
anova <- anova(model, type="II")
plot(model)
dir.create(paste(path, "flag_leaf_area", sep = "/"))
emmeans_loc_year <- CLD(emmeans(model, ~  year * loc),  type = "response", Letters = c(LETTERS))
write.table(anova, paste(path, "/flag_leaf_area/flag_leaf_area_anova.txt", sep = ""), quote = F, row.names = T, col.names = NA, sep = "\t")
write.table(emmeans_loc_year, paste(path, "/flag_leaf_area/flag_leaf_area_emmeans_loc_year.txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

# anthesis_score
model <- lmer(anthesis_score ~ famID * loc * year + (1|loc:rep) + (1|loc:famID:plantID3), data = dat)
anova <- anova(model, type="II")
plot(model)
dir.create(paste(path, "anthesis_score", sep = "/"))
emmeans_loc_year <- CLD(emmeans(model, ~  year * loc), Letters = c(LETTERS))
write.table(anova, paste(path, "/anthesis_score/anthesis_score_anova.txt", sep = ""), quote = F, row.names = T, col.names = NA, sep = "\t")
write.table(emmeans_loc_year, paste(path, "/anthesis_score/anthesis_score_emmeans_loc_year.txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")



#### Step 4: Analysis of Variance within Environments #### 

# filter data by location
stp17 <- filter(dat, loc == "STP" & year == "2017")
stp18 <- filter(dat, loc == "STP" & year == "2018")
tli17 <- filter(dat, loc == "TLI" & year == "2017")
tli18 <- filter(dat, loc == "TLI" & year == "2018")

loc_list <- list(stp17, stp18, tli17, tli18)
loc_names <- c("stp17", "stp18", "tli17", "tli18")

# remove reproductive tiller number and full plant yield, since that trait is missing tli18 environment
normal_traits <- traits[-c(2, 6)]

for (j in 1:length(loc_list)) {
  for (i in 1:length(normal_traits)) {
    formula <- paste0(normal_traits[i], " ~ famID/plantID3 + rep", sep = "")
    model <- lm(formula, data = loc_list[[j]])
    an <- as.data.frame(anova(model))
    emmeans_fam <- as.data.frame(CLD(emmeans(model, ~ famID), Letters = c(LETTERS)))
    emmeans_genet <- as.data.frame(emmeans(model, ~ plantID3|famID))
    write.table(an, paste(path, "/", normal_traits[i], "/", normal_traits[i], "_anova_", loc_names[j], ".txt", sep = ""), quote = F, row.names = T, col.names = NA, sep = "\t")
    write.table(emmeans_fam, paste(path, "/", normal_traits[i], "/", normal_traits[i], "_emmeans_fam_", loc_names[j], ".txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")
    write.table(emmeans_genet, paste(path, "/", normal_traits[i], "/", normal_traits[i], "_emmeans_genet_", loc_names[j], ".txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")
    print(paste(Sys.time(), "done with", normal_traits[i], sep = " "))
  }
}

# for reproductive tiller ct
for (j in 1:length(loc_list)) { 
  model <- lm(sqrt(reproductive_tiller_ct) ~ famID/plantID3 + rep, data = loc_list[[j]])
  anova <- anova(model)
  emmeans_fam <- as.data.frame(CLD(emmeans(model, ~ famID, type = "response"), Letters = c(LETTERS)))
  emmeans_genet <- as.data.frame(emmeans(model,  ~ plantID3|famID, type = "response"))
  emmeans_genet_transformed_scale <- as.data.frame(emmeans(model,  ~ plantID3|famID))
  write.table(anova, paste(path, "/", traits[2], "/", traits[2], "_anova_", loc_names[j], ".txt", sep = ""), quote = F, row.names = T, col.names = NA, sep = "\t")
  write.table(emmeans_fam, paste(path, "/", traits[2], "/", traits[2], "_emmeans_fam_", loc_names[j], ".txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")
  write.table(emmeans_genet, paste(path, "/", traits[2], "/", traits[2], "_emmeans_genet_", loc_names[j], ".txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")
  write.table(emmeans_genet_transformed_scale, paste(path, "/", traits[2], "/", traits[2], "_emmeans_genet_transformed_scale_", loc_names[j], ".txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")
}

# for yield per plant, not measured in all four envs
loc_list_ypp <- loc_list[-4]
loc_names_ypp <- loc_names[-4]

for (j in 1:length(loc_list_ypp)) { 
  formula <- paste(traits[6], "~ famID/plantID3 + rep", sep ="")
  model <- lm(formula, data = loc_list_ypp[[j]])
  anova <- anova(model)
  emmeans_fam <- as.data.frame(CLD(emmeans(model, ~ famID, Letters = c(LETTERS))))
  emmeans_genet <- as.data.frame(emmeans(model,  ~ plantID3|famID))
  write.table(anova, paste(path, "/", traits[6], "/", traits[6], "_anova_", loc_names_ypp[j], ".txt", sep = ""), quote = F, row.names = T, col.names = NA, sep = "\t")
  write.table(emmeans_fam, paste(path, "/", traits[6], "/", traits[6], "_emmeans_fam_", loc_names_ypp[j], ".txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")
  write.table(emmeans_genet, paste(path, "/", traits[6], "/", traits[6], "_emmeans_genet_", loc_names_ypp[j], ".txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")
}
