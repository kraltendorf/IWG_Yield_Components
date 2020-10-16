# Project: IWG_Yield_Components
# Script 2 - Calculating Heritability and Correlations Between Traits
# Author: Kayla R. Altendorf
# Date: 05/25/2020

# Required Packages:
library("dplyr") 
library("tidyr")
library("tibble")
library("lme4")
library("emmeans")
library("stringr")
library("reshape")
library("plyr")
library("multcompView")
library("Hmisc")
library("gtools")
library("cowplot")

# Declare where you want the output to go 
path <- c("/users/kayla.altendorf/OneDrive - USDA/Publications/IWG Yield Components/Scripts for GitHub/output/")

# Get directory (trait) names
dirs <- list.dirs(path = path, full.names = TRUE, recursive = FALSE)

# Get the id_frame from data.csv
dat <- read.csv(paste(path, "/data.csv", sep = ""), header = T)
id_frame <- dat %>% 
  select(famID, plantID3, longID) %>% 
  mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% 
  distinct() %>%
  filter(! is.na(longID))

# Set locations in order
env <- c("stp17", "stp18", "tli17", "tli18")

# Set pub names, or versions of trait names that will go into the publication
# In the same order as the traits in the dataframe
pub_names <- c("Anthesis (1-10)", "Spike Emergence (%)", "Flag Leaf Area (cm2)", "Floret Site Utilization (%)", 
               "Florets Spikelet-1 (ct)", "Height (cm)", "Reproductive Tillers (ct)", "Spikelet Density (cm)",
               "Spikelets Spike-1 (ct)", "Stem Diameter (mm)", "Thousand Grain Weight (mg)",  "Yield Plant-1 (g)", 
               "Yield Spike-1 (g)")



#### Step 1: Acquire Emmeans from Script 1 Output ####

# create empty dataframe
emmeans_env <- data.frame(matrix(NA, nrow = 1295, ncol = 13))
emmeans_all <- list()

# since tli18 yield per plant is missing, we will create a filler using tli17 of NAs so the code runs fine
filler <- read.table(paste(path, "/yield_per_plant/yield_per_plant_emmeans_genet_tli17.txt", sep = ""), head = T)
filler$emmean <- NA
write.table(filler, paste(path, "/yield_per_plant/yield_per_plant_emmeans_genet_tli18.txt", sep = ""),  quote = F, row.names = F, col.names = T, sep = "\t")

for (j in 1:length(env)) {
  for (i in 1:length(dirs)) {
    file <- list.files(dirs[i], pattern = paste("emmeans_genet_", env[j], sep = ""), full.names = TRUE)
    emmeans <- read.table(file, head = T)
    emmeans1 <- emmeans %>% mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% dplyr::select(8, 3)
    emmeans2 <- left_join(id_frame, emmeans1, by = "famID_plantID3")
    emmeans_env[,i] <- emmeans2[,5] 
    colnames(emmeans_env)[i] <- str_split(dirs[i], "/")[[1]][8] 
    }
  emmeans_all[[j]] <- emmeans_env
}


#### Step 2: Run Correlation Analysis ####
cor_out <- list()

for (i in 1:length(emmeans_all)){
  colnames(emmeans_all[[i]]) <- pub_names
  cor <- rcorr(as.matrix(emmeans_all[[i]]), type = "pearson")
  cor.p <- as.data.frame(cor$P)
  cor.r <- as.data.frame(cor$r)
  cor.r <- cor.r %>% mutate_all(funs(round(., 2))) # round r to two digits
  cor.p.stars <- cor.p
  for (k in 1:ncol(cor.p)) {
    for (j in 1:nrow(cor.p)) {
      cor.p.stars[k,j] <- stars.pval(cor.p[k,j])
    }
  }
  cor.r.p <- cor.r
  for (k in 1:ncol(cor.p)) {
    for (j in 1:nrow(cor.r)) {
      cor.r.p[k,j] <- paste(cor.r[k,j], cor.p.stars[k,j], sep = "")
    }
  }
  cor.r.p[upper.tri(cor.r.p)] <- NA
  row.names(cor.r.p) <- pub_names
  cor_out[[i]] <- cor.r.p
}

# create output directory
dir.create(paste(path, "/trait_correlations", sep = ""))

# write out output
for (j in 1:length(env)) {
  write.table(cor_out[[j]], paste(path, "/trait_correlations/cor_sig", env[j], ".txt", sep = ""), quote = F, row.names = T, col.names = NA, sep = "\t")
}



#### Step 3: Calculate Broad and Narrow Sense Heritability ####

# create a dataframe for the output
herit_df <- data.frame(env = c("STP", "TLI"), broad = NA, narrow = NA, trait = NA)
herit_df_list <- replicate(13, herit_df, simplify = FALSE)

# extract trait names
all_traits <- list()
for (i in 1:length(dirs)) {
  all_traits[[i]] <- str_split(dirs[i], "/")[[1]][10]
}

all_traits <- unlist(all_traits)

# set vectors for iterating through years and locations
loc <- c("STP", "TLI")

# make important terms as factor
dat$year <- as.factor(dat$year)
dat$rep <- as.factor(dat$rep)

for (i in 1:length(all_traits)) {
  for (j in 1:length(loc)) {

    # extract data 
    dat_loc_year <- dat[dat$loc == loc[j],]
    
    # calculate broad sense on a genet mean basis
    # if the trait is reproductive tiller number, it requires a transformation, if not, proceed with regular equation
    if (all_traits[i] == "reproductive_tiller_ct") {
      formula <- paste0("sqrt(reproductive_tiller_ct)", " ~ (1|famID:plantID3) + (1|year) + (1|famID:plantID3:year) + (1|rep) + (1|famID:plantID3:rep)" , sep = "")}
    if (! all_traits[i] == "reproductive_tiller_ct") {formula <- paste0(all_traits[i], " ~ (1|famID:plantID3) + (1|year) + (1|famID:plantID3:year) + (1|rep) + (1|famID:plantID3:rep)", sep = "")}
    
    # if we're onto yield per plant at tli in 2018, change values to NA and skip everything else, or it fails
    if (all_traits[i] == "yield_per_plant") {
      heritability_broad <- NA
    }
    
    else {
      model <- lmer(formula, data = dat_loc_year)
      
      
      output <- as.data.frame(VarCorr(model))
      vg <- output$vcov[3] # extracting appropriate variance components
      vgy <- output$vcov[1]
      vgr <- output$vcov[6]
      heritability_broad <- vg / ((vg) + (vgy / 2) + (vgr / 4))
    
      # clean out variables before next iteration
      model <- NA
      output <- NA
      vg <- NA
      vgy <- NA
      vgr <- NA
    
      # calculate narrow sense according to falconer, again if reproductive tiller, requires trans
      # in this case, family and rep are random
      if (all_traits[i] == "reproductive_tiller_ct") {
        formula <- paste0("sqrt(reproductive_tiller_ct)", " ~ (1|famID) + (1|year) + (1|famID:year) + (1|rep) + (1|famID:rep)" , sep = "")}
      if (all_traits[i] == "yield_per_plant") {
        formula <- paste0("sqrt(reproductive_tiller_ct)", " ~ (1|famID) + (1|year) + (1|famID:year)" , sep = "")} # remove
      
      }
      if (! all_traits[i] == "reproductive_tiller_ct") {formula <- paste0(all_traits[i], "~ (1|famID) + (1|year) + (1|famID:year) + (1|rep) + (1|famID:rep)" , sep = "")}
    if (all_traits[i] == "yield_per_plant") {
      heritability_narrow <- NA
    }
    else {
      model <- lmer(formula, data = dat_loc_year)
      output <- as.data.frame(VarCorr(model))
      vf <- output$vcov[3] # extracting appropriate variance components
      ve <- output$vcov[6]
      heritability_narrow <- (vf * 4) / ((vf * 4) + (ve))
    
      # clean out variables before next iteration
      model <- NA
      output <- NA
      vf <- NA
      ve <- NA
  
    }
      # output result into heritability dataframe
      herit_df_list[[i]][j,2] <- heritability_broad
      herit_df_list[[i]][j,3] <- heritability_narrow
      herit_df_list[[i]][j,4] <- all_traits[i]
  }
}


#### now calculate yield_per_plant separately ####
dat_loc_year <- dat[dat$loc == "STP",]


# for STP
# broad sense
formula <- yield_per_plant ~ (1|famID:plantID3) + (1|year) + (1|famID:plantID3:year) + (1|rep) + (1|famID:plantID3:rep)
model <- lmer(formula, data = dat_loc_year)
output <- as.data.frame(VarCorr(model))
vg <- output$vcov[3] # extracting appropriate variance components
vgy <- output$vcov[1]
vgr <- output$vcov[6]
heritability_broad <- vg / ((vg) + (vgy / 2) + (vgr / 4))


# narrow
formula <-  yield_per_plant ~ (1|famID) + (1|year) + (1|famID:year) + (1|rep) + (1|famID:rep)
model <- lmer(formula, data = dat_loc_year)
output <- as.data.frame(VarCorr(model))
vf <- output$vcov[3] # extracting appropriate variance components
ve <- output$vcov[6]
heritability_narrow <- (vf * 4) / ((vf * 4) + (ve))

# edit the STP portion of this dataframe within the list
herit_df_list[[12]][1,2] <- heritability_broad
herit_df_list[[12]][1,3] <- heritability_narrow

  
  

# format the output into a nice table with means
# and sub in pub_names

herit_df <- do.call("rbind", herit_df_list)
herit_df_wide <- herit_df %>% 
  pivot_wider(values_from = c("broad", "narrow"), names_from = c("env")) %>%
  select(trait, broad_STP, narrow_STP, broad_TLI, narrow_TLI) %>%
  mutate(trait = pub_names)

# calculate and append averages
herit_df_wide_avg <- herit_df_wide %>% 
  mutate_if(is.numeric, funs(round(., 2)))

# write out result
write.table(herit_df_wide_avg, paste(path, "/heritabilities.txt", sep = ""),  quote = F, row.names = F, sep = "\t")




#### Step 4: Make a Phenotypic Data Summary Table ####
trait_summary <- list()

for (i in 1:length(emmeans_all)){
  trait_min <- as.data.frame(sapply(emmeans_all[[i]], min, na.rm = T))
  trait_max <- as.data.frame(sapply(emmeans_all[[i]], max, na.rm = T))
  trait_mean <- as.data.frame(sapply(emmeans_all[[i]], mean, na.rm = T))
  trait_sd <- as.data.frame(sapply(emmeans_all[[i]], function(x)sd(x, na.rm = T)))
  trait_sum <- cbind(trait_min, trait_max[,1], trait_mean[,1], trait_sd[,1])
  trait_sum1 <- rownames_to_column(trait_sum, "trait")
  trait_sum2 <- trait_sum1 %>% mutate_if(is.numeric, funs(round(., 2)))
  trait_summary[[i]] <- trait_sum2
}

# make TLI 2018 full plant yield NAs
trait_summary[[4]][12,2:5] <- c(NA, NA, NA, NA)

# rbind all environments together and rename column headers
trait_summary <- cbind(trait_summary[[1]][,1], trait_summary[[1]][,-1], trait_summary[[2]][,-1], trait_summary[[3]][,-1], trait_summary[[4]][,-1])
colnames(trait_summary) <- c("trait", "min", "max", "mean", "sd", "min", "max", "mean", "sd", "min", "max", "mean", "sd", "min", "max", "mean", "sd")

# write output
write.table(trait_summary, paste(path, "/trait_summary.txt", sep = ""), quote = F, row.names = F, sep = "\t")




#### Step 5: Make Boxplots of Pheotypic Data ####
loc <- c("STP", "STP", "TLI", "TLI")
year <- c("2017", "2018", "2017", "2018")
# add year and location information to the emmeans

for (j in 1:length(emmeans_all)) {
  colnames(emmeans_all[[j]]) <- all_traits # return to non pub names
  emmeans_all[[j]]$loc <- loc[j]
  emmeans_all[[j]]$year <- year[j]
}

# rbind into one dataframe
emmeans <- do.call("rbind", emmeans_all)

# we'll go through and create a boxplot for each trait individually
# hardcoding in the average heritabilities as we calculated previously
all_traits
colnames(emmeans)

# anthesis_score
anthesis_score <- ggplot(emmeans, aes(x = loc, y = anthesis_score)) + 
  geom_boxplot(aes(fill=interaction(loc, year))) +
  scale_fill_manual(values=c("#0072B2", "#C19417", "#56B4E9", "#fce6a4")) +
  theme_bw() +
  ylab(pub_names[1]) + 
  theme(legend.position="none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank()) +
  ylim(NA, 10.5)  +
  annotate("text",  x=Inf, y = Inf, label = "H = 0.82; 0.78", vjust= 1.5, hjust=1.05) +
  annotate("text",  x=Inf, y = Inf, label = expression(paste(italic(h)^2, " = 0.63; 0.53")), vjust= 2.5, hjust=1.05)

# emergence_percent
emergence_percent <- ggplot(emmeans, aes(x = loc, y = emergence_percent)) + 
  geom_boxplot(aes(fill=interaction(loc, year))) +
  scale_fill_manual(values=c("#0072B2", "#C19417", "#56B4E9", "#fce6a4")) +
  theme_bw() +
  ylab(pub_names[2]) +
  theme(legend.position="none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank()) + 
  annotate("text",  x=Inf, y = Inf, label = "H = 0.78; 0.72", vjust= 1.5, hjust= 3.3) +
  annotate("text",  x=Inf, y = Inf, label = expression(paste(italic(h)^2, " = 0.52, 0.32")), vjust= 2.5, hjust=3.1)

# flag_leaf_area
flag_leaf_area <- ggplot(emmeans, aes(x = loc, y = flag_leaf_area)) + 
  geom_boxplot(aes(fill=interaction(loc, year))) +
  scale_fill_manual(values=c("#0072B2", "#C19417", "#56B4E9", "#fce6a4")) +
  theme_bw() +
  scale_y_continuous(breaks=c(0, 25, 50, 75, 100))+
  #ylab(pub_names[3]) +
  ylab(expression('Flag Leaf Area (cm'^2*')')) +
  theme(legend.position="none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank()) +
  annotate("text",  x=Inf, y = Inf, label = "H = 0.57; 0.44", vjust= 1.5, hjust=1.05) +
  annotate("text",  x=Inf, y = Inf, label = expression(paste(italic(h)^2, " = 0.13; 0.18")), vjust= 2.5, hjust=1.05)

# floret_site_utilization
floret_site_utilization <- ggplot(emmeans, aes(x = loc, y = floret_site_utilization)) + 
  geom_boxplot(aes(fill=interaction(loc, year))) +
  scale_fill_manual(values=c("#0072B2", "#C19417", "#56B4E9", "#fce6a4")) +
  scale_y_continuous(breaks=c(0, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50)) + 
  theme_bw() +
  ylab(pub_names[4]) +
  theme(legend.position="none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank()) +
  annotate("text",  x=Inf, y = Inf, label = "H = 0.61; 0.50", vjust= 1.5, hjust=1.05) +
  annotate("text",  x=Inf, y = Inf, label = expression(paste(italic(h)^2, " = 0.12; 0.00")), vjust= 2.5, hjust=1.05)

# florets_per_spikelet
florets_per_spikelet <- ggplot(emmeans, aes(x = loc, y = florets_per_spikelet)) + 
  geom_boxplot(aes(fill=interaction(loc, year))) +
  scale_fill_manual(values=c("#0072B2", "#C19417", "#56B4E9", "#fce6a4")) +
  theme_bw() +
  ylab(expression('Florets Spikelet'^-1*'(ct)')) +
  #ylab(pub_names[5]) + 
  theme(legend.position="none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank()) +
  annotate("text",  x=Inf, y = Inf, label = "H = 0.55", vjust= 1.5, hjust=1.05) +
  annotate("text",  x=Inf, y = Inf, label = expression(paste(italic(h)^2, " = 0.52")), vjust= 2.5, hjust=1.05)

# height
height <- ggplot(emmeans, aes(x = loc, y = height)) + 
  geom_boxplot(aes(fill=interaction(loc, year))) +
  scale_fill_manual(values=c("#0072B2", "#C19417", "#56B4E9", "#fce6a4")) +
  theme_bw() +
  scale_y_continuous(breaks=c(25, 50, 75, 100, 125, 150)) + 
  ylab(pub_names[6]) +
  theme(legend.position="none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank()) +
  annotate("text",  x=Inf, y = Inf, label = "H = 0.85; 0.71", vjust= 1.5, hjust=1.05) +
  annotate("text",  x=Inf, y = Inf, label = expression(paste(italic(h)^2, " = 0.63; 0.46")), vjust= 2.5, hjust=1.05)

# reproductive_tiller_ct
reproductive_tiller_ct <- ggplot(emmeans, aes(x = loc, y = reproductive_tiller_ct)) + 
  geom_boxplot(aes(fill=interaction(loc, year))) + 
  scale_fill_manual(values=c("#0072B2", "#C19417", "#56B4E9", "#fce6a4")) +
  theme_bw() +
  ylab(pub_names[7]) +
  theme(legend.position="none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank()) +
  annotate("text",  x=Inf, y = Inf, label = "H = 0.88; 0.47", vjust= 1.5, hjust=1.05) +
  annotate("text",  x=Inf, y = Inf, label = expression(paste(italic(h)^2, " = 0.54; 0.19")), vjust= 2.5, hjust=1.05) +
  annotate("text", x = c(0.8, 1.2, 1.8, 2.2), y = c(85, 185, 110, 95), label = c("bc", "a", "b", "c"))

# spikelet_density
spikelet_density <- ggplot(emmeans, aes(x = loc, y = spikelet_density)) + 
  geom_boxplot(aes(fill=interaction(loc, year))) +
  scale_fill_manual(values=c("#0072B2", "#C19417", "#56B4E9", "#fce6a4")) +
  scale_y_continuous(breaks=c(.50, .75, 1.00, 1.25, 1.5)) + 
  theme_bw() +
  ylab(pub_names[8]) +
  theme(legend.position="none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank()) +
  annotate("text",  x=Inf, y = Inf, label = "H = 0.85; 0.74", vjust= 1.5, hjust=3.2) +
  annotate("text",  x=Inf, y = Inf, label = expression(paste(italic(h)^2, " = 0.66; 0.47")), vjust= 2.5, hjust=3.05)

# spikelets_per_spike
spikelets_per_spike <- ggplot(emmeans, aes(x = loc, y = spikelets_per_spike)) + 
  geom_boxplot(aes(fill=interaction(loc, year))) +
  scale_fill_manual(values=c("#0072B2", "#C19417", "#56B4E9", "#fce6a4")) +
  theme_bw() +
  #ylab(pub_names[9]) +
  ylab(expression('Spikelets Spike'^-1*'(ct)')) +
  theme(legend.position="none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank()) +
  annotate("text",  x=Inf, y = Inf, label = "H = 0.72; 0.65", vjust= 1.5, hjust=3.1) +
  annotate("text",  x=Inf, y = Inf, label = expression(paste(italic(h)^2, " = 0.44; 0.40")), vjust= 2.5, hjust = 3)

# stem_diameter
stem_diameter <- ggplot(emmeans, aes(x = loc, y = stem_diameter)) + 
  geom_boxplot(aes(fill=interaction(loc, year))) +
  scale_fill_manual(values=c("#0072B2", "#C19417", "#56B4E9", "#fce6a4")) +
  theme_bw() +
  ylab(pub_names[10]) +
  theme(legend.position="none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank()) + 
  annotate("text",  x=Inf, y = Inf, label = "H = 0.63; 0.65", vjust= 1.5, hjust=1.05) +
  annotate("text",  x=Inf, y = Inf, label = expression(paste(italic(h)^2, " = 0.27; 0.22")), vjust= 2.5, hjust=1.05)

# thousand_grain_weight
thousand_grain_weight <- ggplot(emmeans, aes(x = loc, y = thousand_grain_weight)) + 
  geom_boxplot(aes(fill=interaction(loc, year))) +
  scale_fill_manual(values=c("#0072B2", "#C19417", "#56B4E9", "#fce6a4")) +
  theme_bw() +
  ylab(pub_names[11]) +
  theme(legend.position="none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank()) +
  annotate("text",  x=Inf, y = Inf, label = "H = 0.86; 0.85", vjust= 1.5, hjust=1.05) +
  annotate("text",  x=Inf, y = Inf, label = expression(paste(italic(h)^2, " = 0.55, 0.50")), vjust= 2.5, hjust=1.05)

# yield_per_plant
# edit this dataframe to ensure a space appears for 2018 TLI
emmeans[5170,]$yield_per_plant <- 500

yield_per_plant <- ggplot(emmeans, aes(x = loc, y = yield_per_plant)) + 
  geom_boxplot(aes(fill=interaction(loc, year))) +
  scale_fill_manual(breaks = c("STP.2017", "STP.2018", "TLI.2017", "TLI.2018"), values=c("#0072B2", "#56B4E9", "#C19417", "#fce6a4"), labels = c("STP 2017", "STP 2018", "TLI 2017", "TLI 2018")) +
  coord_cartesian(ylim = c(0, 130), default = TRUE)   + 
  scale_y_continuous(breaks=c(0, 25, 50, 75, 100, 125)) + 
  theme_bw() +
  #ylab(pub_names[12]) + 
  ylab(expression('Yield Plant'^-1*'(g)')) + 
  theme(legend.position="none",
        legend.text=element_text(size=25),
        legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank()) +
  annotate("text",  x=Inf, y = Inf, label = "H = 0.81; NA", vjust= 1.5, hjust=1.05) +
  annotate("text",  x=Inf, y = Inf, label = expression(paste(italic(h)^2, " = 0.45; NA")), vjust= 2.5, hjust=1.05)+
  annotate("text", x = c(0.8, 1.2, 1.8, 2.2), y = c(36, 64, 35, 5), label = c("b", "a", "b", "NA"), size = 4)

# yield_per_spike
yield_per_spike <- ggplot(emmeans, aes(x = loc, y = yield_per_spike)) + 
  geom_boxplot(aes(fill=interaction(loc, year))) +
  scale_fill_manual(values=c("#0072B2", "#C19417", "#56B4E9", "#fce6a4")) +
  scale_y_continuous(breaks=c(0, 0.25, 0.50, 0.75, 1.00, 1.25)) + 
  theme_bw() +
  #ylab(pub_names[13]) +
  ylab(expression('Yield Spike'^-1*'(g)')) +
  theme(legend.position="none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank()) +
  annotate("text",  x=Inf, y = Inf, label = "H = 0.82; 0.55", vjust= 1.5, hjust=1.05) +
  annotate("text",  x=Inf, y = Inf, label = expression(paste(italic(h)^2, " = 0.21; 0.02")), vjust= 2.5, hjust=1.05)


# make a fake figure to extract the legend
fake_fig <- ggplot(emmeans, aes(x = loc, y = yield_per_plant)) + 
  geom_boxplot(aes(fill=interaction(loc, year))) +
  theme_bw() +
  scale_fill_manual(breaks = c("STP.2017", "STP.2018", "TLI.2017", "TLI.2018"), values=c("#0072B2", "#56B4E9", "#C19417", "#fce6a4"), labels = c("STP 2017", "STP 2018", "TLI 2017", "TLI 2018")) +
  theme(legend.text=element_text(size=25),
        legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

leg <- get_legend(fake_fig)
plot(leg)

# cowplot them all together
# save at 1500, 1000 tiff

# depending on the order you want them (alphabetical)
#plot_grid(anthesis_score, emergence_percent, flag_leaf_area, floret_site_utilization, 
          #florets_per_spikelet, height, reproductive_tiller_ct, spikelet_density, 
          #stem_diameter, thousand_grain_weight, yield_per_plant, yield_per_spike, leg)


plot_grid(yield_per_plant, yield_per_spike, reproductive_tiller_ct, 
          thousand_grain_weight, floret_site_utilization, spikelets_per_spike, florets_per_spikelet,
          height, flag_leaf_area, stem_diameter, spikelet_density, emergence_percent, anthesis_score, leg)

