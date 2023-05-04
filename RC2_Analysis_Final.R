### RC2 Data Analysis

library(plyr)
library(tidyverse)
library(reshape2)
library(ggthemes)
library(ggpubr)
library(gridExtra)
library(picante)
library(vegan)
library(SYNCSA)

source("D:/PNNL Files/Code Information (ESOM and 16S)/FTICR Scripts/R Functions/bMNTD_Fast.R")

# Functions
cond.aro = function(mol){
  mol$aromatics = "Not Aromatic"
  mol$aromatics[which(mol$AI_Mod > 0.5 & mol$AI_Mod <= 0.66)] = "Aromatics"
  mol$aromatics[which(mol$AI_Mod > 0.66)] = "Condensed Aromatics"
  mol$aromatics[which(mol$C > 15 & mol$AI_Mod > 0.66)] = "Polycyclic Aromatics"
  
  return(mol)
}

dist.to.cram = function(mol){
  if(length(which(is.na(mol$NOSC))) > 0){
    stop("You haven't filtered out unassigned peaks.")
  }
  
  # Load DB and clean up
  db = read.csv("D:/PNNL Files/Databases/Hertkorn_CRAM_DB.csv")
  row.names(db) = paste0("CRAM-DB: ", seq(1:nrow(db)))
  db = db[,-1]
  
  # Match mol to DB
  temp.mol = mol[,colnames(db)]
  
  # Merge mol into DB
  db = rbind(db, temp.mol)
  
  # Euclidean distance
  db = vegdist(db, method = "euclidean")
  
  # Convert to matrix
  db = as.matrix(db)
  diag(db) = NA
  
  # Select only DB compounds as comparisons
  db = db[,grep("CRAM", colnames(db))]
  
  # Row averages
  db = rowMeans(db, na.rm = T)
  
  # Remove CRAM compounds
  db = db[-grep("CRAM", names(db))]
  
  # Ensuring identical order with input
  if(!identical(names(db), row.names(temp.mol))){
    stop("Big ol' error text here.")
  }
  
  # Returning result
  return(db)
}

raoQ = function(abund, trait, Hill = TRUE, scale = FALSE, method = "default") {
  
  abund <- as.matrix(abund)
  
  anames <- colnames(abund)[order(colnames(abund))]
  
  abund <- abund[, anames, drop = FALSE]
  
  trait <- as.matrix(trait)
  
  trait <- trait[anames, anames]
  
  if(ncol(abund) != ncol(trait))
    
    stop("Not all species in the abundance matrix appear in the trait matrix!") 
  
  abund <- abund / rowSums(abund)
  
  if(method == "default") 
    
    Q <- apply(abund, 1, function(x) crossprod(x, trait %*% x))
  
  if(method == "divc")
    
    Q <- apply(abund, 1, function(x) x %*% trait^2 %*% (x/2/sum(x)^2))
  
  if(Hill == TRUE) Q <- 1/(1 - Q)
  
  if(scale == TRUE) Q <- Q / max(Q)
  
  names(Q) <- rownames(abund)
  
  return(Q)
  
} # Copied from https://gist.github.com/jslefche/09756ff84afc7b6a82ea0582e663d098; by Jon Lefcheck

rao.pipeline = function(data, mol, package = "base", dist_method = "euclidean", scale = F, compound_type){
  if(missing(compound_type)){
    stop("This function requires that you specify the type of compounds you are working with")
  }
  
  if(package == "SYNCSA"){
    # Keeping user informed to Rao version
    print("You are running the SYNCSA version of Rao's Diversity")
    
    # NOSC Rao
    nosc.rao = as.matrix(vegdist(mol$NOSC, method = dist_method))
    dimnames(nosc.rao) = list(row.names(mol), row.names(mol))
    nosc.rao = rao.diversity(t(data), phylodist = as.matrix(nosc.rao), standardize = F)
    
    # H/C Rao
    hc.rao = as.matrix(vegdist(mol$HtoC_ratio, method = dist_method))
    dimnames(hc.rao) = list(row.names(mol), row.names(mol))
    hc.rao = rao.diversity(t(data), phylodist = as.matrix(hc.rao), standardize = F)
    
    # Carbon Rao
    c.rao = as.matrix(vegdist(mol$C, method = dist_method))
    dimnames(c.rao) = list(row.names(mol), row.names(mol))
    c.rao = rao.diversity(t(data), phylodist = as.matrix(c.rao), standardize = F)
    
    # Multivariate Rao
    multi.rao = scale(mol[,c("C", "HtoC_ratio", "NOSC")])
    multi.rao = as.matrix(vegdist(multi.rao, method = dist_method))
    dimnames(multi.rao) = list(row.names(mol), row.names(mol))
    multi.rao = rao.diversity(t(data), phylodist = as.matrix(multi.rao), standardize = F)
    
    # CRAM Rao
    cram.rao = as.matrix(vegdist(mol$CRAM_Dist, method = dist_method))
    dimnames(cram.rao) = list(row.names(mol), row.names(mol))
    cram.rao = rao.diversity(t(data), phylodist = as.matrix(cram.rao), standardize = F)
    
    # Craft output
    output = rbind(data.frame(Sample = names(nosc.rao$PhyRao), Rao = nosc.rao$PhyRao, Trait = "NOSC", Compound_Type = compound_type), 
                   data.frame(Sample = names(hc.rao$PhyRao), Rao = hc.rao$PhyRao, Trait = "H/C", Compound_Type = compound_type), 
                   data.frame(Sample = names(c.rao$PhyRao), Rao = c.rao$PhyRao, Trait = "Carbon #", Compound_Type = compound_type), 
                   data.frame(Sample = names(multi.rao$PhyRao), Rao = multi.rao$PhyRao, Trait = "Multivariate", Compound_Type = compound_type), 
                   data.frame(Sample = names(cram.rao$PhyRao), Rao = cram.rao$PhyRao, Trait = "Distance to CRAM", Compound_Type = compound_type))
    
  } else {
    # Keeping user informed to Rao version
    print("You are running the Jon Lefcheck's raoQ.")
    
    # NOSC Rao
    nosc.rao = as.matrix(vegdist(mol$NOSC, method = dist_method))
    dimnames(nosc.rao) = list(row.names(mol), row.names(mol))
    nosc.rao = raoQ(t(data), as.matrix(nosc.rao), Hill = F, scale = scale, method = "default")
    
    # H/C Rao
    hc.rao = as.matrix(vegdist(mol$HtoC_ratio, method = dist_method))
    dimnames(hc.rao) = list(row.names(mol), row.names(mol))
    hc.rao = raoQ(t(data), as.matrix(hc.rao), Hill = F, scale = scale, method = "default")
    
    # Carbon Rao
    c.rao = as.matrix(vegdist(mol$C, method = dist_method))
    dimnames(c.rao) = list(row.names(mol), row.names(mol))
    c.rao = raoQ(t(data), as.matrix(c.rao), Hill = F, scale = scale, method = "default")
    
    # Multivariate Rao
    multi.rao = scale(mol[,c("C", "HtoC_ratio", "NOSC")])
    multi.rao = as.matrix(vegdist(multi.rao, method = dist_method))
    dimnames(multi.rao) = list(row.names(mol), row.names(mol))
    multi.rao = raoQ(t(data), as.matrix(multi.rao), Hill = F, scale = scale, method = "default")
    
    # CRAM Rao
    cram.rao = as.matrix(vegdist(mol$CRAM_Dist, method = dist_method))
    dimnames(cram.rao) = list(row.names(mol), row.names(mol))
    cram.rao = raoQ(t(data), as.matrix(cram.rao), Hill = F, scale = scale, method = "default")
    
    # Craft output
    output = rbind(data.frame(Sample = names(nosc.rao), Rao = nosc.rao, Trait = "NOSC", Compound_Type = compound_type), 
                   data.frame(Sample = names(hc.rao), Rao = hc.rao, Trait = "H/C", Compound_Type = compound_type), 
                   data.frame(Sample = names(c.rao), Rao = c.rao, Trait = "Carbon #", Compound_Type = compound_type), 
                   data.frame(Sample = names(multi.rao), Rao = multi.rao, Trait = "Multivariate", Compound_Type = compound_type), 
                   data.frame(Sample = names(cram.rao), Rao = cram.rao, Trait = "Distance to CRAM", Compound_Type = compound_type))
  } # Loop controls which version of the Rao calculation to use - the non-SYNCSA one is significantly faster
  
  return(output)
}
# Quick note about Rao Q usage: the SYNCSA version uses the square root of the Gower distance in its
# calculation of FD. The faster raoQ from Jon Lefcheck's GitHub is agnostic to distance input; I
# prefer usage of Euclidean distances due to "easier" (read: more practical) math.

pull_legend <- function(plot){ 
  temp <- ggplot_gtable(ggplot_build(plot)) 
  leg <- which(sapply(temp$grobs, function(x) x$name) == "guide-box") 
  legend <- temp$grobs[[leg]] 
  legend
} # Copied from stackoverflow

stan.theme = theme(axis.text = element_text(size = 12, color = "black"),
                   axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                   axis.title = element_text(size = 14, color = "black"),
                   axis.ticks = element_line(color = "black"),
                   panel.grid = element_blank(),
                   panel.background = element_blank())

# ############################### #
#### Loading and cleaning data ####
# ############################### #

# Load in data
setwd("D:/PNNL Files/Other 48hr Sampling Campaigns/RC2 Data/Hawkes List/")
data = read.csv("Processed_RC2_Hawkes_Data.csv", row.names = 1, check.names = F)
mol = read.csv("Processed_RC2_Hawkes_Mol.csv", row.names = 1, check.names = F)

# Additional information
meta = read.csv("RC2 Temporal Study Responses (Responses) - Form Responses 1.csv")
catch = read.csv("Rc2_allsites_0717.csv")
npoc = read.csv("NPOC-TN_RC2_Diversity_Paper_11162021.csv")

# Load in transformations
trans = read.csv("RC2_Hawkes_Trans_Profiles.csv", row.names = 1)
trans = trans[,-1]
colnames(trans) = gsub("Sample_", "", colnames(trans))

# Load in trees
mcd = read.tree("RC2_Hawkes_MCD_UPGMA.tre")
td = read.tree("RC2_Hawkes_All-Trans_UPGMA.tre")

# Matching datasets for phylogenetic analyses
mcd.phylo = match.phylo.data(mcd, data)
td.phylo = match.phylo.data(td, data)
rm(td, mcd)

# Setting to pres/abs
data[data > 0] = 1

# Setting up factors
factors = data.frame(Samples = colnames(data), Site = str_extract(colnames(data), "00[0-9][0-9]"), Replicate = str_extract(colnames(data), "ICR.[0-9]"))

# Padding numbers in metadata
meta$Site_Vial_ID_.4_digit_numeric_code. = str_pad(meta$Site_Vial_ID_.4_digit_numeric_code., 4, "left", "0")

# Clean up catch-data
catch = catch[which(catch$Study == "Temporal"), c('Study',"Name","TotDASqKM","D50_m")]
catch$Name[grep(catch$Name,pattern = "American")] = "American River"
catch$Name[grep(catch$Name,pattern = "Union Gap")] = "Union Gap"
catch$Name[grep(catch$Name,pattern = "Little Naches")] = "Little Naches"
catch$Name[grep(catch$Name,pattern = "Mabton")] = "Mabton"
catch$Name[grep(catch$Name,pattern = "Kiona")] = "Kiona"
catch$Name[grep(catch$Name,pattern = "Craig Road 1")] = "Naches- Craig Road 1"
catch$Name[grep(catch$Name,pattern = "Craig Road 2")] = "Naches- Craig Road 2"

# Merge catchment into metadata
meta = merge(x = meta, y = catch, by.x = "Location", by.y = "Name")

# Merge meta into factors
factors = merge(x = factors, y = meta, by.x = "Site", by.y = "Site_Vial_ID_.4_digit_numeric_code.")
rm(meta, catch)

# Adjusting NPOC/TN sample names and merging into factors
npoc$Samples = gsub("_OCN-", "_ICR\\.", npoc$Sample_ID)
npoc$Samples = gsub("$", "_p08", npoc$Samples)
factors = merge(x = factors, y = npoc, by = "Samples")
rm(npoc)

# Reorder data
data = data[,order(colnames(data))]
factors = factors[order(factors$Samples),]
trans = trans[,order(colnames(trans))]

# Confirming identity
if(!identical(colnames(data), factors$Samples)){
  stop("Things don't seem to be matching...")
}

if(!identical(colnames(data), colnames(trans))){
  stop("Things don't seem to be matching...")
}

# Renaming sites to match downstream code
factors$Site = gsub("^", "RC2_", factors$Site)

# Subset object and assigned aromatic status
assigned.mol = mol[!is.na(mol$bs1_class),]
assigned.data = data[!is.na(mol$bs1_class),]
assigned.mol = cond.aro(assigned.mol)
assigned.mol$CRAM_Dist = dist.to.cram(assigned.mol)

# Adding in catchment area object
factors$Catchment_Area = round(factors$TotDASqKM, digits = 0)
factors$Catchment_Area = gsub("$", " sq. km", factors$Catchment_Area)
factors$Catchment_Area = factor(factors$Catchment_Area, levels = c("206 sq. km", "384 sq. km", "2466 sq. km", 
                                                             "8977 sq. km", "13461 sq. km", "14145 sq. km"))

# Add in site IDs
factors$site_ID = case_when(grepl("American River", factors$Location) ~ "T06",
                            grepl("Kiona", factors$Location) ~ "T07",
                            grepl("Little Naches", factors$Location) ~ "T05P",
                            grepl("Mabton", factors$Location) ~ "T02",
                            grepl("Naches- Craig Road 1", factors$Location) ~ "T41",
                            grepl("Naches- Craig Road 2", factors$Location) ~ "T42",
                            grepl("Union Gap", factors$Location) ~ "T03")

# List geospatial data
geo.files = list.files(path = "../Spatial Analyses/Geospatial Data/", pattern = "updated", full.names = T)

# Load in each geospatial data file and merge
geospat = read.csv(geo.files[1]) %>% select(-COMID) # Load initial file to append to

for(f in geo.files[-1]){
  temp = read.csv(f) %>% select(-COMID)
  geospat = geospat %>% left_join(temp, by = "site_ID")
}

# Selecting parameters of interest
geospat = geospat %>% select(site_ID, CAT_AET2015_ANN, TOT_AET2015_ANN, CAT_PET2015_ANN, TOT_PET2015_ANN, CAT_PPT2015_ANN, TOT_PPT2015_ANN, 
                             CAT_TAV2015_ANN, TOT_TAV2015_ANN, CAT_CONTACT, TOT_CONTACT, CAT_RECHG, TOT_RECHG, CAT_urban16, TOT_urban16, 
                             CAT_forest16, TOT_forest16, CAT_wetland16, TOT_wetland16, CAT_agrc16, TOT_agrc16, CAT_shrub16, TOT_shrub16,
                             CAT_BASIN_AREA, TOT_BASIN_AREA, CAT_STREAM_SLOPE, TOT_STREAM_SLOPE, StreamOrde)

# Renaming columns
colnames(geospat)=c("site_ID", "Actual Evapo. '15 (Cat.)", "Actual Evapo. '15 (Tot.)",
                    "Potent. Evapo '15 (Cat.)", "Potent. Evapo '15 (Tot.)", 
                    "Precip. '15 (Cat.)", "Precip. '15 (Tot.)", "Avg. Temp '15 (Cat.)",
                    "Avg. Temp '15 (Tot.)", "Contact (Cat.)", "Contact (Tot.)", 
                    "Recharge (Cat.)", "Recharge (Tot.)", "Urban % '16 (Cat.)", 
                    "Urban % '16 (Tot.)", "Forest % '16 (Cat.)", "Forest % '16 (Tot.)", 
                    "Wetland % '16 (Cat.)", "Wetland % '16 (Tot.)", "Agric. % '16 (Cat.)", 
                    "Agric. % '16 (Tot.)", "Shrub % '16 (Cat.)", "Shrub % '16 (Tot.)", 
                    "Basin Area (Cat.)", "Basin Area (Tot.)", "Stream Slope (Cat.)",
                    "Stream Slope (Tot.)", "Stream Order")
lanuse.names = c("Urban % '16 (Tot.)", "Forest % '16 (Tot.)", "Wetland % '16 (Tot.)",
                 "Agric. % '16 (Tot.)", "Shrub % '16 (Tot.)")

# Merging in geospatial data
factors = factors %>% left_join(geospat, by = "site_ID")


# ############################# #
#### Geospatial Pre-analysis ####
# ############################# #

supp.fig.1 = factors[!duplicated(factors$Catchment_Area),] %>% 
  select(colnames(geospat), TotDASqKM, -site_ID) %>%
  gather(variable, value, -TotDASqKM) %>%
  ggplot(aes(x = TotDASqKM, y = value))+
  geom_point() + geom_smooth(method = "lm")+
  xlab("Catchment Area (km^2)") + ylab("Geospatial Value")+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  facet_wrap(variable~., scales = "free_y")+
  theme_bw() + stan.theme + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))


# ########################### #
#### Analyzing NPOC and TN ####
# ########################### #

# Plotting NPOC and TN by catchment area (Figure 1)
plot1 = factors %>% filter(!(Outlier %in% "TRUE")) %>% select(c(Site, TotDASqKM, NPOC_mg_C_per_L, TN_mg_N_per_L)) %>%
  melt(id.vars = c("Site", "TotDASqKM")) %>%
  group_by(TotDASqKM, variable) %>%
  dplyr::summarise(Mean_Value = mean(value), Median_Value = median(value), StanDev = sd(value)) %>%
  mutate(variable = gsub("NPOC_mg_C_per_L", "NPOC (mg C/L)", variable)) %>%
  mutate(variable = gsub("TN_mg_N_per_L", "TN (mg N/L)", variable)) %>%
  ggplot(aes(x = TotDASqKM, y = Mean_Value))+
  geom_point() + geom_errorbar(aes(ymin=Mean_Value-StanDev, ymax=Mean_Value+StanDev), width = 0)+
  geom_smooth(method = "lm") + 
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  facet_grid(variable~., scales = "free_y") +  
  xlab("Catchment Area (km^2)") + ylab("Average Concentration")+
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"),
                     axis.title = element_text(size = 14, color = "black"),
                     axis.ticks = element_line(color = "black"),
                     panel.border = element_rect(color = "black"),
                     panel.grid = element_blank(),
                     panel.background = element_blank())

# Plotting NPOC and TN by land-cover variables
supp.fig.2 = factors %>% filter(!(Outlier %in% "TRUE")) %>% 
  select(lanuse.names, Catchment_Area, NPOC_mg_C_per_L, TN_mg_N_per_L) %>%
  group_by(Catchment_Area) %>% summarise_all(.funs = "mean") %>%
  gather(Land_Cover, `Land Cover Value`, -Catchment_Area, -NPOC_mg_C_per_L, -TN_mg_N_per_L) %>%
  gather(Chemistry, `Concentration`, -Catchment_Area, -Land_Cover, -`Land Cover Value`) %>%
  mutate(Chemistry = gsub("NPOC_mg_C_per_L", "NPOC (mg C/L)", Chemistry)) %>%
  mutate(Chemistry = gsub("TN_mg_N_per_L", "TN (mg N/L)", Chemistry)) %>%
  ggplot(aes(x = `Land Cover Value`, y = `Concentration`))+
  geom_point() + geom_smooth(method = "lm")+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  facet_grid(Chemistry~Land_Cover, scales = "free")+
  theme_bw() + stan.theme + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))


# ########################## #
#### Molecular properties ####
# ########################## #

# Identifying molecular characteristics
char = data.frame(Samples = colnames(data), `Average Mass` = rep(NA, ncol(data)), 
                     `Average NOSC` = NA, `Average AI-Mod` = NA, `Average DBE` = NA,
                     `Average H/C` = NA, `Average O/C` = NA, stringsAsFactors = F, check.names = F)

for(i in 1:ncol(data)){
  temp = data[which(data[,i] > 0), i, drop = F] # Need to keep names, looking at columns
  temp = mol[which(row.names(mol) %in% row.names(temp)),]
  
  char$`Average Mass`[i] = mean(as.numeric(row.names(temp)), na.rm = T) # Each of these lines is averaging a given characteristic across a sample
  char$`Average NOSC`[i] = mean(temp$NOSC, na.rm = T)
  char$`Average AI-Mod`[i] = mean(temp$AI_Mod, na.rm = T)
  char$`Average DBE`[i] = mean(temp$DBE, na.rm = T)
  char$`Average H/C`[i] = mean(temp$HtoC_ratio, na.rm = T)
  char$`Average O/C`[i] = mean(temp$OtoC_ratio, na.rm = T)
  
  rm("temp")
} # I'm not sure how to do this without the for-loop, but I'm simply just finding the mean for peak stats

# Pairwise statistics for characteristics
char.melt = char %>% gather(Property, value, -Samples) %>%
  left_join(factors, by = "Samples")
char.stat = NULL
uniq.prop = unique(char.melt$Property)
uniq.catc = unique(char.melt$Catchment_Area)

for(curr.prop in uniq.prop){
  temp = char.melt[which(char.melt$Property %in% curr.prop),]
  
  for(i in 1:(length(uniq.catc)-1)){
    for(j in (i+1):length(uniq.catc)){
      w = which(temp$Catchment_Area %in% uniq.catc[c(i,j)])
      stat = wilcox.test(value~Catchment_Area, data = temp[w,])
      stat = data.frame(Property = curr.prop,
                        Comparison = paste0(uniq.catc[i], " vs. ", uniq.catc[j]),
                        W = stat$statistic, p = stat$p.value)
      char.stat = rbind(char.stat, stat)
      
    }
  } # pairwise loop
} # property loop

char.stat$Significant = case_when(char.stat$p > 0.05 ~ "No",
                                  char.stat$p <= 0.05 ~ "Yes")

# Plotting properties by catchment
plot2 = char.melt %>%
  ggplot(aes(x = Catchment_Area, y = value))+
  geom_boxplot() + geom_jitter()+
  stat_compare_means(method = "kruskal.test", label = "p.signif")+
  facet_wrap(.~Property, scales = "free_y")+
  ylab("Average Property Value") + xlab("Catchment Area (km^2)")+
  theme_bw() + stan.theme

# Making figure 1
fig.1 = grid.arrange(ggarrange(plot1, labels = "A"),
                     ggarrange(plot2, labels = "B"),
                     layout_matrix = rbind(c(1,2,2)))
  
rm(plot1, plot2, uniq.prop, uniq.catc, char.melt)

# ############################## #
#### Alpha-diversity analyses ####
# ############################## #

# Counting peaks and transformations
div = NULL
div = cbind(SR = colSums(data), div)
div = as.data.frame(cbind(TR = colSums(trans), div))
div$Normalized_TR = div$TR/div$SR
row.names(div) = gsub("Sample_", "", row.names(div))

# MCD PD
faith = pd(t(mcd.phylo$data), mcd.phylo$phy, include.root = F)
div = data.frame(div, MCD_PD = faith$PD)

# TD PD
faith = pd(t(td.phylo$data), td.phylo$phy, include.root = F)
div = data.frame(div, TD_PD = faith$PD)

# Preparing data for plotting
div = melt(as.matrix(div))
div$Site = gsub("_ICR.*", "", div$Var1)

# Adjusting variable names
div$Variables = case_when(div$Var2  == "TR" ~ "Trans. Count",
                          div$Var2 == "SR" ~ "# of Peaks",
                          div$Var2 == "Normalized_TR" ~ "Trans. Count (Norm.)",
                          div$Var2 == "MCD_PD" ~ "Faith's PD (MCD)",
                          div$Var2 == "TD_PD" ~ "Faith's PD (TD)")

# Ordering factors
div$Variables = factor(div$Variables, levels = c("# of Peaks", "Trans. Count", "Trans. Count (Norm.)", "Faith's PD (MCD)", "Faith's PD (TD)"))

# Plotting mean values against catchment area
alpha.plot = div %>% left_join(factors[, c("Site", "TotDASqKM")], by = "Site") %>%
  mutate(TotDASqKM = as.character(TotDASqKM)) %>%
  group_by(Variables, TotDASqKM) %>% 
  dplyr::summarise(Mean_Value = mean(value), Median_Value = median(value), StanDev = sd(value),
                   StanErr = sd(value)/sqrt(length(value))) %>%
  mutate(TotDASqKM = as.numeric(TotDASqKM)) %>%
  ggplot(aes(x = TotDASqKM, y = Mean_Value))+
  geom_point() + 
  geom_errorbar(aes(ymin=Mean_Value-StanErr, ymax=Mean_Value+StanErr), width = 0)+
  geom_smooth(method = lm) + 
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  facet_grid(Variables~., scales = "free_y")+
  xlab("Catchment Area (km^2)") + ylab("Average Alpha Diversity")+
  theme_bw() + stan.theme + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# Plotting diversity values against NPOC/TN
alpha.npoc = div %>% 
  left_join(factors[, c("Samples", "NPOC_mg_C_per_L", "TN_mg_N_per_L", "Outlier")], by = c("Var1" = "Samples")) %>%
  filter(!(Outlier %in% "TRUE")) %>% select(-Outlier) %>%
  dplyr::rename(Diversity = value) %>% melt(id.vars = c("Var1", "Var2", "Diversity", "Site", "Variables")) %>%
  mutate(variable = gsub("NPOC_mg_C_per_L", "NPOC (mg C/L)", variable)) %>%
  mutate(variable = gsub("TN_mg_N_per_L", "TN (mg N/L)", variable)) %>%
  ggplot(aes(x = value, y = Diversity))+
  geom_point() + geom_smooth(method = "lm") + 
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  facet_grid(Variables~variable, scales = "free")+
  xlab("Concentration") + ylab("Alpha Diversity")+
  theme_bw() + stan.theme

# Plotting diversity values against land cover values
alpha.lanuse = div %>% 
  left_join(factors[, c("Samples", "Catchment_Area")], 
            by = c("Var1" = "Samples")) %>%
  select(-Var1) %>% group_by(Catchment_Area, Variables) %>%
  summarise(Mean_Value = mean(value), Median_Value = median(value), StanDev = sd(value),
            StanErr = sd(value)/sqrt(length(value))) %>%
  left_join(factors[!duplicated(factors$Catchment_Area),
                    c("Catchment_Area", lanuse.names)], by = "Catchment_Area") %>%
  gather(Land_Cover, `Land Cover Value`, -Catchment_Area, 
         -Variables, -Mean_Value, -Median_Value, -StanDev, -StanErr) %>%
  ggplot(aes(x = `Land Cover Value`, y = Mean_Value))+
  geom_point() + 
  geom_errorbar(aes(ymin=Mean_Value-StanErr, ymax=Mean_Value+StanErr), width = 0)+
  geom_smooth(method = "lm") + 
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  facet_grid(Variables~Land_Cover, scales = "free")+
  xlab("Land Cover %") + ylab("Average Alpha Diversity")+
  theme_bw() + stan.theme + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# Generating alpha diversity figure (Figure 2)
fig.2 = grid.arrange(ggarrange(alpha.plot, labels = "A"),
                     ggarrange(alpha.lanuse, labels = "B"),
                     layout_matrix = rbind(c(1,2,2,2),
                                           c(1,2,2,2)))
supp.fig.3 = alpha.npoc

rm(alpha.npoc, alpha.lanuse, alpha.plot)


# ############################# #
#### Beta-diversity analyses ####
# ############################# #

## Generate Bray-Curtis PCoA
dist = vegdist(t(data), method = "jaccard")
beta = pcoa(dist)

# Arrange data
beta.sig = beta$values
beta = data.frame(Samples = row.names(beta$vectors), beta$vectors) %>% left_join(factors, by = "Samples")

# Plot temporary PCoA
bray.pcoa = ggplot(beta, aes(x = Axis.1, y = Axis.2))+
  geom_point(aes(fill = Catchment_Area), color = "black", pch = 21, size = 5)+
  stat_ellipse(aes(color = Catchment_Area))+
  scale_fill_tableau(name = "Catchment Area")+
  scale_color_tableau(name = "Catchment Area")+
  xlab(paste0("Axis 1 (", round(beta.sig$Relative_eig[1]*100, digits = 2), "%)"))+ 
  ylab(paste0("Axis 2 (", round(beta.sig$Relative_eig[2]*100, digits = 2), "%)"))+
  theme_bw() + theme(text = element_text(size = 14),
                     axis.text = element_text(color = "black", size = 14),
                     axis.title = element_text(color = "black", size = 16),
                     axis.ticks = element_line(color = "black"),
                     panel.border = element_rect(size = 1, color = "black"),
                     panel.background = element_blank(),
                     panel.grid = element_blank())

# Pull legend
leg.obj = pull_legend(bray.pcoa)

# Plot real PCoA
bray.pcoa = ggplot(beta, aes(x = Axis.1, y = Axis.2))+
  geom_point(aes(fill = Catchment_Area), color = "black", pch = 21, size = 5)+
  stat_ellipse(aes(color = Catchment_Area))+
  scale_fill_tableau(name = "Catchment Area")+
  scale_color_tableau(name = "Catchment Area")+
  xlab(paste0("Axis 1 (", round(beta.sig$Relative_eig[1]*100, digits = 2), "%)"))+ 
  ylab(paste0("Axis 2 (", round(beta.sig$Relative_eig[2]*100, digits = 2), "%)"))+
  theme_bw() + theme(legend.position = "none",
                     text = element_text(size = 14),
                     axis.text = element_text(color = "black", size = 14),
                     axis.title = element_text(color = "black", size = 16),
                     axis.ticks = element_line(color = "black"),
                     panel.border = element_rect(size = 1, color = "black"),
                     panel.background = element_blank(),
                     panel.grid = element_blank())

## Generate MCD bMNTD PCoA
# Generate PCoA
dist = comdistnt_fast(comm = t(mcd.phylo$data), dis = cophenetic(mcd.phylo$phy), abundance.weighted = F)
beta = pcoa(dist)

# Arrange data
beta.sig = beta$values
beta = data.frame(Samples = row.names(beta$vectors), beta$vectors) %>% left_join(factors, by = "Samples")

# Plot PCoA
mcd.pcoa = ggplot(beta, aes(x = Axis.1, y = Axis.2))+
  geom_point(aes(fill = Catchment_Area), color = "black", pch = 21, size = 5)+
  stat_ellipse(aes(color = Catchment_Area))+
  scale_fill_tableau(name = "Catchment Area")+
  scale_color_tableau(name = "Catchment Area")+
  xlab(paste0("Axis 1 (", round(beta.sig$Rel_corr_eig[1]*100, digits = 2), "%)"))+ 
  ylab(paste0("Axis 2 (", round(beta.sig$Rel_corr_eig[2]*100, digits = 2), "%)"))+
  theme_bw() + theme(legend.position = "none",
                     text = element_text(size = 14),
                     axis.text = element_text(color = "black", size = 14),
                     axis.title = element_text(color = "black", size = 16),
                     axis.ticks = element_line(color = "black"),
                     panel.border = element_rect(size = 1, color = "black"),
                     panel.background = element_blank(),
                     panel.grid = element_blank())

## Generate MCD bMNTD PCoA
# Generate PCoA
dist = comdistnt_fast(comm = t(td.phylo$data), dis = cophenetic(td.phylo$phy), abundance.weighted = F)
beta = pcoa(dist)

# Arrange data
beta.sig = beta$values
beta = data.frame(Samples = row.names(beta$vectors), beta$vectors) %>% left_join(factors, by = "Samples")

# Plot PCoA
td.pcoa = ggplot(beta, aes(x = Axis.1, y = Axis.2))+
  geom_point(aes(fill = Catchment_Area), color = "black", pch = 21, size = 5)+
  stat_ellipse(aes(color = Catchment_Area))+
  scale_fill_tableau(name = "Catchment Area")+
  scale_color_tableau(name = "Catchment Area")+
  xlab(paste0("Axis 1 (", round(beta.sig$Rel_corr_eig[1]*100, digits = 2), "%)"))+ 
  ylab(paste0("Axis 2 (", round(beta.sig$Rel_corr_eig[2]*100, digits = 2), "%)"))+
  theme_bw() + theme(legend.position = "none",
                     text = element_text(size = 14),
                     axis.text = element_text(color = "black", size = 14),
                     axis.title = element_text(color = "black", size = 16),
                     axis.ticks = element_line(color = "black"),
                     panel.border = element_rect(size = 1, color = "black"),
                     panel.background = element_blank(),
                     panel.grid = element_blank())


# ######################## #
#### Arranging Figure 4 ####
# ######################## #
supp.fig.4 = grid.arrange(ggarrange(bray.pcoa, labels = "A"),
                     ggarrange(mcd.pcoa, labels = "B"),
                     ggarrange(td.pcoa, labels = "C"),
                     leg.obj,
                     layout_matrix = rbind(c(1, 2), 
                                           c(3, 4)))


# ################################# #
#### Testing out Rao's Diversity ####
# ################################# #

### Whole community Rao
# Total Rao
total.rao = rao.pipeline(assigned.data, assigned.mol, scale = T, compound_type = "Total")

# Converting columns
total.rao$Site = str_extract(row.names(total.rao), "RC2_00[0-9][0-9]")
total.rao$Sample = row.names(total.rao)
total.rao$Sample = gsub("_p081|_p082|_p083|_p084", "_p08", total.rao$Sample)

# Ordering factors
total.rao$Trait = factor(total.rao$Trait, levels = c("Carbon #", "H/C", "NOSC", "Multivariate", "Distance to CRAM"))

# Plotting NPOC/TN values against Rao's diversity
tot.rao.conc = total.rao %>% 
  left_join(factors[, c("Samples", "NPOC_mg_C_per_L", "TN_mg_N_per_L", "Outlier")], by = c("Sample" = "Samples")) %>%
  filter(!(Outlier %in% "TRUE")) %>% select(-Outlier) %>%
  melt(id.vars = c("Sample", "Rao", "Trait", "Compound_Type", "Site")) %>%
  mutate(variable = gsub("NPOC_mg_C_per_L", "NPOC (mg C/L)", variable)) %>%
  mutate(variable = gsub("TN_mg_N_per_L", "TN (mg N/L)", variable)) %>%
  ggplot(aes(x = value, y = Rao))+
  geom_point(aes(color = Trait)) + geom_smooth(method = "lm") + 
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  facet_grid(Trait~variable, scales = "free") + scale_color_solarized()+
  xlab("Concentration") + ylab("Rao's FD")+
  theme_bw() + theme(legend.position = "none",
                     text = element_text(size = 14),
                     axis.text.x = element_text(colour = "black"),
                     axis.text.y = element_text(colour = "black"),
                     panel.border = element_rect(size = 1, colour = "black"),
                     panel.background = element_blank(),
                     panel.grid = element_blank())


# Plotting average rao by catchment area
tot.rao.catch = total.rao %>% left_join(factors[, c("Site", "TotDASqKM")], by = "Site") %>%
  # distinct(Sample, .keep_all = T) %>%
  group_by(Trait, TotDASqKM) %>%
  dplyr::summarise(Median_Rao = median(Rao), Mean_Rao = mean(Rao), StanDev = sd(Rao),
                   StanErr = sd(Rao)/sqrt(length(Rao))) %>%
  ggplot(aes(x = TotDASqKM, y = Mean_Rao))+
  geom_point(aes(color = Trait))+ 
  geom_errorbar(aes(ymin=Mean_Rao-StanErr, ymax=Mean_Rao+StanErr), width = 0)+
  geom_smooth(method = "lm")+ 
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  facet_grid(Trait~., scales = "free_y") + scale_color_solarized()+
  xlab("Catchment Area (km^2)") + ylab("Average Rao's FD")+
  theme_bw() + theme(legend.position = "none",
                     text = element_text(size = 14),
                     axis.text.x = element_text(colour = "black"),
                     axis.text.y = element_text(colour = "black"),
                     panel.border = element_rect(size = 1, colour = "black"),
                     panel.background = element_blank(),
                     panel.grid = element_blank())

# Plotting average rao by land cover
tot.rao.lanuse = total.rao %>% left_join(factors[, c("Site", "Catchment_Area")], by = "Site") %>%
  group_by(Trait, Catchment_Area) %>%
  summarise(Median_Rao = median(Rao), Mean_Rao = mean(Rao), StanDev = sd(Rao), 
            StanErr = sd(Rao)/sqrt(length(Rao))) %>%
  left_join(factors[,c(lanuse.names, "Catchment_Area")], by = "Catchment_Area") %>%
  group_by(Trait, Catchment_Area) %>% summarise_all(.funs = "mean") %>%
  gather(Land_Cover, `Land Cover Value`, -Trait, -Catchment_Area, -Median_Rao, 
         -Mean_Rao, -StanDev, -StanErr) %>%
  ggplot(aes(x = `Land Cover Value`, y = Mean_Rao))+
  geom_point(aes(color = Trait))+ 
  geom_errorbar(aes(ymin=Mean_Rao-StanErr, ymax=Mean_Rao+StanErr), width = 0)+
  geom_smooth(method = "lm")+ 
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  facet_grid(Trait~Land_Cover, scales = "free") + scale_color_solarized()+
  xlab("Land Cover (%)") + ylab("Average Rao's FD")+
  theme_bw() + theme(legend.position = "none",
                     text = element_text(size = 14),
                     axis.text.x = element_text(colour = "black"),
                     axis.text.y = element_text(colour = "black"),
                     panel.border = element_rect(size = 1, colour = "black"),
                     panel.background = element_blank(),
                     panel.grid = element_blank())



### Merge plots
# sup.fig2 = ggarrange(plot1, plot2, common.legend = T, labels = c("A", "B"))
# sup.figx = ggarrange(plot.a, plot.b, labels = c("A", "B"))
fig.3 = grid.arrange(ggarrange(tot.rao.catch, labels = "A"),
                     ggarrange(tot.rao.lanuse, labels = "B"),
                     layout_matrix = rbind(c(1,2,2,2),
                                           c(1,2,2,2)))
supp.fig.5 = tot.rao.conc

### Plotting transformation counts vs. Rao's
# Generating average datasets for direct comparisons
mean.rao = total.rao %>% left_join(factors[, c("Site", "Catchment_Area")], by = "Site") %>%
  group_by(Trait, Catchment_Area) %>%
  summarise(Median_Rao = median(Rao), Mean_Rao = mean(Rao), StanDev = sd(Rao), 
            StanErr = sd(Rao)/sqrt(length(Rao)))
mean.div = div %>% left_join(factors[, c("Site", "Catchment_Area")], by = "Site") %>%
  group_by(Variables, Catchment_Area) %>% 
  summarise(Mean_Value = mean(value), Median_Value = median(value), StanDev = sd(value),
            StanErr = sd(value)/sqrt(length(value)))

# Plotting direct comparisons between trans and Rao
mean.rao %>% left_join(mean.div, by = "Catchment_Area") %>% 
  filter(Trait %in% c("H/C", "NOSC")) %>%
  filter(Variables %in% c("Trans. Count", "Trans. Count (Norm.)")) %>%
  ggplot(aes(x = Mean_Rao, y = Mean_Value))+
  geom_smooth(method = "lm")+ 
  geom_point(color = "black", fill = "white", size = 4, pch = 21)+
  geom_errorbar(aes(ymin=Mean_Value-StanErr.y, ymax=Mean_Value+StanErr.y), width = 0)+
  geom_errorbarh(aes(xmin=Mean_Rao-StanErr.x, xmax=Mean_Rao+StanErr.x), height = 0)+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  facet_grid(Variables~Trait, scales = "free")+ 
  xlab("Average Rao's FD") + ylab("Average Alpha Diversity")+
  theme_bw() + theme(legend.position = "none",
                     text = element_text(size = 14),
                     axis.text.x = element_text(colour = "black"),
                     axis.text.y = element_text(colour = "black"),
                     panel.border = element_rect(size = 1, colour = "black"),
                     panel.background = element_blank(),
                     panel.grid = element_blank())
