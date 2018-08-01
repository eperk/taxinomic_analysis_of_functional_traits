#loading all required packages-------------------------------------------------------------------------
require(BIEN)
require(tidyr)
require(dplyr)
require(stringr)
require(phangorn)
require(ape)
require(phytools)
require(phylolm)
require(MPSEM)
library(picante)
library(readr)
library(ggplot2)
require(mosaic)
require(lattice)
require(nlme)
library(geiger)
library(sp)
require(foreach)
require(phylobase)
require(tibble)
require(miscTools)
require(lme4)

#create original WSLA dataset-------------------------------------------------------------------------
WSLA <- BIEN_trait_trait("leaf area per leaf dry mass")
WSLA_raw <- WSLA
GlobalLeafPhenologyDatabase <- read.csv("~/Documents/BEIN data R/data/raw/GlobalLeafPhenologyDatabase.csv")
Zanne <- GlobalLeafPhenologyDatabase
colnames(Zanne)[colnames(Zanne)=="Binomial"] <- "binomial"
Zanne$binomial <- str_replace_all(Zanne$binomial,"\\s+","_")
colnames(WSLA)[colnames(WSLA)=="scrubbed_species_binomial"] <- "binomial"
WSLA$binomial <- str_replace_all(WSLA$binomial,"\\s+","_")

#create final_WLSA dataset-------------------------------------------------------------------------
WSLA$LMA <- (1/as.numeric(WSLA$trait_value))*1000
#####  is.numeric = asks r if a column in a dataframe if somethings numeric
#####  as.numeric = makes r say a column is numeric
#####  $ = adds something / defines a column in a dataframe example: dataframe$column
#####  g per m2
WSLA <- left_join(WSLA, Zanne, by="binomial")
##### Joins stuff in Zanne by column "Binomial"

binomial <- WSLA_raw[,1, drop=FALSE]
#####Binomial data only from WSLA dataframe
WSLA_species_by_family <-c (binomial$scrubbed_species_binomial)
#####gives vector of all species in dataframe binomial
taxonomy_data<- BIEN_taxonomy_species(WSLA_species_by_family)
##### Gives all taxonomic data of species in dataframe binomial
family_binomial <- taxonomy_data[ -c(1:4,6,8:9) ] 
##### leaves me with just species binomal and family binomial
colnames(family_binomial)[colnames(family_binomial)=="scrubbed_species_binomial"] <- "binomial"
family_binomial$binomial <- str_replace_all(family_binomial$binomial,"\\s+","_")

WSLA <- left_join(WSLA, family_binomial, by="binomial")

##### deletes all unreasonable values-------------------------------------------------------------------------
WSLA_fixed <- subset(WSLA, as.numeric(trait_value)>1 & as.numeric(trait_value)<100)



##### gives dataframe with count of species occurances-------------------------------------------------------------------------
WSLA_species_count <- WSLA_fixed %>% 
  group_by(binomial) %>% 
  summarize(count=n())
WSLA_fixed <- left_join(WSLA_fixed, WSLA_species_count, by="binomial")
####joins species count thus species with only 1 count can be identified
####This is such a janky fix 
WSLA_species_family_count <- WSLA_fixed%>% 
  group_by(binomial, scrubbed_family) %>% 
  summarize(count=n())
##### gives dataframe with binomials by family with count
WSLA_fixed_lrgcount <- subset(WSLA_fixed, as.numeric(count)>1)
##### makes a datatable with only species with a count over 1

#creation of clean datasets for use in fossil tree integraiton-------------------------------------------------------------------------
florissant_fossil <- read_csv("~/Documents/BEIN data R/data/raw/florissant_fossil.csv")
renova_fossil <- read_csv("~/Documents/BEIN data R/data/raw/renova_fossil.csv")
bridgecreek_fossil <- read_csv("~/Documents/BEIN data R/data/raw/renova_fossil.csv")

florissant_fossil <- florissant_fossil[-c(1,7)]
renova_fossil <- renova_fossil[-c(1,7,10:12)]
bridgecreek_fossil <- bridgecreek_fossil[-c(1,7,10:12)]

florissant_fossil$species [is.na(florissant_fossil$species)] <- "sp."
florissant_fossil <- na.omit(florissant_fossil)

renova_fossil$species [is.na(renova_fossil$species)] <- "sp."
renova_fossil <- na.omit(renova_fossil)

bridgecreek_fossil$species [is.na(bridgecreek_fossil$species)] <- "sp."
bridgecreek_fossil <- na.omit(bridgecreek_fossil)

saveRDS(florissant_fossil, file = "~/Documents/BEIN data R/data/processed/04_florissant_fossil_clean.rds")
saveRDS(renova_fossil, file = "~/Documents/BEIN data R/data/processed/04_renova_fossil_clean.rds")
saveRDS(bridgecreek_fossil, file = "~/Documents/BEIN data R/data/processed/04_bridgecreek_fossil_clean.rds")

#fossil datasheet for phylogeny integration-----------
florissant_fossil_phylo <- readRDS("~/Documents/BEIN data R/data/processed/04_florissant_fossil_clean.rds")
renova_fossil_phylo <- readRDS("~/Documents/BEIN data R/data/processed/04_renova_fossil_clean.rds")
bridgecreek_fossil_phylo <- readRDS("~/Documents/BEIN data R/data/processed/04_bridgecreek_fossil_clean.rds")

renova_fossil_phylo <- renova_fossil_phylo[-c(7)]
bridgecreek_fossil_phylo <- bridgecreek_fossil_phylo[-c(7)]

all_fossil_phyloint <- rbind(renova_fossil_phylo, bridgecreek_fossil_phylo, florissant_fossil_phylo)

all_fossil_phyloint$binomial <- paste(all_fossil_phyloint$Genus, all_fossil_phyloint$species)
all_fossil_phyloint$binomial <- str_replace_all(all_fossil_phyloint$binomial,"\\s+","_")

all_fossil_phyloint$`Approximate Date (Ma)` <- str_replace_all(all_fossil_phyloint$`Approximate Date (Ma)`,"Early Oligocene","33900000")
all_fossil_phyloint$`Approximate Date (Ma)` <- str_replace_all(all_fossil_phyloint$`Approximate Date (Ma)`,"Late Eocene","37800000")
all_fossil_phyloint$`Approximate Date (Ma)` <- str_replace_all(all_fossil_phyloint$`Approximate Date (Ma)`,"34.07","34070000")

colnames(all_fossil_phyloint)[colnames(all_fossil_phyloint)=="Approximate Date (Ma)"] <- "date"

saveRDS(all_fossil_phyloint, file = "~/Documents/BEIN data R/data/processed/04_fossil_phylo_integration.rds")

all_fossil_phyloint$date <- as.numeric(all_fossil_phyloint$date)
#creation of matched dataset where i have species with full taxonomy info------
missing_genus <- subset(all_fossil_phyloint, !(all_fossil_phyloint$Genus %in% fossil_tax$scrubbed_genus))
all_fossil_phyloint_match<-subset(all_fossil_phyloint, 
  !(all_fossil_phyloint$Genus %in% missing_genus$Genus))

#integrate fossil to all data-------------------------------------------------------------------------
florissant_fossil_int <- readRDS("~/Documents/BEIN data R/data/processed/04_florissant_fossil_clean.rds")
renova_fossil_int <- readRDS("~/Documents/BEIN data R/data/processed/04_renova_fossil_clean.rds")
bridgecreek_fossil_int <- readRDS("~/Documents/BEIN data R/data/processed/04_bridgecreek_fossil_clean.rds")

renova_fossil_int <- renova_fossil_int[-c(6,7)]
bridgecreek_fossil_int <- bridgecreek_fossil_int[-c(6,7)]
florissant_fossil_int <- florissant_fossil_int[-c(6)]

all_fossil <- rbind(florissant_fossil_int, renova_fossil_int, bridgecreek_fossil_int)
all_fossil$binomial <- paste(all_fossil$Genus, all_fossil$species)
all_fossil$binomial <- str_replace_all(all_fossil$binomial,"\\s+","_")
all_fossil <- all_fossil[,c(6,1:5)]
all_fossil <- all_fossil[-c(2,3)]

royer_data_fossil_int <- royer_data_sub[-c(3,5)]

colnames(all_fossil)[colnames(all_fossil)=="Petiole Width (cm)"] <- "avg_petiole_width"
colnames(all_fossil)[colnames(all_fossil)=="Leaf Area (cm^2)"] <- "avg_LA"
colnames(all_fossil)[colnames(all_fossil)=="PW^2/A"] <- "log_pet_leafarea"
all_fossil$log_pet_leafarea <- log(all_fossil$log_pet_leafarea)

royer_data_fossil_int <- rbind(all_fossil,royer_data_fossil_int)

royer_data_fossil_int <- aggregate(royer_data_fossil_int[,-1], by=list(royer_data_fossil_int$binomial), mean)
colnames(royer_data_fossil_int)[colnames(royer_data_fossil_int)=="Group.1"] <- "binomial"
rownames(royer_data_fossil_int)<-royer_data_fossil_int$binomial
royer_data_fossil_int <- royer_data_fossil_int[-c(1)]

#lookup table creation for intfossil function----------------------------------------------------------------------------------------------
fossil_comb <- rbind(florissant_fossil_int,renova_fossil_int,bridgecreek_fossil_int)
fossil_gen <- fossil_comb$Genus
fossil_gen <- unique(fossil_gen)
fossil_gen_df <- as.data.frame(fossil_gen)
fossil_tax <- BIEN_taxonomy_genus(fossil_gen)
fossil_tax <- fossil_tax[-c(1,7:9)]
fossil_tax <- unique(fossil_tax)

#final dataframe creation-------------------------------------------------------------------------------------------------
final_WSLA_DF <- WSLA_fixed_lrgcount[ -c(2,4:13) ]
colnames(final_WSLA_DF)[colnames(final_WSLA_DF)=="trait_value"] <- "SLA"

#creating of trees-------------------------------------------------------------------------------------------------
tree_plant <- read.tree('~/Documents/BEIN data R/data/raw/phylodata/Vascular_Plants_rooted.dated.tre')
#from https://datadryad.org//resource/doi:10.5061/dryad.63q27

tree_tips <- tree_plant$tip.label
###### all tips in og tree
tree_tips_df <- as.data.frame(tree_tips) 
###### makes a dataframe
colnames(tree_tips_df)[colnames(tree_tips_df)=="tree.tips"] <- "binomial"

all_species_df <- as.data.frame(unique(final_WSLA_DF$binomial))
colnames(all_species_df)[colnames(all_species_df)=="unique(final_WSLA_DF$binomial)"] <- "binomial"

not_in_WSLA <- subset(tree_tips, !(tree_tips %in% all_species_df$binomial))
####subset of stuff that isn't in all.species.df but is in the tree
#########this subseting function is so important

tree_WSLA_species <- drop.tip(tree_plant, not_in_WSLA)
WSLA_tree_tips <- tree_WSLA_species$tip.label
WSLA_tree_tips_df <- as.data.frame(WSLA_tree_tips)
#####plots Zanne tree
#####drops non WSLA bionomials

tree_massive <-read.tree('~/Documents/BEIN data R/data/raw/phylodata/ALLMB.tre')
#from https://github.com/FePhyFoFum/big_seed_plant_trees/releases



#cleaning data for inital PEM---------------
royer_data <- read_csv("~/Documents/BEIN data R/data/raw/royer_data.csv")

royer_data$binomial <- str_replace_all(royer_data$binomial,"\\s+","_")

royer_data_sub<-
  royer_data %>% 
  dplyr::filter(binomial%in%tree_plant$tip.label) %>% 
  group_by(binomial) %>% 
  summarise(avg_petiole_width=mean(petiole_width, na.rm=TRUE),
            avg_LMA=mean(LMA, na.rm=TRUE), avg_LA=mean(leaf_area)) %>% 
  as.data.frame()



royer_data_sub$log_lma <- log(royer_data_sub$avg_LMA)
royer_data_sub$log_pet_leafarea <- log(royer_data_sub$avg_petiole_width^2/royer_data_sub$avg_LA)

rownames(royer_data_sub)<-royer_data_sub$binomial

indx<-which(tree_plant$tip.label%in%royer_data_sub$binomial==FALSE)

royer_tree <- drop.tip(tree_plant, tree_plant$tip.label[indx])
royer_tips <- royer_tree$tip.label


spmatch <- match(royer_tree$tip.label,
                 royer_data_sub[,1L])
royer_match_data <- royer_data_sub[spmatch,]


royer_phylo_lma<-match.phylo.data(royer_tree,royer_data_sub)
royer_phylo_lma$data$avg_petiole_width<-as.numeric(royer_phylo_lma$data$avg_petiole_width)
royer_phylo_lma$data$avg_LMA<-as.numeric(royer_phylo_lma$data$avg_LMA)
royer_phylo_lma$data$log_lma<-as.numeric(royer_phylo_lma$data$log_lma)
royer_phylo_lma$data$log_pet_leafarea<-as.numeric(royer_phylo_lma$data$log_pet_leafarea)

#royer data cleaning for LME4------------
royer_data_LME4<-
  royer_data %>% 
  group_by(binomial) %>% 
  summarise(avg_petiole_width=mean(petiole_width, na.rm=TRUE),
            avg_LMA=mean(LMA, na.rm=TRUE), avg_LA=mean(leaf_area)) %>% 
  as.data.frame()

royer_data_LME4$log_lma <- log(royer_data_LME4$avg_LMA)
royer_data_LME4$log_pet_leafarea <- log(royer_data_LME4$avg_petiole_width^2/royer_data_LME4$avg_LA)

royer_data_LME4$binomial <- as.character(gsub("\\?", "", royer_data_LME4$binomial))
royer_data_LME4$binomial <- as.character(gsub("\\#", "", royer_data_LME4$binomial))
royer_data_LME4$binomial <- as.character(gsub("\\_1", "1", royer_data_LME4$binomial))
royer_data_LME4$binomial <- as.character(gsub("\\_2", "2", royer_data_LME4$binomial))
royer_data_LME4$binomial <- as.character(gsub("\\_3", "3", royer_data_LME4$binomial))
royer_data_LME4$binomial <- as.character(gsub("\\_4", "4", royer_data_LME4$binomial))
royer_data_LME4$binomial <- as.character(gsub("\\_5", "5", royer_data_LME4$binomial))
royer_data_LME4$binomial <- as.character(gsub("\\_cf.", "", royer_data_LME4$binomial))
royer_data_LME4$binomial <- as.character(gsub("\\_\\(trilobum)", "", royer_data_LME4$binomial))
royer_data_LME4$binomial <- as.character(gsub("\\-herculis", "", royer_data_LME4$binomial))
royer_data_LME4$binomial <- as.character(gsub("\\_\\(nitidissima)", "", royer_data_LME4$binomial))
royer_data_LME4$binomial <- as.character(gsub("\\_var_remyi", "", royer_data_LME4$binomial))


royer_data_LME4<-
  royer_data_LME4 %>% 
  filter(!grepl('unknown', binomial))

royer_data_LME4<-
  royer_data_LME4 %>% 
  filter(!grepl('unkown', binomial))

royer_data_LME4<-
  royer_data_LME4 %>% 
  filter(!grepl('cf._Lasiopetalum', binomial))


royer_data_LME4_split <-
  royer_data_LME4 %>% separate(binomial, 
                c("genus", "species"))

colnames(royer_data_LME4_split)[colnames(royer_data_LME4_split)=="genus"] <- "scrubbed_genus"

royer_tax_data <- BIEN_taxonomy_genus(royer_data_LME4_split$genus)
royer_tax_data_sub <- subset(royer_tax_data, royer_data_LME4_split$genus %in% royer_tax_data$scrubbed_genus)
royer_tax_data_sub <- royer_tax_data_sub [-c(1,7,8:9)]
royer_tax_data_sub <- unique(royer_tax_data_sub)
royer_tax_data_sub <- na.omit(royer_tax_data_sub)
royer_tax_full <- left_join(royer_data_LME4_split, royer_tax_data_sub, by= "scrubbed_genus")

#creation of fossil newdata for preidciton-------
 all_fossil_LMEpred <- all_fossil %>%
  separate(binomial, 
    c("scrubbed_genus", "species"))
all_fossil_royer_pred <- left_join(all_fossil_LMEpred, fossil_tax, by = "scrubbed_genus")

all_fossil_royer_pred <- na.omit(all_fossil_royer_pred)

#clean data saving-------------------------------------------------------------------------------------------------
saveRDS(WSLA_raw, file="~/Documents/BEIN data R/data/processed/00_WSLAraw.rds")
saveRDS(WSLA, file="~/Documents/BEIN data R/data/processed/02_WSLAPhen.rds")
saveRDS(WSLA_fixed_lrgcount, file="~/Documents/BEIN data R/data/processed/02_family_count.rds")
saveRDS(final_WSLA_DF, file = "~/Documents/BEIN data R/data/processed/03_finalWSLADF.rds")
saveRDS(royer_data_fossil_int, file="~/Documents/BEIN data R/data/processed/04_royer_data_fossil_int.rds")
write.tree(tree_WSLA_species, file = "~/Documents/BEIN data R/data/processed/03_WSLA_species.tre")
saveRDS(royer_phylo_lma, "~/Documents/BEIN data R/data/processed/05_royer_clean_df_phylo.rds")
saveRDS(royer_tax_full, "~/Documents/BEIN data R/data/processed/07_lm4_royer")
