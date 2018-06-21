require(BIEN)
require(tidyr)
require(dplyr)
require(stringr)

#create original WSLA dataset-------------------
WSLA <- BIEN_trait_trait("leaf area per leaf dry mass")
WSLA_raw <- WSLA
GlobalLeafPhenologyDatabase <- read.csv("~/Documents/BEIN data R/GlobalLeafPhenologyDatabase.csv")
Zanne <- GlobalLeafPhenologyDatabase
colnames(Zanne)[colnames(Zanne)=="Binomial"] <- "binomial"
Zanne$binomial <- str_replace_all(Zanne$binomial,"\\s+","_")
colnames(WSLA)[colnames(WSLA)=="scrubbed_species_binomial"] <- "binomial"
WSLA$binomial <- str_replace_all(WSLA$binomial,"\\s+","_")

#########################################################create final_WLSA dataset
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

WSLA_fixed <- subset(WSLA, as.numeric(trait_value)>1 & as.numeric(trait_value)<100)
##### deletes all unreasonable values, theres still unreasonable values in in tho


WSLA_species_count <- WSLA_fixed %>% 
  group_by(binomial) %>% 
  summarize(count=n())
##### gives dataframe with count of species occurances
WSLA_fixed <- left_join(WSLA_fixed, WSLA_species_count, by="binomial")
####joins species count thus species with only 1 count can be identified
####This is such a janky fix 
WSLA_species_family_count <- WSLA_fixed%>% 
  group_by(binomial, scrubbed_family) %>% 
  summarize(count=n())
##### gives dataframe with binomials by family with count
WSLA_fixed_lrgcount <- subset(WSLA_fixed, as.numeric(count)>1)
##### makes a datatable with only species with a count over 1


final_WSLA_DF <- WSLA_fixed_lrgcount[ -c(2,4:13) ]
colnames(final_WSLA_DF)[colnames(final_WSLA_DF)=="trait_value"] <- "SLA"








write.csv(WSLA, file="./data/WSLAPhen.csv")
write.csv(fam.phen.SLA, file="./data/fam_phen_SLA.csv")
write.csv(WSLA_fixed_lrgcount, file="./data/largecount.csv")
saveRDS(final_WSLA_DF, file = "./data/finalWSLADF.rds")
##### creates new dataframe with data, allows me to keep data 



