require(BIEN)
require(tidyr)
require(dplyr)

BIEN_trait_list()
trait_vector <-c("whole plant height","leaf dry mass per leaf fresh mass")
BIEN_trait_trait(trait_vector)
WSLA <- BIEN_trait_trait("whole plant leaf area per whole plant leaf dry mass")
View(WSLA)
Zanne <- GlobalLeafPhenologyDatabase

write.csv(WSLA, file="./data/WSLAPhen.csv")

#####creates new dataframe with data, allows me to keep data 

WSLA_raw <- WSLA
##### Saves me an original version of the data


WSLA$LMA <- 1/as.numeric(WSLA$trait_value)
#####  is.numeric = asks r if a column in a dataframe if somethings numeric
#####  as.numeric = makes r say a column is numeric
#####  $ = adds something / defines a column in a dataframe example: dataframe$column

colnames(WSLA)[colnames(WSLA)=="scrubbed_species_binomial"] <- "Binomial"
##### Renames a column

WSLA <- left_join(WSLA, Zanne, by="Binomial")
##### Joins stuff in Zanne by column "Binomial"

PhenEv <- subset(WSLA, WSLA$Phenology == "EV", Drop=TRUE)
PhenD <- subset(WSLA, WSLA$Phenology == "D", Drop=TRUE)
##### Subsets out of WSLA only Plants with Phenology of  "EV" or "D"