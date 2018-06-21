#######creation of clean datasets for use in fossil tree integraiton
florissant_fossil <- read_csv("./data/florissant_fossil.csv")
renova_fossil <- read_csv("./data/renova_fossil.csv")
bridgecreek_fossil <- read_csv("./data/renova_fossil.csv")

florissant_fossil <- florissant_fossil[-c(1,7)]
renova_fossil <- renova_fossil[-c(1,7,10:12)]
bridgecreek_fossil <- bridgecreek_fossil[-c(1,7,10:12)]

florissant_fossil$species [is.na(florissant_fossil$species)] <- "sp."
florissant_fossil <- na.omit(florissant_fossil)

renova_fossil$species [is.na(renova_fossil$species)] <- "sp."
renova_fossil <- na.omit(renova_fossil)

bridgecreek_fossil$species [is.na(bridgecreek_fossil$species)] <- "sp."
bridgecreek_fossil <- na.omit(bridgecreek_fossil)

saveRDS(florissant_fossil, file = "./data/florissant_fossil_clean.rds")
saveRDS(renova_fossil, file = "./data/renova_fossil_clean.rds")
saveRDS(bridgecreek_fossil, file = "./data/bridgecreek_fossil_clean.rds")


######################phytools fossil integration
florissant_fossil_int <- readRDS("./data/florissant_fossil_clean.rds")
renova_fossil_int <- readRDS("./data/renova_fossil_clean.rds")
bridgecreek_fossil_int <- readRDS("./data/bridgecreek_fossil_clean.rds")

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

saveRDS(royer_data_fossil_int, file="./data/foyer_data_fossil_int")


r_fossilrix <- as.matrix(royer_data_fossil_int)
r_fossilrix <- r_fossilrix[,1:3]

locate.fossil(royer_tree, r_fossilrix)






#####temptesting
royer_data_temp <- royer_data_sub[-c(3,5)]
temp_fossil <- all_fossil[1,]
temp_full <- rbind(royer_data_temp, temp_fossil)
temp_full <- temp_full[-c(1)]

temp_full <- as.matrix(temp_full)

fossil_mean <- aggregate(all_fossil[,-1], by=list(all_fossil$binomial), mean)

temptree <- locate.fossil(royer_tree, temp_full, search="exhaustive",plot=TRUE)

add.arrow(tree=temptree,tip="Fagopsis_longifolia",col="red",lwd=3,hedl=0.06,angle=50)


######testingcontall_fossil$binomial
######located.fossil w/ date constraints and edge constrants?
###### should just work
###### ask revell for help
#### get_mrca
##### figtree
tree_quercus <- add.species.to.genus(tree_plant, "Quercus_sp", where=c("root"))
tree_royer_q <- add.species.to.genus(royer_tree, "Quercus_sp", where = c("root"))
tree_fossil <- add.species.to.genus(tree_plant, fossil_mean$Group.1, where=c("root"))



plot(tree_quercus, type = "fan", show.tip.label = FALSE)
plot(tree_royer_q, type = "fan")
add.arrow(tree=tree_royer_q,tip="Quercus_sp",col="blue")

