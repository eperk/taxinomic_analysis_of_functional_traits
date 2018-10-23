#fossil matrix creation-------------------------------------------------------------------------------------------------
r_fossilrix <- as.matrix(royer_data_fossil_int)
r_fossilrix <- r_fossilrix[,1:3]

#fossil location N.B. broken-------------------------------------------------------------------------------------------------
#locate.fossil(royer_tree, r_fossilrix)

#addfunction function usage-------------------------------------------------------------------------------------------------

fossil_tree <- addfossil(tree_rosidae, mintime = 0, maxtime = NA, name = "rosidae sp.", edge = NA)
plot(fossil_tree, show.tip.label = FALSE, type = "fan")
add.arrow(tree = fossil_tree, tip = "rosidae sp.", col="red",lwd=3,hedl=0.06,angle=50)

          
#temptesting-------------------------------------------------------------------------------------------------
royer_data_temp <- royer_data_sub[-c(3,5)]
temp_fossil <- all_fossil[1,]
temp_full <- rbind(royer_data_temp, temp_fossil)
temp_full <- temp_full[-c(1)]

temp_full <- as.matrix(temp_full)

fossil_mean <- aggregate(all_fossil[,-1], by=list(all_fossil$binomial), mean)

temptree <- locate.fossil(royer_tree, temp_full, search="exhaustive",plot=TRUE)

add.arrow(tree=temptree,tip="Fagopsis_longifolia",col="red",lwd=3,hedl=0.06,angle=50)

similax_insert <- intfossil(tree_plant, mintime = 0, maxtime = 33900000, name = "Smilax sp.", edge = NA, genus = "Smilax")
plot(similax_insert, type="fan", show.tip.label=FALSE)
add.arrow(tree = prunus_insert, tip = "Smiliax sp.", col="red",lwd=3,hedl=0.06,angle=90)


######testingcontall_fossil$binomial
######located.fossil w/ date constraints and edge constrants?
###### should just work
###### ask revell for help
#### get_mrca
##### figtree
###taxize
#use ppgm functions to conduct my phlyo testing - see notebook

#add species to genus testing-------------------------------------------------------------------------------------------------
tree_quercus <- add.species.to.genus(tree_plant, "Quercus_sp", where=c("root"))
tree_royer_q <- add.species.to.genus(royer_tree, "Quercus_sp", where = c("root"))
tree_fossil <- add.species.to.genus(tree_plant, fossil_mean$Group.1, where=c("root"))



plot(tree_quercus, type = "fan", show.tip.label = FALSE)
plot(tree_royer_q, type = "fan")
add.arrow(tree=tree_royer_q,tip="Quercus_sp",col="blue")

#creation of dataframe to be run through for loop----------
all_fossil_phyloint<- readRDS("~/Documents/BEIN data R/data/processed/04_fossil_phylo_integration.rds")
all_fossil_phyloint_match<-subset(all_fossil_phyloint, 
                  !(all_fossil_phyloint$Genus %in% missing_genus$Genus))
all_fossil_phyloint_match$date <- as.numeric(all_fossil_phyloint_match$date)



#####current for loop, still has issues---------------


  for(i in 1:nrow(all_fossil_phyloint_match))
  {
    intfossil(
      tree = tree_plant, 
      mintime = 0,
      maxtime = as.numeric(all_fossil_phyloint_match$date[i]), 
      name = all_fossil_phyloint_match$binomial[i], 
      edge = NA,  
      genus = all_fossil_phyloint_match$Genus[i],
      fossil_tax = fossil_tax)
  }



 foreach(i:nrow(all_fossil_phyloint_match)) %do%
 {
   intfossil(
     tree = tree_plant, 
     mintime = 0,
     maxtime = as.numeric(all_fossil_phyloint_match$date[i]), 
     name = all_fossil_phyloint_match$binomial[i], 
     edge = NA,  
     genus = all_fossil_phyloint_match$Genus[i]
     
   )
 }
   


intfossil_test <- intfossil(
  tree = tree_plant, 
  mintime = 0,
  maxtime = 33900000, 
  name = "Koelreuteria_sp.", 
  edge = NA,  
  genus = "Koelreuteria")


