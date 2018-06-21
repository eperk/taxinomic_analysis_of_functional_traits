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
library(foreach)

royer_phylo_df<-readRDS("./outputs/Royer_clean_df_phylo.rds")

## Function
PEM.royer <- function(sp_name, df, tree){
  
  royer_phylo_rmv <- drop.tip(tree, sp_name)
  royer.train <- df[df$data$binomial!=sp_name,]
  royer.test <- df[df$data$binomial==sp_name,]
  royer_loc <- getGraphLocations(df$phy,
                                 targets = "Piper_reticulatum")
  

  predictions<-list(data=royer.train,phylo=royer_phylo_rmv, royer.test, royer_loc)
  
}


## Execute function
sp_list<-as.character(royer_phylo_df$data$binomial)


#foreach
tmp_obj<-
foreach(i=1:length(sp_list))%do%
{
  
  PEM.royer(sp_name =sp_list[i], df=royer_phylo_df$data, tree=royer_phylo_df$phy)
}

names(tmp_obj)<-sp_list


##for loop

#tmp_obj<-NULL

#for (i in 1:length(sp_list)){
  
 # sp1<-test1(sp_name =sp_list[i], df=royer_phylo_df$data, tree=royer_phylo_df$phy)
  #tmp_obj[[i]] <- sp1
  
  #tmp_obj=rbind(tmp_obj,sp1)
  
}
#names(tmp_obj)<-sp_list


