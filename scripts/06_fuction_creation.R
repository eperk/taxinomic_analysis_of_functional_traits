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


#addfossil function---------------------------------
addfossil<- function(tree,mintime=0,maxtime=NA,name="fossil",edge=NA) {
  #function depends on ape
  require(ape)
  if(is.na(maxtime)){maxtime=max(dist.nodes(tree))/2}
  tree$node.label<-((length(tree$tip)+1):((length(tree$tip)*2)-1))
  treeage<-max(dist.nodes(tree))/2
  M<-dist.nodes(tree)
  maxedge<-(as.numeric(treeage - M[tree$edge[,1],tree$edge[1,1]]))
  minedge<-(as.numeric(treeage - M[tree$edge[,2],tree$edge[1,1]]))
  if(!is.na(edge)){edgesample<-edge}
  if(is.na(edge)){edgesample<-sample(which(maxedge>mintime & minedge<maxtime),1)}
  dedge<-tree$edge[edgesample,2]
  place<-runif(1,max(c(minedge[edgesample],mintime)),min(c(maxtime,maxedge[edgesample])))
  fossil<-list(edge=matrix(c(2,1),1,2), tip.label=name, edge.length=runif(1,min=0.0000000001,max=(place-max(c(minedge[edgesample],mintime)))), Nnode=1)
  class(fossil)<-"phylo"
  tree<-bind.tree(tree,fossil,where=dedge,position=place-minedge[edgesample])
  tree$node.label<-as.numeric(tree$node.label)+1
  newnode=which(is.na(tree$node.label))
  tree$node.label[(newnode+1):length(tree$node.label)]<-as.numeric(tree$node.label[(newnode+1):length(tree$node.label)])+1
  tree$node.label[newnode]<-as.numeric(tree$node.label[newnode-1])+1
  return(tree)
}
#intfossil function-----------------------------------------
intfossil <- function(tree, mintime=0,maxtime=NA, name="fossil", edge=NA, genus="genus")
  {
  require(ape)
  lookup <- match(genus, fossil_tax$scrubbed_genus)
  taxonomy <- fossil_tax[na.omit(lookup), ]
  cladetree <- extract.clade(tree, taxonomy$order)
  
  if(is.na(maxtime)){maxtime=max(dist.nodes(cladetree))/2}
  cladetree$node.label<-((length(cladetree$tip)+1):((length(cladetree$tip)*2)-1))
  treeage<-max(dist.nodes(cladetree))/2
  M<-dist.nodes(cladetree)
  maxedge<-(as.numeric(treeage - M[cladetree$edge[,1],cladetree$edge[1,1]]))
  minedge<-(as.numeric(treeage - M[cladetree$edge[,2],cladetree$edge[1,1]]))
  if(!is.na(edge)){edgesample<-edge}
  if(is.na(edge)){edgesample<-sample(which(maxedge>mintime & minedge<maxtime),1)}
  #dropclade-----
  clade <- cladetree$tip.label
  cull_tree <- drop.tip(tree_plant, clade, trim.internal = TRUE, collapse.singles = TRUE)
  dedge<-cladetree$edge[edgesample,2]
  place<-runif(1,max(c(minedge[edgesample],mintime)),min(c(maxtime,maxedge[edgesample])))
  fossil<-list(edge=matrix(c(2,1),1,2), tip.label=name, edge.length=runif(1,min=0.0000000001,max=(place-max(c(minedge[edgesample],mintime)))), Nnode=1)
  class(fossil)<-"phylo"
  tree<-bind.tree(cladetree,fossil,where=dedge,position=place-minedge[edgesample])
  cladetree$node.label<-as.numeric(tree$node.label)+1
  newnode=which(is.na(tree$node.label))
  cladetree$node.label[(newnode+1):length(tree$node.label)]<-as.numeric(tree$node.label[(newnode+1):length(tree$node.label)])+1
  cladetree$node.label[newnode]<-as.numeric(tree$node.label[newnode-1])+1
  #insertclade busted here to fix, label node that the tree is cut at with an NA, so it should work?-----
  tree_full <- bind.tree(tree, cull_tree, where = which(tree$node.label == genus))
  #return(tree_full)
  return (tree)
  
}
#testing below-------
rosidae_insert <- intfossil(tree_plant, mintime = 0, maxtime = 35000000, name = "Rosa sp.", edge = NA, genus = "Rosa")
plot(rosidae_insert, show.tip.label = FALSE, type = "fan")
add.arrow(tree = rosidae_insert, tip = "Rosa sp.", col="red",lwd=3,hedl=0.06,angle=90)

tree_bind_test <- bind.tree(tree_plant, tree_rosidae, where = which(tree_plant$node.label == "Rosales"))

cull_tree <- drop.tip(tree_plant, rosidae_tips, trim.internal = TRUE, collapse.singles = TRUE)

drop_clade_temp <- rosidea_tips

#tree_temp <- add.species.to.genus(tree_rosidae, "Rosa sp.", genus=NULL, where = c("root"))

##for loop

#tmp_obj<-NULL

#for (i in 1:length(sp_list)){
  
 # sp1<-test1(sp_name =sp_list[i], df=royer_phylo_df$data, tree=royer_phylo_df$phy)
  #tmp_obj[[i]] <- sp1
  
  #tmp_obj=rbind(tmp_obj,sp1)
  
#names(tmp_obj)<-sp_list


