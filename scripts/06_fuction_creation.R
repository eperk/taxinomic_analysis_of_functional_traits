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

#intfossil function-----------------------------------------
intfossil <- function(tree, mintime=0,maxtime=NA, name="fossil", edge=NA, genus="genus")
  {
  require(ape)
  lookup <- match(genus, fossil_tax$scrubbed_genus)
  taxonomy <- fossil_tax[na.omit(lookup), ]
  order <- taxonomy$order
  cladetree <- extract.clade(tree, order)
  if(is.na(maxtime)){maxtime=max(dist.nodes(cladetree))/2}
  treeage<-max(dist.nodes(cladetree))/2
  M<-dist.nodes(cladetree)
  maxedge<-(as.numeric(treeage - M[cladetree$edge[,1],cladetree$edge[1,1]]))
  minedge<-(as.numeric(treeage - M[cladetree$edge[,2],cladetree$edge[1,1]]))
  if(!is.na(edge)){edgesample<-edge}
  if(is.na(edge)){edgesample<-sample(which(maxedge>mintime & minedge<maxtime),1)}
  #dropclade-----
  clade <- cladetree$tip.label
  cull_tree <- drop.clade(tree_plant, clade)
  #cladetree$node.label<-((length(cladetree$tip)+1):((length(cladetree$tip)*2)-1))
  dedge<-cladetree$edge[edgesample,2]
  place<-runif(1,max(c(minedge[edgesample],mintime)),min(c(maxtime,maxedge[edgesample])))
  #removed to see if bindtip works better--------------
  fossil<-list(edge=matrix(c(2,1),1,2), tip.label=name, edge.length=runif(1,min=0.0000000001,max=(place-max(c(minedge[edgesample],mintime)))), Nnode=1)
  class(fossil)<-"phylo"
  tree<-bind.tree(cladetree,fossil,where=dedge,position=place-minedge[edgesample])
  #cladetree$node.label<-as.numeric(tree$node.label)+1
  #newnode=which(is.na(tree$node.label))
  #cladetree$node.label[(newnode+1):length(tree$node.label)]<-as.numeric(tree$node.label[(newnode+1):length(tree$node.label)])+1
  #cladetree$node.label[newnode]<-as.numeric(tree$node.label[newnode-1])+1
  #insertclade busted here to fix, label node that the tree is cut at with an NA, so it should work?-----
  tree_full <- bind.tree(tree, cull_tree, where = which(tree$node.label == order))
  return(tree_full)
  #return (tree)
  
}
#testing below-------
rosidae_insert <- intfossil(tree_plant, mintime = 0, maxtime = 33900000, name = "Rosa sp.", edge = NA, genus = "Rosa")
plot(rosidae_insert, show.tip.label = FALSE, type = "fan")
add.arrow(tree = rosidae_insert, tip = "Rosa sp.", col="red",lwd=3,hedl=0.06,angle=90)

prunus_insert <- intfossil(tree_plant, mintime = 0, maxtime = 33900000, name = "Prunus_scottii", edge = NA, genus = "Prunus")
plot(prunus_insert, type="fan", show.tip.label=FALSE)
add.arrow(tree = prunus_insert, tip = "Prunus_scottii", col="red",lwd=3,hedl=0.06,angle=90)

hydrangea_insert <- intfossil(tree_plant, mintime = 0, maxtime = 33900000, name = "hydrangea sp.", edge = NA, genus = "Hydrangea")
plot(hydrangea_insert, type="fan", show.tip.label=FALSE)
add.arrow(tree = hydrangea_insert, tip = "hydrangea sp.", col="red",lwd=3,hedl=0.06,angle=90)


#full manual running of function THERE IS NO REASON FOR THIS TO BE DIFFERENT FROM THE OTHER ONE BUT IT IS ANYWAY------
maxtime <- 33900000
mintime <- 0
lookup <- match("Rosa", fossil_tax$scrubbed_genus)
taxonomy <- fossil_tax[na.omit(lookup), ]
order <- taxonomy$order
cladetree <- extract.clade(tree_plant, order)
M<-dist.nodes(cladetree)
treeage<-max(dist.nodes(cladetree))/2
maxedge<-(as.numeric(treeage - M[cladetree$edge[,1],cladetree$edge[1,1]]))
minedge<-(as.numeric(treeage - M[cladetree$edge[,2],cladetree$edge[1,1]]))
edgesample<-sample(which(maxedge>mintime & minedge<maxtime),1)
#dropclade-----
dedge<-cladetree$edge[edgesample,2]
place<-runif(1,max(c(minedge[edgesample],mintime)),min(c(maxtime,maxedge[edgesample])))
fossil<-list(edge=matrix(c(2,1),1,2), tip.label="rosa sp.", edge.length=runif(1,min=0.0000000001,max=(place-max(c(minedge[edgesample],mintime)))), Nnode=1)
class(fossil)<-"phylo"
tree_temptest<-bind.tree(cladetree,fossil,where=dedge,position=place-minedge[edgesample])
dropclade <- cladetree$tip.label
plant_sisters_temp <- getSisters(tree_plant,"Rosales", mode = c("number"))
cull_tree <- drop.clade(tree_plant, dropclade)
#note broken, need to figure out how to retain tree when creating culltree-------
tree_bind_test <- bind.tree(cull_tree, tree_temptest, where = which(tree_plant$node.label == order))
plot(tree_bind_test, type="fan", show.tip.label=FALSE)
add.arrow(tree = hydrangea_insert, tip = "prunus sp.", col="red",lwd=3,hedl=0.06,angle=90)


#tree_temp <- add.species.to.genus(tree_rosidae, "Rosa sp.", genus=NULL, where = c("root"))

##for loop

#tmp_obj<-NULL

for (i in 1:length(all_fossil_phyloint))
  {
  
  sp1<-test1(sp_name =sp_list[i], df=royer_phylo_df$data, tree=royer_phylo_df$phy)
  tmp_obj[[i]] <- sp1
  
  tmp_obj=rbind(tmp_obj,sp1)
  
names(tmp_obj)<-sp_list

}
