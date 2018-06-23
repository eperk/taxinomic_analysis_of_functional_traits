#fossil matrix creation-------------------------------------------------------------------------------------------------
r_fossilrix <- as.matrix(royer_data_fossil_int)
r_fossilrix <- r_fossilrix[,1:3]

#fossil location N.B. broken-------------------------------------------------------------------------------------------------
#locate.fossil(royer_tree, r_fossilrix)

#ppgm funciton edit-------------------------------------------------------------------------------------------------

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

