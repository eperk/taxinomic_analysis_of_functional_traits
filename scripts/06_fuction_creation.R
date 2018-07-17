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
#extract.clade.label-----
#depedency within drop.clade.label------
#https://github.com/MarioJose/r-functions
extract.clade.label <- function(tree, node){
  if(!is.vector(node, mode = "character") | length(node) > 1)
    stop("'node' parameter must be a character vector of length 1")
  
  if(sum(node %in% tree$node.label) == 0)
    stop("tree has not node labels defined in 'node' parameter")
  
  if(!is.rooted(tree))
    stop("tree must be rooted")
  
  # Solve bug reading newick in phytools 0.6
  if(tree$Nnode > length(tree$node.label))
    tree$node.label <- c(tree$node.label, rep('', tree$Nnode - length(tree$node.label)))
  
  # Reorder tree
  #tree <- reorder(tree)
  
  # number of tips
  not <- length(tree$tip.label)
  # node number
  nn <- which(tree$node.label %in% node) + not
  
  edge <- tree$edge
  
  # nodes descendants
  nd <- tmp <- edge[edge[ ,1] %in% nn , 2]
  while(length(tmp) > 0){
    tmp <- edge[edge[ ,1] %in% tmp , 2]
    if(length(tmp) > 0) nd <- c(nd, tmp)
  }
  # tips desdendants
  td <- nd[nd <= not]
  # only nodes descendats
  nd <- nd[nd > not] 
  
  # New tips positions
  newt <- matrix(NA, nrow = not, ncol = 2, 
                 dimnames = list(1:not, c('tips', 'new')))
  newt[ ,1] <- 1:not
  newt <- newt[newt[ ,1] %in% td, ]
  newt[ ,2] <- rank(newt[ ,1])
  
  # New nodes positions
  newn <- matrix(NA, nrow = tree$Nnode, ncol = 2, 
                 dimnames = list(1:tree$Nnode, c('nodes', 'new')))
  newn [ ,1] <- 1:tree$Nnode + not
  newn <- newn[newn[ ,1] %in% c(nd, nn), ]
  if(!is.matrix(newn))
    newn <- matrix(newn, nrow = 1, dimnames = list(1, names(newn)))
  newn[ ,2] <- rank(newn[ ,1]) + length(td)
  
  # Reference table to new nodes positions
  refedge <- rbind(newt, newn)
  
  newedge <- edge[edge[ ,1] %in% c(nd, nn), ]
  newedge[ ,1] <- refedge[match(newedge[ ,1], refedge[ ,1]) ,2]
  newedge[ ,2] <- refedge[match(newedge[ ,2], refedge[ ,1]) ,2]
  
  # Update tree  
  tree$edge <- newedge
  tree$tip.label <- tree$tip.label[newt[ ,1]]
  tree$node.label <- tree$node.label[newn[ ,1] - not]
  tree$Nnode <- length(tree$node.label)
  
  return(tree)
}
#drop.clade.label function, dependency for intfossil---------
#https://github.com/MarioJose/r-functions
drop.clade.label <- function(tree, node){
  # Remove all tips from 'node', keeping 'node' label as a tip
  
  if(!is.vector(node, mode = "character") | length(node) > 1)
    stop("'node' parameter must be a character vector of length 1")
  
  if(sum(node %in% tree$node.label) == 0)
    stop("tree has not node labels defined in 'node' parameter")
  
  if(!is.rooted(tree))
    stop("tree must be rooted")
  
  # Solve bug reading newick in phytools 0.6
  if(tree$Nnode > length(tree$node.label))
    tree$node.label <- c(tree$node.label, rep('', tree$Nnode - length(tree$node.label)))
  
  # Reorder tree
  #tree <- reorder(tree)
  
  # number of tips
  not <- length(tree$tip.label)
  # node number
  nn <- which(tree$node.label %in% node) + not
  
  # 'node' tree extracted
  extr <- extract.clade.label(tree, node)
  # number of tips of extracted 'node' tree
  extrnot <- length(extr$tip.label)
  
  # new number of tips
  nnot <- not - extrnot
  
  edge <- tree$edge
  
  # nodes descendants
  nd <- tmp <- edge[edge[ ,1] %in% nn , 2]
  while(length(tmp) > 0){
    tmp <- edge[edge[ ,1] %in% tmp , 2]
    if(length(tmp) > 0) nd <- c(nd, tmp)
  }
  # tips desdendants
  td <- nd[nd <= not]
  # only nodes descendats
  nd <- nd[nd > not] 
  
  # New tips positions
  newt <- matrix(NA, nrow = not, ncol = 2, 
                 dimnames = list(1:not, c('tips', 'new')))
  newt[ ,1] <- 1:not
  newt <- newt[!(newt[ ,1] %in% td), ]
  newt[ ,2] <- rank(newt[ ,1])
  
  # New nodes positions
  newn <- matrix(NA, nrow = tree$Nnode, ncol = 2, 
                 dimnames = list(1:tree$Nnode, c('nodes', 'new')))
  newn [ ,1] <- 1:tree$Nnode + not
  newn <- newn[!(newn[ ,1] %in% c(nd, nn)), ]
  # 'nnot' plus 1L to new tip 'node'
  newn [ ,2] <- rank(newn[ ,1]) + nnot + 1L
  
  # Reference table to new nodes positions
  refedge <- rbind(newt, newn, matrix(c(nn, (nnot + 1)), nrow = 1))
  
  newedge <- edge[!(edge[ ,1] %in% c(nd, nn)), ]
  newedge[ ,1] <- refedge[match(newedge[ ,1], refedge[ ,1]) ,2]
  newedge[ ,2] <- refedge[match(newedge[ ,2], refedge[ ,1]) ,2]
  
  # Update tree  
  tree$edge <- newedge
  tree$tip.label <- c(tree$tip.label[newt[ ,1]], node)
  tree$node.label <- tree$node.label[newn[ ,1] - not]
  tree$Nnode <- length(tree$node.label)
  
  return(tree)
}

#extract.node.label and drop.node.lable fromp------
#https://github.com/MarioJose/r-functions/blob/master/drop.clade.label/extract.clade.label.r

#intfossil function-----------------------------------------
#created through addfossil function https://github.com/michellelawing/ppgm
#dependencies: ape, phytools, drop.clade.labels
intfossil <- function(tree, mintime=0,maxtime=NA, name="fossil", edge=NA, genus="genus")
  {
  require(ape)
  lookup <- match(genus, fossil_tax$scrubbed_genus)
  taxonomy <- fossil_tax[na.omit(lookup), ]
  order <- taxonomy$order
  cladetree <- extract.clade.label(tree, order)
  cladetree$edge.length<- extract.clade(tree, order)$edge.length
  if(is.na(maxtime)){maxtime=max(dist.nodes(cladetree))/2}
  treeage<-max(dist.nodes(cladetree))/2
  M<-dist.nodes(cladetree)
  maxedge<-(as.numeric(treeage - M[cladetree$edge[,1],cladetree$edge[1,1]]))
  minedge<-(as.numeric(treeage - M[cladetree$edge[,2],cladetree$edge[1,1]]))
  if(!is.na(edge)){edgesample<-edge}
  if(is.na(edge)){edgesample<-sample(which(maxedge>mintime & minedge<maxtime),1)}
  #dropclade-----
  clade <- cladetree$tip.label
  cull_tree <- drop.clade.label(tree, order)
  cull_tree$edge.length <- drop.tip(tree, clade, collapse.singles = FALSE)$edge.length
  #cladetree$node.label<-((length(cladetree$tip)+1):((length(cladetree$tip)*2)-1))
  dedge<-cladetree$edge[edgesample,2]
  place<-runif(1,max(c(minedge[edgesample],mintime)),min(c(maxtime,maxedge[edgesample])))
  #removed to see if bindtip works better--------------
  #if(is_not_going_to_work)
  #then(dont)
  fossil<-list(edge=matrix(c(2,1),1,2), tip.label=name, edge.length=runif(1,min=0.0000000001,max=(place-max(c(minedge[edgesample],mintime)))), Nnode=1)
  class(fossil)<-"phylo"
  fossiltree<-bind.tree(cladetree,fossil,where=dedge,position=place-minedge[edgesample])
  
  edgenum <- as.numeric(nrow(cull_tree$edge))
  lengthmatrix <- as.matrix(tree$edge.length)
  missinglength <- as.matrix(lengthmatrix[edgenum,])
  nedge <- rbind(as.matrix(cull_tree$edge.length), missinglength)
  cull_tree$edge.length <- nedge
  
  #cladetree$node.label<-as.numeric(tree$node.label)+1
  #newnode=which(is.na(tree$node.label))
  #cladetree$node.label[(newnode+1):length(tree$node.label)]<-as.numeric(tree$node.label[(newnode+1):length(tree$node.label)])+1
  #cladetree$node.label[newnode]<-as.numeric(tree$node.label[newnode-1])+1-----
  tree_full <- bind.tree(cull_tree, fossiltree, where = which(cull_tree$tip.label == order))
  return(tree_full)

}
#testing below--------------------------------------------------------------------------------------------------
rosidae_insert <- intfossil(tree_plant, mintime = 0, maxtime = 33900000, name = "Rosa sp.", edge = NA, genus = "Rosa")
plot(rosidae_insert, show.tip.label = FALSE)
add.arrow(tree = rosidae_insert, tip = "Rosa sp.", col="red",lwd=3,hedl=0.06,angle=90)

prunus_insert <- intfossil(rosidae_insert, mintime = 0, maxtime = 33900000, name = "Prunus_scottii", edge = NA, genus = "Prunus")
plot(prunus_insert, type="fan", show.tip.label=FALSE)
add.arrow(tree = prunus_insert, tip = "Prunus_scottii", col="red",lwd=3,hedl=0.06,angle=90)
add.arrow(tree = prunus_insert, tip = "Rosa sp.", col="blue",lwd=3,hedl=0.06,angle=90)


hydrangea_insert <- intfossil(prunus_insert, mintime = 0, maxtime = 33900000, name = "hydrangea sp.", edge = NA, genus = "Hydrangea")
plot(hydrangea_insert, type="fan", show.tip.label=FALSE)
add.arrow(tree = hydrangea_insert, tip = "hydrangea sp.", col="green",lwd=3,hedl=0.06,angle=90)
add.arrow(tree = hydrangea_insert, tip = "Prunus_scottii", col="red",lwd=3,hedl=0.06,angle=90)
add.arrow(tree = hydrangea_insert, tip = "Rosa sp.", col="blue",lwd=3,hedl=0.06,angle=90)

#full manual running of function, now identical, use for contiuned function testing---------
maxtime <- 33900000
mintime <- 0
lookup <- match("Rosa", fossil_tax$scrubbed_genus)
taxonomy <- fossil_tax[na.omit(lookup), ]
order <- taxonomy$order

#renumbers edges works fine?-------
cull_tree <- drop.clade.label(tree_plant, order)
cladetree <- extract.clade.label(tree_plant, order)
cladetree$edge.length<- extract.clade(tree_plant, order)$edge.length
cull_temp_edge <- drop.tip(tree_plant, cladetree$tip.label,  collapse.singles = FALSE)
cull_tree$edge.length <- cull_temp_edge$edge.length
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

#CRITICAL ADDS FIXED EDGES------

edgenum <- as.numeric(nrow(cull_tree$edge))
lengthmatrix <- as.matrix(tree_plant$edge.length)
missinglength <- as.matrix(lengthmatrix[edgenum,])
nedge <- rbind(as.matrix(cull_tree$edge.length), missinglength)
cull_tree$edge.length <- nedge

tree_bind_test <- bind.tree(cull_tree, tree_temptest, where = which(cull_tree$tip.label == order))
plot(tree_bind_test, show.tip.label=FALSE)
