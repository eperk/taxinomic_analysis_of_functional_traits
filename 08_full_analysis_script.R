#full script for modeling and predicting LMA and phenology from fossil data
#required packages----
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
require(RColorBrewer)

#Initial phylogenetic approach----
#This approach inserts a species of your choice randomly into a phylogeny based on taxonomy and its most recent fossil
#this is a custom function
#not exactly relevent for the lm and glm models but might come in handly

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
intfossil <- function(tree, mintime=0,maxtime=NA, name="fossil", edge=NA, genus="genus", fossil_tax)
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
  if(is.na(edge)){edgesample<-sample(which(maxedge>mintime & minedge<maxtime),1)}
  cull_tree <- drop.clade.label(tree, order)
  cull_tree$edge.length <- drop.tip(tree, cladetree$tip.label, collapse.singles = FALSE)$edge.length
  
  tree_4 <- as(tree, "phylo4")
  reqedge <- getEdge(tree_4, order)
  reqlength <- edgeLength(tree_4)[reqedge]
  length <- tree$edge.length
  reqloc <- as.numeric(which(grepl(reqlength, length)))
  nedge <- as.matrix(cull_tree$edge.length)
  fixedge <- insertRow(nedge, reqloc, reqlength)
  
  cull_tree$edge.length <- fixedge
  
  dedge<-cladetree$edge[edgesample,2]
  place<-runif(1,max(c(minedge[edgesample],mintime)),min(c(maxtime,maxedge[edgesample])))
  #if(is_not_going_to_work)
  #then(dont)
  fossil<-list(edge=matrix(c(2,1),1,2), tip.label=name, edge.length=runif(1,min=0.0000000001,max=(place-max(c(minedge[edgesample],mintime)))), Nnode=1)
  class(fossil)<-"phylo"
  fossiltree<-bind.tree(cladetree,fossil,where=dedge,position=place-minedge[edgesample])
  tree_full <- bind.tree(cull_tree, fossiltree, where = which(cull_tree$tip.label == order))
  tree_full <- force.ultrametric(tree_full, method = "extend")
  return(tree_full)
  
}



#mixed-effects modelling-----
#this section is where we predict lma from petiole^2/la and family
#data creation

royer_tax_full <- readRDS("~/Documents/BEIN data R/data/processed/07_lm4_royer")
royer_fossil_dropped <-   all_fossil_royer_pred %>% 
  filter(!grepl('Cercidiphyllaceae', scrubbed_family)) %>%
  filter(!grepl('Smilacaceae', scrubbed_family)) %>%
  filter(!grepl('Staphyleaceae', scrubbed_family))

royer_tax_fam_count <- royer_tax_full%>% 
  group_by( scrubbed_family) %>% 
  summarize(count=n()) %>% 
  filter(!grepl('Unknown', scrubbed_family))

royer_tax_count <- left_join(royer_tax_full, royer_tax_fam_count, by="scrubbed_family")

#creation of royer_tax_top_5/ royer_tax_over_10---------------
#These are imporatant as these determine the top 5 extant families with the most samples, and the top 10 families respecitively------
royer_tax_top5 <- subset(royer_tax_count, as.numeric(count)>23)
royer_tax_over10 <- subset(royer_tax_count, as.numeric(count)>9)

#creation of the royer_lmfam_top5 and royer_lm_top5, these are the linear models for comparing log lma and contains predictions for log(lma)----
#for 5 fossil famlies Ericaceae, Fabaceae, Fagaceae, Myrtaceae, and Proteaceae
royer_lm_top5 <- lm(log_lma~log_pet_leafarea, data = royer_tax_top5)
pred_lmfam_top5 <-  as.data.frame(predict(royer_lmfam_top5, interval="prediction"))
royer_lmfam_top5_bound <-  cbind(royer_tax_top5, pred_lmfam_top5)


#The following two sections plot log_lma vs log(petiole^2/leafarea) for the top 5 families---------
ggplot(royer_lmfam_top5_bound, aes(log_pet_leafarea, log_lma))+
  aes(color=royer_lmfam_top5_bound$scrubbed_family)+
  geom_point() +
  geom_smooth(method=lm, se=TRUE)+
  labs(x="log(petiole^2/leafarea)")+
  labs(y="log(lma)")+
  labs(color = "Family")

ggplot(royer_lmfam_top5_bound, aes(log_pet_leafarea, log_lma))+
  geom_point() +
  geom_smooth(method=lm, se=TRUE)+
  labs(x="log(petiole^2/leafarea)")+
  labs(y="log(lma)")+
  labs(color = "Family")

##This is the logisitic regression predicting phenlogoy off of log_lma-----
ggplot(glm_pred_data, aes(log_lma, phenology))+ 
  geom_point() + 
  stat_smooth(method="glm", method.args=list(family="binomial"), se=FALSE)

ggplot(glm_pred_data, aes(log_lma, phenology, group=(scrubbed_family), color=scrubbed_family)) + 
  geom_point() + 
  stat_smooth(method="glm", method.args=list(family="binomial"), se=FALSE)




