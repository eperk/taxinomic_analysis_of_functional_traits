require(tidyr)
require(dplyr)
require(phangorn)
require(ape)
require(phytools)
require(phylolm)


######fabaceae fastanc analysis
WSLA.tree.nodes <- tree.WSLA.species$node.label
WLSA.tree.nodes.df <- as.data.frame(WSLA.tree.nodes)
fabaceae <- as.vector("Fabaceae")
tree.fabaceae <- extract.clade(tree.WSLA.species, "Fabaceae")
fabaceae.tips<- tree.fabaceae$tip.label
fabaceae.SLA.mean.binomials <- subset(SLA.only.means, SLA.only.means$binomial %in% fabaceae.tips)
fabaceae.labels <- fabaceae.SLA.mean.binomials[,-1]
fabaceae.labels <- as.data.frame(fabaceae.labels)
rownames(fabaceae.labels) <- fabaceae.SLA.mean.binomials[,1]
fabaceae.matrix <-as.matrix(fabaceae.labels)[,1]
fabaceae.fastAnc<-fastAnc(tree.fabaceae,fabaceae.matrix,vars=TRUE,CI=TRUE)
fabaceae.mapping<-contMap(tree.fabaceae,fabaceae.matrix,plot=FALSE)
plot(fabaceae.mapping,legend=0.7*max(nodeHeights(tree.fabaceae)), fsize=c(0.7,0.9), ftype="off")
plot(tree.fabaceae, type="fan", show.tip.label = FALSE)


#######simmap fabaceae
species.state.fabaceae <- subset(deomited.states, (deomited.states$binomial %in% fabaceae.tips))
species.pheno.fabaceae<-species.state.fabaceae$Phenology

####simmap tree creation 
not.in.fabaceae <- subset(fabaceae.tips, !(fabaceae.tips %in% species.state.fabaceae$binomial))
not.in.fabaceae.df <-as.data.frame(not.in.fabaceae)
tree.simmap.fabaceae <- drop.tip(tree.fabaceae, not.in.fabaceae)
names(species.pheno.fabaceae) <- species.state.fabaceae$binomial

###########
fabaceae.state.matrix <-as.matrix(species.state.fabaceae)[,2]
fabaceae.simmap <- make.simmap(tree.simmap.fabaceae, species.pheno.fabaceae, nsim=1)
plotSimmap(fabaceae.simmap, ftype="off")




################
######Myrtales fastanc analysis
WSLA.tree.nodes <- tree.WSLA.species$node.label
Myrtales <- as.vector("Myrtales")
tree.myrtales <- extract.clade(tree.WSLA.species, "Myrtales")
myrtales.tips<- tree.myrtales$tip.label
myrtales.SLA.mean.binomials <- subset(SLA.only.means, SLA.only.means$binomial %in% myrtales.tips)
myrtales.labels <- myrtales.SLA.mean.binomials[,-1]
myrtales.labels <- as.data.frame(myrtales.labels)
rownames(myrtales.labels) <- myrtales.SLA.mean.binomials[,1]
myrtales.matrix <-as.matrix(myrtales.labels)[,1]
myrtales.fastAnc<-fastAnc(tree.myrtales,myrtales.matrix,vars=TRUE,CI=TRUE)
myrtales.mapping<-contMap(tree.myrtales,myrtales.matrix,plot=FALSE)
plot(myrtales.mapping,legend=0.7*max(nodeHeights(tree.myrtales)), fsize=c(0.7,0.9), ftype="off")


#######simmap
species.state <- final_WLSA_DF [-c(2,3,5,6)]
deduped.states <- unique( species.state[ , 1:2 ] )
deomited.states <- na.omit(deduped.states)
species.state.myrtales <- subset(deomited.states, (deomited.states$binomial %in% myrtales.tips))
species.pheno.myrtales<-species.state.myrtales$Phenology

####simmap tree creation something broke between wednesday and sunday and I cant figure it out
not.in.myrtales <- subset(myrtales.tips, !(myrtales.tips %in% species.state.myrtales$binomial))
not.in.myrtales.df <-as.data.frame(not.in.myrtales)
tree.simmap.myrtales <- drop.tip(tree.myrtales, not.in.myrtales)
names(species.pheno.myrtales) <- species.state.myrtales$binomial

###########myrtales.state.matrix <-as.matrix(species.state.myrtales)[,2]
myrtales.simmap <- make.simmap(tree.simmap.myrtales, species.pheno.myrtales, nsim=1)
plotSimmap(myrtales.simmap)

