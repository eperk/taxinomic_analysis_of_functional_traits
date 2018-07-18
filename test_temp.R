tr <- extract.clade(tree_plant, "Acorales")
tr$node.label <- (1:6)
plot(tr)
tiplabels()
nodelabels()
edgelabels(tr$edge.length)

edgelookup <- cbind(tr$edge,tr$edge.length)

specficnode<-sapply("4",function(x,y) which(y==x),y="4")
## first get the node numbers of the tips
nodes<-sapply("Acorales",function(x,y) which(y==x),y="Acorales")
## then get the edge lengths for those nodes
edge.lengths<-setNames(tr$edge.length[sapply(nodes,
  function(x,y) which(y==x),y=tr$edge[,2])],names(nodes))

tr_4 <- as(tr, "phylo4")
plot(tr_4)
getEdge(tr_4, "4")
edgeLength(tr_4)["10-11"]


reqedge <- (tr)


nodelabels <- as.data.frame(tr$node.label)
colnames(nodelabels)[colnames(nodelabels)=="tr$node.label"] <- "nodes"
getnode <- which(nodelabels$nodes == 4)


node <- 11
nodes<-sapply(node,function(x,y) which(y==x),y=tr$node.label)
edge.lengths<-setNames(tr$edge.length[sapply(nodes, function(x,y) which(y==x),y=tr$edge[,2])],names(nodes))

#THIS IS SO GODDAMN IMPORTANT---------------
#nevermind, doesnt work

edgenum <- as.numeric(nrow(tr_dr2$edge))
lengthmatrix <- as.matrix(tr$edge.length)
missinglength <- as.matrix(lengthmatrix[edgenum,])
nedge <- rbind(as.matrix(tr_dr2$edge.length), missinglength)
tr_dr2$edge.length <- nedge


tr_cull <- extract.clade.label(tr, "2")
plot(collapse.singles(tr_cull), node.depth = 2)
snodelabels(collapse.singles(tr_cull)$node.label, cex=0.8, adj = 0.5)
tr_cull$edge.length<- extract.clade(tr, "2")$edge.length


tr_dr2 <- drop.clade.label(tr, "2")
plot(collapse.singles(tr_dr2), node.depth = 2)
nodelabels(collapse.singles(tr_dr2)$node.label, cex=0.8, adj = 0.5)
tr_test <- drop.tip(tr, tr_cull$tip.label,  collapse.singles = FALSE)
tr_dr2$edge.length <- fixedge

tr_4 <- as(tr, "phylo4")
plot(tr_4)
reqedge <- getEdge(tr_4, "2")
reqlength <- edgeLength(tr_4)[reqedge]
reqloc <- as.numeric(which(grepl(reqlength, lenght_temp)))

nedge <- as.matrix(tr_dr2$edge.length)

fixedge <- insertRow(nedge, reqloc, reqlength)

tr_test <- drop.tip(tr, tr_cull$tip.label,  collapse.singles = FALSE)
plot(tr_test)
nodelabels(collapse.singles(tr_test)$node.label, cex=0.8, adj = 0.5)
plot(tr_test)



tr_full_temp <- bind.tree(tr_dr2,tr_cull, where = which(tr_dr2$tip.label == "2"))
par(mar=c(1,1,1,1))
plot(tr_full_temp)
nodelabels(collapse.singles(tr_full_temp)$node.label, cex=0.8, adj = 0.5)



###these functions, when used in conjunction work perfectly!! lit