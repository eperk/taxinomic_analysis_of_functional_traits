#PEM analysis-------------------------------------------------------------------------
#creation of all components required for PEM analysis-----------

royer_data <- read_csv("~/Documents/BEIN data R/data/raw/royer_data.csv")

royer_data$binomial <- str_replace_all(royer_data$binomial,"\\s+","_")

royer_data_sub<-
  royer_data %>% 
  dplyr::filter(binomial%in%tree_plant$tip.label) %>% 
  group_by(binomial) %>% 
  summarise(avg_petiole_width=mean(petiole_width, na.rm=TRUE),
            avg_LMA=mean(LMA, na.rm=TRUE), avg_LA=mean(leaf_area)) %>% 
            as.data.frame()

royer_data_sub$log_lma <- log(royer_data_sub$avg_LMA)
royer_data_sub$log_pet_leafarea <- log(royer_data_sub$avg_petiole_width^2/royer_data_sub$avg_LA)

rownames(royer_data_sub)<-royer_data_sub$binomial

indx<-which(tree_plant$tip.label%in%royer_data_sub$binomial==FALSE)

royer_tree <- drop.tip(tree_plant, tree_plant$tip.label[indx])
royer_tips <- royer_tree$tip.label


spmatch <- match(royer_tree$tip.label,
                 royer_data_sub[,1L])
royer_match_data <- royer_data_sub[spmatch,]


royer_phylo_lma<-match.phylo.data(royer_tree,royer_data_sub)
royer_phylo_lma$data$avg_petiole_width<-as.numeric(royer_phylo_lma$data$avg_petiole_width)
royer_phylo_lma$data$avg_LMA<-as.numeric(royer_phylo_lma$data$avg_LMA)
royer_phylo_lma$data$log_lma<-as.numeric(royer_phylo_lma$data$log_lma)
royer_phylo_lma$data$log_pet_leafarea<-as.numeric(royer_phylo_lma$data$log_pet_leafarea)
write_rds(royer_phylo_lma, "./outputs/Royer_clean_df_phylo.rds")

royer_train <- royer_match_data[-(1:2),,drop=FALSE]
royer_phylo_rmv <- drop.tip(royer_phylo_lma$phy, c("Piper_reticulatum","Piper_arboreum"))
royer_test <- royer_match_data[1:2,,drop=FALSE]

royer_lm <-  glm(data = royer_data_sub, log_lma~log_pet_leafarea)
lm_pred_royer <- predict.glm(royer_lm, newdata = royer_test, se.fit = TRUE)%>%
  as.data.frame()


royer_phy <- royer_phylo_lma$phy
royer_phy$node.label <- paste("N",1L:royer_phy$Nnode)
plot(royer_phy,show.node.label=TRUE)
edgelabels(1L:nrow(royer_phy$edge),
           edge=1L:nrow(royer_phy$edge),bg="white",cex=0.75)





#creation of PEM-------
#NB. manual testing using the first integrated fossil tree



royer_pgraph <- Phylo2DirectedGraph(intfossil_test)

steepness <- rep(0,attr(royer_pgraph,"ev")[1L])
evol_rate <- rep(1,attr(royer_pgraph,"ev")[1L])

royer_PEM <- PEM.build(royer_pgraph,
                       d="distance", sp="species",
                       a=steepness,psi=evol_rate)

royer_PEM_opt2 <- PEM.fitSimple(
  y = royer_train[,"log_lma"],
  x = royer_train[,"log_pet_leafarea"],
  w = royer_pgraph,
  d = "distance", sp="species",
  lower = 0, upper = 1)


royer_model_pet_lma <- lmforwardsequentialAICc(
  y = royer_train[,"log_lma"],
  x = royer_train[,"log_pet_leafarea",drop=FALSE],
  object = royer_PEM_opt2)

summary(royer_model_pet_lma)



royer_loc <- getGraphLocations(royer_phylo_lma$phy,
        targets = c("Piper_reticulatum","Piper_arboreum"))


pred_royer <- predict.PEM(
                object=royer_PEM_opt2,
                targets=royer_loc,
                lmobject=royer_model_pet_lma,
                newdata=royer_test,
                "prediction",0.95) %>% 
                as.data.frame()




#testing area for first phylo with fossil integration----
#NB. intfossil_test tree is too big and creates integer overflow---------------
test_pgraph <- Phylo2DirectedGraph(intfossil_test)

steepness <- rep(0,attr(test_pgraph,"ev")[1L])
evol_rate <- rep(1,attr(test_pgraph,"ev")[1L])

test_PEM <- PEM.build(test_pgraph,
                       d="distance", sp="species",
                       a=steepness,psi=evol_rate)

#####always make a regular lm of the log(LMA) and log((PW^2)/LA)





