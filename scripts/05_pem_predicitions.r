#PEM analysis-------------------------------------------------------------------------
#creation of all components required for PEM analysis-----------


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



royer_pgraph <- Phylo2DirectedGraph(royer_phy)

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





