#load tree and files-------------------------------------------------------------------------------------------------
tree_WSLA_species <- read.tree("~/Documents/BEIN data R/data/processed/03_WSLA_species.tre")
final_WSLA_DF <- readRDS("~/Documents/BEIN data R/data/processed/03_finalWSLADF.rds")

#sla anc analysis-------------------------------------------------------------------------------------------------
SLA_only <- final_WSLA_DF [-c(3:6)]
SLA_only$SLA<- as.numeric(SLA_only$SLA)
SLA_only_means <- aggregate(SLA_only[,-1], by=list(SLA_only$binomial), mean)
colnames(SLA_only_means)[colnames(SLA_only_means)=="Group.1"] <- "binomial"
colnames(SLA_only_means)[colnames(SLA_only_means)=="x"] <- "SLA"
SLA_mean_tree_binomials <- subset(SLA_only_means, SLA_only_means$binomial %in% WSLA_tree_tips)
SLA_labels <- SLA_mean_tree_binomials[,-1]
SLA_labels <- as.data.frame(SLA_labels)
rownames(SLA_labels) <- SLA_mean_tree_binomials[,1]
SLA <-as.matrix(SLA_labels)[,1]
SLA_fastAnc<-fastAnc(tree_WSLA_species,SLA,vars=TRUE,CI=TRUE)
SLA_mapping<-contMap(tree_WSLA_species,SLA,plot=FALSE)
plot(SLA_mapping,type="fan",legend=0.7*max(nodeHeights(tree_WSLA_species)), fsize=c(0.7,0.9))


#rosidae fastanc analysis-------------------------------------------------------------------------------------------------
WSLA_tree_nodes <- tree_WSLA_species$node.label
Rosidae <- as.vector("Rosidae")
tree_rosidae <- extract.clade(tree_WSLA_species, "Rosidae")
rosidae_tips<- tree_rosidae$tip.label
rosidae_SLA_mean_binomials <- subset(SLA_only_means, SLA_only_means$binomial %in% rosidae_tips)
rosidae_labels <- rosidae_SLA_mean_binomials[,-1]
rosidae_labels <- as.data.frame(rosidae_labels)
rownames(rosidae_labels) <- rosidae_SLA_mean_binomials[,1]
rosidae_matrix <-as.matrix(rosidae_labels)[,1]
rosidae_fastAnc<-fastAnc(tree_rosidae,rosidae_matrix,vars=TRUE,CI=TRUE)
rosidae_mapping<-contMap(tree_rosidae,rosidae_matrix,plot=FALSE)
plot(rosidae_mapping,legend=0.7*max(nodeHeights(tree_rosidae)), fsize=c(0.7,0.9), ftype="off")
plot(tree_rosidae, type="fan", show.tip.label = FALSE)


#simmap rosidae-------------------------------------------------------------------------------------------------
species_state <- final_WSLA_DF [-c(2,3,5,6)]
deduped_states <- unique( species_state[ , 1:2 ] )
deomited_states <- na.omit(deduped_states)
species_state_rosidae <- subset(deomited_states, (deomited_states$binomial %in% rosidae_tips))
species_pheno_rosidae<-species_state_rosidae$Phenology

#simmap tree creation-------------------------------------------------------------------------------------------------
not_in_rosidae <- subset(rosidae_tips, !(rosidae_tips %in% species_state_rosidae$binomial))
not_in_rosidae_df <-as.data.frame(not_in_rosidae)
tree_simmap_rosidae <- drop.tip(tree_rosidae, not_in_rosidae)
names(species_pheno_rosidae) <- species_state_rosidae$binomial

#rosidae simmap-------------------------------------------------------------------------------------------------
rosidae_simmap <- make.simmap(tree_simmap_rosidae, species_pheno_rosidae, nsim=1)
plotSimmap(rosidae_simmap, ftype="off")

#binary vs continous trait phylogenetic logistic regressions-------------------------------------------------------------------------------------------------
SLA_pheno_binomal <- final_WSLA_DF [-c(5:6)]
SLA_pheno_binomal$SLA<- as.numeric(SLA_pheno_binomal$SLA)
SLA_LMA_pheno_means <- aggregate(SLA_pheno_binomal[,-1], by=list(SLA_pheno_binomal$binomial), mean)
colnames(SLA_LMA_pheno_means)[colnames(SLA_LMA_pheno_means)=="Group.1"] <- "binomial"
SLA_LMA_pheno_means <- left_join(SLA_LMA_pheno_means, Zanne, by= "binomial")
SLA_LMA_pheno_means <- SLA_LMA_pheno_means [-c(4)]
colnames(SLA_LMA_pheno_means)[colnames(SLA_LMA_pheno_means)=="Phenology.y"] <- "Phenology"
SLA_LMA_pheno_means <- na.omit(SLA_LMA_pheno_means)

#rosidae binary vs continous trait phylogenetic logistic regressions-------------------------------------------------------------------------------------------------
rosidae_SLA_LMA_pheno_means <- subset(SLA_LMA_pheno_means, SLA_LMA_pheno_means$binomial %in% rosidae_tips)
rosidae_labels_SLA_LMA <- rosidae_SLA_LMA_pheno_means[,1]
rosidae_labels_SLA_LMA <- as.data.frame(rosidae_labels_SLA_LMA)
rownames(rosidae_SLA_LMA_pheno_means) <- rosidae_labels_SLA_LMA[,1]

# rosidae Dropping D_EV so binary trait is achieved-------------------------------------------------------------------------------------------------
rosidae_SLA_LMA_pheno_means <- rosidae_SLA_LMA_pheno_means[!rosidae_SLA_LMA_pheno_means$Phenology == "D_EV", ]

#rosiade tree creation for phylo logistic regression-------------------------------------------------------------------------------------------------
not_in_rosidae_na_D_EV <- subset(rosidae_tips, !(rosidae_tips %in% rosidae_SLA_LMA_pheno_means$binomial))
tree_simmap_rosidae_na_D_EV <- drop.tip(tree_rosidae, not_in_rosidae_na_D_EV)
tree_simmap_rosidae_na_D_EV_tips <- tree_simmap_rosidae_na_D_EV$tip.label
tree_simmap_rosidae_na_D_EV_tips <- as.data.frame(tree_simmap_rosidae_na_D_EV_tips)
rosidae_SLA_LMA_pheno_means <- rosidae_SLA_LMA_pheno_means[ order(match(rosidae_SLA_LMA_pheno_means$binomial, 
              tree_simmap_rosidae_na_D_EV_tips$tree_simmap_rosidae_na_D_EV_tips)), ]
rosidae_SLA_LMA_pheno_means$Phenology <- as.character(rosidae_SLA_LMA_pheno_means$Phenology)

#rosidae Decidous=1 Evergreen=0 creation of binary factor-------------------------------------------------------------------------------------------------
rosidae_SLA_LMA_pheno_means$Phenology[rosidae_SLA_LMA_pheno_means$Phenology == "D"] <- "1"
rosidae_SLA_LMA_pheno_means$Phenology[rosidae_SLA_LMA_pheno_means$Phenology == "EV"] <- "0"

#rosidae phylogenetic linear regression-------------------------------------------------------------------------------------------------
phylin <- phylolm(SLA~Phenology, data = rosidae_SLA_LMA_pheno_means, phy = tree_simmap_rosidae_na_D_EV)
summary(phylin)
plot.phylolm(phylin)


#rosidae phylogenetic logistic regression-------------------------------------------------------------------------------------------------
phylog <- phyloglm(as.numeric(Phenology)~SLA, data = rosidae_SLA_LMA_pheno_means, 
                   phy = tree_simmap_rosidae_na_D_EV, method = c("logistic_IG10"))
summary(phylog)
plot.phyloglm(phylog)

#rosidae ggplot of phylog-------------------------------------------------------------------------------------------------
coefs <- coef(phylog)
x_plot <- seq(-50, 50, by = 0.1)
y_plot <- plogis(coefs[1] + coefs[2] * x_plot)
plot_data <- data.frame(x_plot, y_plot)

ggplot(plot_data) + 
  geom_line(aes(x_plot, 
                y_plot), 
            col = "red") + 
  xlab("x") + 
  ylab("p(y | x)")

#binary vs continous trait phylogenetic logistic regressions for all species-------------------------------------------------------------------------------------------------
SLA_pheno_binomal <- final_WSLA_DF [-c(5:6)]
SLA_pheno_binomal$SLA<- as.numeric(SLA_pheno_binomal$SLA)
SLA_LMA_pheno_means <- aggregate(SLA_pheno_binomal[,-1], by=list(SLA_pheno_binomal$binomial), mean)
colnames(SLA_LMA_pheno_means)[colnames(SLA_LMA_pheno_means)=="Group.1"] <- "binomial"
SLA_LMA_pheno_means <- left_join(SLA_LMA_pheno_means, Zanne, by= "binomial")
SLA_LMA_pheno_means <- SLA_LMA_pheno_means [-c(4)]
colnames(SLA_LMA_pheno_means)[colnames(SLA_LMA_pheno_means)=="Phenology.y"] <- "Phenology"
SLA_LMA_pheno_means <- na.omit(SLA_LMA_pheno_means)
all_SLA_LMA_pheno_means <- subset(SLA_LMA_pheno_means, SLA_LMA_pheno_means$binomial %in% WSLA.tree.tips)
all_SLA_LMA_pheno_means_labels <- all_SLA_LMA_pheno_means[,1]
all_SLA_LMA_pheno_means_labels <- as.data.frame(all_SLA_LMA_pheno_means_labels)
rownames(all_SLA_LMA_pheno_means) <- all_SLA_LMA_pheno_means_labels[,1]

#Dropping D_EV so binary trait is achieved-------------------------------------------------------------------------------------------------
indx<-which(!all_SLA_LMA_pheno_means$Phenology == "D_EV")
all_SLA_LMA_pheno_means <- all_SLA_LMA_pheno_means[indx, ]

#allspecies phylo logsitic tree creation-------------------------------------------------------------------------------------------------
not_in_all_na_D_EV <- subset(WSLA.tree.tips, !(WSLA.tree.tips %in% all_SLA_LMA_pheno_means$binomial))
tree_WSLA_na_WSLA <- drop.tip(tree_WSLA_species, not_in_all_na_D_EV)
tree_WSLA_na_WSLA_tips <- tree_WSLA_na_WSLA$tip.label
tree_WSLA_na_WSLA_tips <- as.data.frame(tree_WSLA_na_WSLA_tips)


#matching to tree-------------------------------------------------------------------------------------------------
all_SLA_LMA_pheno_means <- all_SLA_LMA_pheno_means[ order(match(all_SLA_LMA_pheno_means$binomial, 
                                                tree_WSLA_na_WSLA_tips$tree_WSLA_na_WSLA_tips)), ]

#Decidous=1 Evergreen=0 creation of binary factor-------------------------------------------------------------------------------------------------
all_SLA_LMA_pheno_means$Phenology[all_SLA_LMA_pheno_means$Phenology == "D"] <- "0"
all_SLA_LMA_pheno_means$Phenology[all_SLA_LMA_pheno_means$Phenology == "EV"] <- "1"

#phylogenetic linear regression-------------------------------------------------------------------------------------------------
phylin_all <- phylolm(SLA~Phenology, data = all_SLA_LMA_pheno_means, phy = tree_WSLA_na_WSLA)
summary(phylin_all)
plot.phylolm(phylin_all)


##phylogenetic logistic regression-------------------------------------------------------------------------------------------------
phylog_all <- phyloglm(as.numeric(Phenology)~SLA, data = all_SLA_LMA_pheno_means, 
                       phy = tree_WSLA_na_WSLA, method = c("logistic_IG10"))
summary(phylog_all)
plot.phyloglm(phylog_all)


#phylogenetic logistic regression, log-------------------------------------------------------------------------------------------------
phylog_all_log <- phyloglm(as.numeric(Phenology)~log(SLA), data = all_SLA_LMA_pheno_means, phy = tree_WSLA_na_WSLA, method = c("logistic_IG10"))
summary(phylog_all_log)
plot.phyloglm(phylog_all_log)



#ggplot of phylog.all-------------------------------------------------------------------------------------------------
coefs_all <- coef(phylog_all)
x_plot_all <- seq(0, 55, by = 0.1)
y_plot_all <- plogis(coefs_all[1] + coefs_all[2] * x_plot_all)
plot_data_all <- data.frame(x_plot_all, y_plot_all)

ggplot(plot_data_all) + 
  geom_line(aes(x_plot_all, 
                y_plot_all), 
            col = "red") + 
  xlab("SLA") + 
  ylab("p(y | x)") +
  geom_jitter(data = all_SLA_LMA_pheno_means,
              aes(
                x = SLA,
                y = as.numeric(Phenology)
              ))

#log plotting-------------------------------------------------------------------------------------------------
coef_all_log <- coef(phylog_all_log)
x_plot_all_log <- seq(0, 4, by = 0.1)
y_plot_all_log <- plogis(coef_all_log[1] + coef_all_log[2] * x_plot_all_log)
plot_data_all_log <- data.frame(x_plot_all_log, y_plot_all_log)

ggplot(plot_data_all_log) + 
  geom_line(aes(x_plot_all_log, 
                y_plot_all_log), 
            col = "red") + 
  xlab("log SLA") + 
  ylab("p(y | x)") +
  geom_point(aes(
    
  ))



#ABANDON ALL HOPE YE WHO ENTER THIS PART OF THE CODE, TESTING BELOW-------------------------------------------------------------------------------------------------

phylog_all_ape <- compar.gee(as.numeric(Phenology)~SLA, data = all_SLA_LMA_pheno_means, 
                             phy = tree_WSLA_na_WSLA, family = binomial())
phylog_predict <- predict(phylog_all_ape)


########maybe this will work with smaller tree?
#myrtales_SLA_mean_binomials_pheno <- left_join(myrtales.SLA.mean.binomials, Zanne, by = "binomial")
#myrtales_SLA_mean_binomials_pheno <- na.omit(myrtales_SLA_mean_binomials_pheno)

#myrtales_SLA_mean_binomials_pheno

#myrtales_SLA_mean_binomials_pheno$Phenology[myrtales_SLA_mean_binomials_pheno$Phenology == "D"] <- "0"
#myrtales_SLA_mean_binomials_pheno$Phenology[myrtales_SLA_mean_binomials_pheno$Phenology == "EV"] <- "1"

#not_in_myrtales_new <- subset(myrtales.tips, !(myrtales.tips %in% myrtales_SLA_mean_binomials_pheno$binomial))
#tree.phyloglm.myrtales <- drop.tip(tree.myrtales, not.in.myrtales.new)

#tree.phyloglm.myrtales.tips <- tree.phyloglm.myrtales$tip.label

#myrtales.SLA.labels <- myrtales_SLA_mean_binomials_pheno[,1]
#myrtales.SLA.labels <- as.data.frame(myrtales.SLA.labels)
#rownames(myrtales_SLA_mean_binomials_pheno) <- myrtales.SLA.labels[,1]

#myrtales_SLA_mean_binomials_pheno <- myrtales_SLA_mean_binomials_pheno[ order(match(myrtales_SLA_mean_binomials_pheno$binomial, 
    #                                                                                tree.phyloglm.myrtales$tree.phyloglm.myrtales.tips)), ]

#phylog.myrtales<- compar.gee(as.numeric(Phenology)~SLA, data = myrtales_SLA_mean_binomials_pheno, 
 #                           phy = tree.phyloglm.myrtales, family = binomial())




################ residual trait analysis
phylog.all.resid <- residuals(phylog.all)
phylog.all.resid <- as.data.frame(phylog.all.resid)

phylog.resids <-as.matrix(phylog.all.resid)[,1]
resids.fastAnc<-fastAnc(tree_WSLA_na_WSLA,phylog.resids,vars=TRUE,CI=TRUE)
resids.mapping<-contMap(tree_WSLA_na_WSLA,phylog.resids,plot=FALSE)
plot(resids.mapping,type="fan",legend=0.7*max(nodeHeights(tree_WSLA_na_WSLA)), fsize=c(0.7,0.9))

################regular logistic regressions
SLA.LMA.glm.predict.families <- left_join(all_SLA_LMA_pheno_means, family_binomial, by = "binomial")
phylo.glm <- glm(as.numeric(Phenology)~SLA, data= all_SLA_LMA_pheno_means, family = binomial())
summary(phylo.glm)

predict.glm(SLA.LMA.glm.predict.families)

###################### family level predictions?




##maybe phylogenetic signal will be useful?
###phylosig(tree = tree_WSLA_na_WSLA, )

####################################


######## from what i 

#####run regular logistic regression, compare predictions
###### go through Garland and Ives again, look for ways to compare the regular and the phylogenetic
###### how can we compare phylogentic to nonphylogenetic
###### how does it handle making prediciotns? how does it work with phylo vs nonphylo
###### lm predicitons with phylogenetics  how to predcit new data with phylogenetic logistic regression
##### compare predicitons that are generated from the regular regression and the phylogenetic regression
##### look into phyloglm package and how it treats the phylogenetic information, what bearing does the phylogenetic info have on predicitons?
##### is the phylogenetic structure in the resisduals?
########is phylogenetic model more informative? look at literature
############ use coefficents to try and predict the the phyloglm model???????




