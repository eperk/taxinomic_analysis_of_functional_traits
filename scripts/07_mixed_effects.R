#mixed effects modeling-------
#additional lm vs lmfam comparisons


royer_tax_full <- readRDS("~/Documents/BEIN data R/data/processed/07_lm4_royer")

royer_lme <- lmer(log_lma ~ log_pet_leafarea + (1|order/scrubbed_family/scrubbed_genus), data = royer_tax_full)
royer_lme_sum<- summary(royer_lme)

royer_pred <- as.data.frame(predict(royer_lme, newdata = all_fossil_royer_pred, allow.new.levels= TRUE))

#lm predictions and models using family as factor---------
royer_fossil_dropped <-   all_fossil_royer_pred %>% 
  filter(!grepl('Cercidiphyllaceae', scrubbed_family)) %>%
  filter(!grepl('Smilacaceae', scrubbed_family)) %>%
  filter(!grepl('Staphyleaceae', scrubbed_family))
  
#family lm----

#fixing royer_tax_full cause its messed up---------
royer_tax_fam_count <- royer_tax_full%>% 
  group_by( scrubbed_family) %>% 
  summarize(count=n()) %>% 
  filter(!grepl('Unknown', scrubbed_family))

royer_tax_count <- left_join(royer_tax_full, royer_tax_fam_count, by="scrubbed_family")

royer_tax_top5 <- subset(royer_tax_count, as.numeric(count)>23)


royer_lmfam <- lm(log_lma ~ log_pet_leafarea+scrubbed_family, data = royer_tax_na_omit)
royer_lmfam_top5 <- lm(log_lma ~ log_pet_leafarea+scrubbed_family, data = royer_tax_top5)

#predicitons for fossils and family----
pred_lmfam <-as.data.frame(predict(royer_lmfam, newdata=royer_fossil_dropped, interval="prediction"))
royer_lmfam_pred <-  cbind(royer_fossil_dropped, pred_lmfam )

#lm preiciton, just regular-----
royer_lm <- lm(log_lma ~ log_pet_leafarea, data = royer_tax_na_omit)
pred_lm <-as.data.frame(predict(royer_lm, newdata=all_fossil_royer_pred, interval="prediction"))
royer_lm_pred <-  cbind(all_fossil_royer_pred, pred_lm )

#lm prediction but for stuff that exists-------
pred_lm_extant <-as.data.frame(predict(royer_lm, interval="prediction"))
royer_lm_extant <-  cbind(royer_tax_full, pred_lm_extant)

#lm prediciton but for familys and stuff that exists------
pred_lmfam_extant <- as.data.frame(predict(royer_lmfam, interval="prediction"))
royer_lmfam_extant <-  cbind(royer_tax_na_omit, pred_lmfam_extant)

pred_lmfam_top5 <-  as.data.frame(predict(royer_lmfam_top5, interval="prediction"))
royer_lmfam_top5_bound <-  cbind(royer_tax_top5, pred_lmfam_top5)


#ggplot(royer_tax_full, aes(x = log_lma, y = log_pet_leafarea)) + 
 # geom_point(col = "red") +
  #stat_smooth(method = "lm", col = "red") 


#predicition interval for nofam----------

ggplot(royer_lm_pred, aes(log_pet_leafarea, log_LMA))+
  geom_point() +
  geom_line(aes(y=lwr), color = "red", linetype = "dashed")+
  geom_line(aes(y=upr), color = "red", linetype = "dashed")+
  geom_smooth(method=lm, se=TRUE)

#preiction interval for fam-----
ggplot(royer_lmfam_pred, aes(log_pet_leafarea, log_LMA))+
  geom_point() +
  geom_line(aes(y=lwr), color = "red", linetype = "dashed")+
  geom_line(aes(y=upr), color = "red", linetype = "dashed")+
  geom_smooth(method=lm, se=TRUE)

###prediction interval for extant--------
ggplot(royer_lm_extant, aes(log_pet_leafarea, log_lma))+
  geom_point() +
  geom_line(aes(y=lwr), color = "red", linetype = "dashed")+
  geom_line(aes(y=fit), color = "blue", linetype = "dashed")+
  geom_line(aes(y=upr), color = "red", linetype = "dashed")+
  geom_smooth(method=lm, se=TRUE)


#pred interval for extant+family-----
ggplot(royer_lmfam_extant, aes(log_pet_leafarea, log_lma, group=interaction(scrubbed_family)))+
  geom_point() +
  geom_line(aes(y=lwr), color = "red", linetype = "dashed")+
  geom_line(aes(y=fit), color = "blue", linetype = "dashed")+
  geom_line(aes(y=upr), color = "red", linetype = "dashed")+
  geom_smooth(method=lm, se=TRUE)+
  theme(legend.position="none")

#make a color pallete--------------
myColors <- brewer.pal(5,"Set1")
names(myColors) <- levels(royer_lmfam_top5_bound$scrubbed_family)
colScale <- scale_colour_manual(name = "scrubbed_family",values = myColors)

ggplot(royer_lmfam_top5_bound, aes(log_pet_leafarea, log_lma))+
  geom_point() +
  geom_line(aes(y=lwr, group=interaction(scrubbed_family)), colour="red", linetype = "dashed")+
  geom_line(aes(y=fit, group=interaction(scrubbed_family)), color = "blue", linetype = "dashed")+
  geom_line(aes(y=upr, group=interaction(scrubbed_family)), color="red", linetype = "dashed")+
  geom_smooth(method=lm, se=TRUE)+
  theme(legend.position="none")

