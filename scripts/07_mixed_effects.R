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
royer_tax_na_omit <- na.omit(royer_tax_full)
royer_lmfam <- lm(log_lma ~ log_pet_leafarea+scrubbed_family, data = royer_tax_na_omit)

#predicitons for fossils and family----
pred_lmfam <-as.data.frame(predict(royer_lmfam, newdata=royer_fossil_dropped, interval="prediction"))
royer_lmfam_pred <-  cbind(royer_fossil_dropped, pred_lmfam )

#lm preiciton, just regular-----
royer_lm <- lm(log_lma ~ log_pet_leafarea, data = royer_tax_full)
pred_lm <-as.data.frame(predict(royer_lm, newdata=all_fossil_royer_pred, interval="prediction"))
royer_lm_pred <-  cbind(all_fossil_royer_pred, pred_lm )

#lm prediction but for stuff that exists-------
pred_lm_extant <-as.data.frame(predict(royer_lm, interval="prediction"))
royer_lm_extant <-  cbind(royer_tax_full, pred_lm_extant)

#lm prediciton but for familys and stuff that exists------
pred_lmfam_extant <- as.data.frame(predict(royer_lmfam, interval="prediction"))
royer_lmfam_extant <-  cbind(royer_tax_na_omit, pred_lmfam_extant)


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

###prediction interval for extant
ggplot(royer_lm_extant, aes(log_pet_leafarea, log_lma))+
  geom_point() +
  geom_line(aes(y=lwr), color = "red", linetype = "dashed")+
  geom_line(aes(y=fit), color = "blue", linetype = "dashed")+
  geom_line(aes(y=upr), color = "red", linetype = "dashed")+
  geom_smooth(method=lm, se=TRUE)


#pred interval for extant+family

ggplot(royer_lmfam_extant, aes(log_pet_leafarea, log_lma))+
  geom_point() +
  geom_line(aes(y=lwr), color = "red", linetype = "dashed")+
  geom_line(aes(y=fit), color = "blue", linetype = "dashed")+
  geom_line(aes(y=upr), color = "red", linetype = "dashed")+
  geom_smooth(method=lm, se=TRUE)
