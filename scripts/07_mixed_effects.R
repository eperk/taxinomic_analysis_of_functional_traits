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

#creation of royer_top_5
royer_tax_top5 <- subset(royer_tax_count, as.numeric(count)>23)

#linear models that everything is based on
royer_lmfam <- lm(log_lma ~ log_pet_leafarea+scrubbed_family, data = royer_tax_na_omit)


royer_lmfam_top5 <- lm(log_lma ~ log_pet_leafarea+scrubbed_family, data = royer_tax_top5)
royer_lm_nofam_top5 <- lm(log_lma ~ log_pet_leafarea, data = royer_tax_top5)
royer_lmfam_top5_2 <- lm(log_lma ~ log_pet_leafarea*scrubbed_family, data = royer_tax_top5)


anova(royer_lmfam_top5_2, royer_lmfam_top5, royer_lm_nofam_top5)

#predicitons for fossils and family----
pred_lmfam <-as.data.frame(predict(royer_lmfam_top5, newdata=royer_fossil_dropped, interval="prediction"))
names(pred_lmfam)<-c("logLMAfam","lwrfam","uprfam")
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

#lm and prediciton for top 5 famlies
royer_lm_top5 <- lm(log_lma~log_pet_leafarea, data = royer_tax_top5)
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

ggplot(royer_lmfam_top5_bound, aes(log_pet_leafarea, log_lma))+
  geom_point() +
  geom_line(aes(y=lwr, group=interaction(scrubbed_family)), colour="red", linetype = "dashed")+
  geom_line(aes(y=fit, group=interaction(scrubbed_family)), color = "blue", linetype = "dashed")+
  geom_line(aes(y=upr, group=interaction(scrubbed_family)), color="red", linetype = "dashed")+
  geom_smooth(method=lm, se=TRUE)+
  theme(legend.position="none")


ggplot(royer_lmfam_top5_bound, )+
  aes(color=royer_lmfam_top5_bound$scrubbed_family)+
  geom_point() +
  geom_smooth(method=lm, se=TRUE)+
  labs(x="log(petiole^2/leafarea)")+
  labs(y="log(lma)")+
  labs(color = "Family")


##creation of fossil esitmates of top 5 fossil families

fossil_top5 <- royer_fossil_dropped%>% 
  group_by( scrubbed_family) %>% 
  summarize(count=n())
fossil_top5 <- left_join(fossil_top5, royer_fossil_dropped, by="scrubbed_family")
fossil_top5 <- subset(fossil_top5, as.numeric(count)>22)

#creation of extant off of the five top families

fossil_top5_tax <- as.data.frame(unique(fossil_top5$scrubbed_family))

extant_top5_fossil_tax <- subset(royer_tax_na_omit, (royer_tax_na_omit$scrubbed_family %in% fossil_top5_tax))


#prediciton interval nofamily top 5 fossils

extant_top5_fossil_lm <- lm(log_lma ~ log_pet_leafarea, data = extant_top5_fossil_tax)
pred_extant_top5_fossil_lm <-as.data.frame(predict(extant_top5_fossil_lm, newdata=fossil_top5, interval="prediction"))
fossil_top5_lm_predint <-  cbind(fossil_top5, pred_extant_top5_fossil_lm )

#prediction interval family top 5 fossils

extant_top5_fossil_lm_fam <- lm(log_lma ~ log_pet_leafarea+scrubbed_family, data = extant_top5_fossil_tax)
pred_extant_top5_fossil_lm_fam <-as.data.frame(predict(extant_top5_fossil_lm_fam, newdata=fossil_top5, interval="prediction"))
fossil_top5_lm_predint_fam <-  cbind(fossil_top5, pred_extant_top5_fossil_lm_fam )
colnames(fossil_top5_lm_predint_fam)[colnames(fossil_top5_lm_predint_fam)=="lwr"] <- "lwr_fam"
colnames(fossil_top5_lm_predint_fam)[colnames(fossil_top5_lm_predint_fam)=="fit"] <- "fit_fam"
colnames(fossil_top5_lm_predint_fam)[colnames(fossil_top5_lm_predint_fam)=="upr"] <- "upr_fam"


fossil_top5_lm_predint_fam_dropped <- subset(fossil_top5_lm_predint_fam, select = c(avg_LA, fit_fam, upr_fam, lwr_fam))
fossil_top5_all_pred <- left_join(fossil_top5_lm_predint_fam_dropped, fossil_top5_lm_predint, by = "avg_LA")


#loging? maybe----
fossil_top5_all_pred$lwr_fam <- log(fossil_top5_all_pred$lwr_fam)
fossil_top5_all_pred$upr_fam <- log(fossil_top5_all_pred$upr_fam)
fossil_top5_all_pred$lwr <- log(fossil_top5_all_pred$lwr)
fossil_top5_all_pred$upr <- log(fossil_top5_all_pred$upr)


anova(extant_top5_fossil_lm, extant_top5_fossil_lm_fam)

##finished ggplot, this is horrifying, but does show what we expected.
ggplot(fossil_top5_all_pred, aes(fit_fam, fit)) + 
  geom_point()+
  geom_abline()+
  geom_errorbarh(aes(xmin=lwr_fam, xmax=upr_fam))+
  geom_errorbar( aes(ymin=lwr, ymax=upr))+
  scale_x_continuous(name="family", limits=c(2, 7)) +
  scale_y_continuous(name="no family", limits=c(2, 7))
