#mixed effects modeling-------

royer_tax_full <- readRDS("~/Documents/BEIN data R/data/processed/07_lm4_royer")

royer_lme <- lmer(log_lma ~ log_pet_leafarea + (1|order/scrubbed_family/scrubbed_genus), data = royer_tax_full)
royer_lme_sum<- summary(royer_lme)

royer_pred <- preidct(royer_lme, newdata = )