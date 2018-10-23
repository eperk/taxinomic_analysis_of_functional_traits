#read from rds-------------------------------------------------------------------------
WSLA <- readRDS("~/Documents/BEIN data R/data/processed/02_WSLAPhen.rds")
WSLA_fam_count <- readRDS("~/Documents/BEIN data R/data/processed/02_family_count.rds")

#generic visualization and analysis-------------------------------------------------------------------------
histogram(~as.numeric(trait_value) | Phenology, data=subset(WSLA, as.numeric(trait_value)>1 & as.numeric(trait_value)<100))
histogram(~as.numeric(trait_value), data=subset(WSLA, WSLA$scrubbed_family == "Fabaceae"))
bwplot(Phenology~as.numeric(trait_value), data=subset(WSLA, as.numeric(trait_value)>1 & as.numeric(trait_value)<100))
favstats(~as.numeric(trait_value) | Phenology , data=subset(WSLA, as.numeric(trait_value)>1 & as.numeric(trait_value)<100))
TukeyHSD(as.numeric(trait_value)~Phenology, data=subset(WSLA, as.numeric(trait_value)>1 & as.numeric(trait_value)<100))

# takes all my WSLA data and then groups it by family and Phenology through pipes-------------------------------------------------------------------------
fam.phen.LMA <- WSLA_fam_count %>% 
  group_by(scrubbed_family, Phenology) %>%
   summarise(meanLMA = mean(log(as.numeric(LMA))))
fam.phen.SLA <- WSLA_fam_count %>% 
  group_by(scrubbed_family, Phenology) %>% 
   summarize(meanSLA = mean(log(as.numeric(trait_value))))

#clearingg na-------------------------------------------------------------------------
fam.phen.SLA <- na.omit(fam.phen.SLA)
fam.phen.LMA <- na.omit(fam.phen.LMA)

#statistical estimates
favstats(~as.numeric(meanLMA) | Phenology , data=fam.phen.LMA)
TukeyHSD(as.numeric(meanLMA)~Phenology, data=fam.phen.LMA)



#ggplots of phenology and lma by family: shows trends-------------------------------------------------------------------------
  ggplot(fam.phen.LMA,
         aes(x= Phenology, 
             y= meanLMA,
             color=scrubbed_family,
             group=scrubbed_family)) +
  geom_point() +
    geom_line()+
   theme(legend.position="none")+
  labs(y= "Mean Leaf Mass Area (LMA)")

  ggplot(fam.phen.SLA,
         aes(x= Phenology, 
             y= meanSLA,
             color=scrubbed_family,
             group=scrubbed_family)) +
    geom_point() +
    geom_line()+
    theme(legend.position="none")

            
            