require(tidyr)
require(dplyr)
require(mosaic)
require(lattice)
require(ggplot2)



histogram(~as.numeric(trait_value) | Phenology, data=subset(WSLA, as.numeric(trait_value)>1 & as.numeric(trait_value)<100))
#####Histogram of all trait values, subseted data over 1 and under 100 for m2/-1kg, split by phenology

histogram(~as.numeric(trait_value), data=subset(WSLA, WSLA$scrubbed_family == "Fabaceae"))
#####Histogram of all trait values, subseted data over 1 and under 100 for m2/-1kg, split by phenology

bwplot(Phenology~as.numeric(trait_value), data=subset(WSLA, as.numeric(trait_value)>1 & as.numeric(trait_value)<100))
#####bwplot of all traits by phenology between 1 and 100

favstats(~as.numeric(trait_value) | Phenology , data=subset(WSLA, as.numeric(trait_value)>1 & as.numeric(trait_value)<100))
##### favstats of above

TukeyHSD(as.numeric(trait_value)~Phenology, data=subset(WSLA, as.numeric(trait_value)>1 & as.numeric(trait_value)<100))
##### tukey-paired t-test of traits by phenology


fam.phen.LMA <- WSLA_fixed_lrgcount %>% 
  group_by(scrubbed_family, Phenology) %>%
   summarise(meanLMA = mean(log(as.numeric(LMA))))
##### takes all my WSLA data and then groups it by family and Phenology through pipes
fam.phen.SLA <- WSLA_fixed_lrgcount %>% 
  group_by(scrubbed_family, Phenology) %>% 
  summarize(meanSLA = mean(log(as.numeric(trait_value))))

fam.phen.SLA <- na.omit(fam.phen.SLA)
fam.phen.LMA <- na.omit(fam.phen.LMA)
fam.phen.SLA_raw <- fam.phen.SLA
##### Saves me an original version of the data
##### gets rid of N/A values

  ggplot(fam.phen.LMA,
         aes(x= Phenology, 
             y= meanLMA,
             color=scrubbed_family,
             group=scrubbed_family)) +
  geom_point() +
    geom_line()+
   theme(legend.position="none")
  #####Line and point plot, for phenology and LMA
  
  ggplot(fam.phen.SLA,
         aes(x= Phenology, 
             y= meanSLA,
             color=scrubbed_family,
             group=scrubbed_family)) +
    geom_point() +
    geom_line()+
    theme(legend.position="none")
  

  #####Line and point plot, for phenology and SLA
  ###### group= scrubbed family tells ggplot what geom_line should connect by

  
            
            