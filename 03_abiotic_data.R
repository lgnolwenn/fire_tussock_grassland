# ANALYSES OF ABIOTIC DATA FROM DEEP STREAM
## Code for cleaning and analysing the abiotic parameters from Deep Stream grassland.
## Note: paths need to be changed to match your device.

# LOAD PACKAGES
library(tidyverse)
library(vegan)
library(ggbiplot)
library(lme4)
library(emmeans)
library(ggeffects)

# LOAD DATA
raw_abiotic_data<-read.csv("C:/Users/nolwe/Documents/M2/internship/project/abiotic_soil.csv")

#CLEANING

abiotic_data<-raw_abiotic_data

##add treatment variable: unburnt (UB), control (CO), spring (SP) and summer (SU)
#CO for control not burnt in 2001 but in 2019, SU for summer burnt in 2001 and burnt in 2019, 
#SP for spring burnt in 2001 and burnt in 2019, UB for plots that never burnt
abiotic_data$treatment<-substr(abiotic_data$plot_id, 1, 2)

##remove E plots: exclosure subplots are not interesting for our analyses
abiotic_data<-subset(abiotic_data, !grepl("_E_",abiotic_data$subplot_id))

## remove useless variables
variables=colnames(abiotic_data)
### Remove duplicated measures of minerals in different units. 
### Keep the % of Cation Exchange Capacity unit.
maf=grepl("MAF", variables) # variables in MAF units
#### Variables in me/100mg unit
cat=grepl(paste(c("K","Ca","Mg","Na"), collapse='|'), variables) 
bs=grepl("BS", variables)
for(i in 1:length(variables)){if (bs[i]){cat[i]<-FALSE}}
### Remove variables used for calculating others
calc=grepl(paste(c("CEC", "saturation", "volume.weight"), collapse='|'), variables)
### Remove some C and N measures
drop_cn=grepl(paste(c("mineralisable.N", "N.ratio", "total.C", "total.N"), collapse='|'), variables)
### Drop all the variables to remove
abiotic_data_2<-abiotic_data[,!maf&!cat&!calc&!drop_cn]
variables_2=colnames(abiotic_data_2)[-c(1:3)]

# VISUALISATION
abiotic.pca <- princomp(abiotic_data_2[,-c(1:4,14)], cor=TRUE)
ggbiplot(abiotic.pca) +
  geom_point(aes(colour = abiotic_data_2$plot)) +
  labs(colour = "Site", x = "PC1", y = "PC2") +
  theme_bw()

# STATISTICS ON INTERESTING VARIABLES

## Select useful variables
abiotic_data_2$fire=substr(abiotic_data_2$plot_id, 1, 2)
abiotic_data_plot=abiotic_data_2 %>% select(-subplot_id, -time, -treatment, -sample)
abiotic_data_plot=unique(abiotic_data_plot)
### Get mean and standard deviation for each variable
mean_sd_abiotic<-abiotic_data_plot %>% group_by(fire)%>% select(-plot_id) %>% summarise_all(., list(name = mean, sd= sd))
abiotic_data_2$fire=as.factor(abiotic_data_2$fire)

## pH

lm.ph<-lm(pH~fire, data=abiotic_data_2) #linear model
em.ph<-emmeans(lm.ph, "fire") # contrasts between fire histories
pairs(em.ph)
pred_ph <- ggpredict(lm.ph, terms="fire") # predictions
names(pred_ph) <- c("fire", "predicted", "std.error", "conf.low", "conf.high")

plot_ph <- 
  ggplot(pred_ph, aes(x = fire, y = predicted, colour=fire, group=fire)) + #
  geom_point(position=position_dodge(0.6), size=3) +
  geom_point(data = abiotic_data_2, aes(x = fire, y = pH , fill = fire), position = position_dodge(width = 0.6), alpha = 0.5, size=2) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, width = 0.3), position=position_dodge(0.6), linewidth=1) +
  scale_color_manual(values=c("#FF9900", "#00b0f0", "#6633FF", "#04c728"), labels=c("control", "spring", "summer", "unburnt")) +
  scale_x_discrete(limits=c("UB","CO","SP", "SU"), breaks=c("UB","CO","SP", "SU"), labels=c("unburnt", "control", "spring", "summer"))+
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(colour="black", face="italic", size=12)) + 
  theme(axis.text.y = element_text(colour="black", size=12)) +
  ylab("pH") +
  xlab(" Fire history") +
  theme(aspect.ratio=1) +
  theme(legend.position="none") 

## olsen P

lm.p<-lm(Olsen.P~fire, data=abiotic_data_2) #linear model
em.p<-emmeans(lm.p, "fire") # contrasts between fire histories
pairs(em.p)
pred_p <- ggpredict(lm.p, terms="fire") # predictions
names(pred_p) <- c("fire", "predicted", "std.error", "conf.low", "conf.high")

plot_p <- 
  ggplot(pred_p, aes(x = fire, y = predicted, colour=fire, group=fire)) + #
  geom_point(position=position_dodge(0.6), size=3) +
  geom_point(data = abiotic_data_2, aes(x = fire, y = Olsen.P , fill = fire), position = position_dodge(width = 0.6), alpha = 0.5, size=2) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, width = 0.3), position=position_dodge(0.6), linewidth=1) +
  scale_color_manual(values=c("#FF9900", "#00b0f0", "#6633FF", "#04c728"), labels=c("control", "spring", "summer", "unburnt")) +
  scale_x_discrete(limits=c("UB","CO","SP", "SU"), breaks=c("UB","CO","SP", "SU"), labels=c("unburnt", "control", "spring", "summer"))+
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(colour="black", face="italic", size=12)) + 
  theme(axis.text.y = element_text(colour="black", size=12)) +
  labs(x= "Fire history", y=expression(Olsen~Phosphorus~"("~italic(mg. ~ L^-1)~")")) +
  theme(aspect.ratio=1) +
  theme(legend.position="none") 

## Available N

lm.n<-lm(available.N~fire, data=abiotic_data_2) #linear model
em.n<-emmeans(lm.n, "fire")# contrasts between fire histories
pairs(em.n)
pred_n <- ggpredict(lm.n, terms="fire") #predictions
names(pred_n) <- c("fire", "predicted", "std.error", "conf.low", "conf.high")

plot_n <- 
  ggplot(pred_n, aes(x = fire, y = predicted, colour=fire, group=fire)) + #
  geom_point(position=position_dodge(0.6), size=3) +
  geom_point(data = abiotic_data_2, aes(x = fire, y = available.N , fill = fire), position = position_dodge(width = 0.6), alpha = 0.5, size=2) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, width = 0.3), position=position_dodge(0.6), linewidth=1) +
  scale_color_manual(values=c("#FF9900", "#00b0f0", "#6633FF", "#04c728"), labels=c("control", "spring", "summer", "unburnt")) +
  scale_x_discrete(limits=c("UB","CO","SP", "SU"), breaks=c("UB","CO","SP", "SU"), labels=c("unburnt", "control", "spring", "summer"))+
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(colour="black", face="italic", size=12)) + 
  theme(axis.text.y = element_text(colour="black", size=12)) +
  labs(x= "Fire history", y=expression(Available~Nitrogen~"("~italic(kg.~ ha^-1)~")")) +
  theme(aspect.ratio=1) +
  theme(legend.position="none") 

## Organic.matter
lm.om<-lm(organic.matter~fire, data=abiotic_data_2) #linear model
em.om<-emmeans(lm.om, "fire") # contrasts between fire histories
pairs(em.om)
pred_om <- ggpredict(lm.om, terms="fire") #predictions
names(pred_om) <- c("fire", "predicted", "std.error", "conf.low", "conf.high")


plot_om <- 
  ggplot(pred_om, aes(x = fire, y = predicted, colour=fire, group=fire)) + #
  geom_point(position=position_dodge(0.6), size=3) +
  geom_point(data = abiotic_data_2, aes(x = fire, y = organic.matter , fill = fire), position = position_dodge(width = 0.6), alpha = 0.5, size=2) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, width = 0.3), position=position_dodge(0.6), linewidth=1) +
  scale_color_manual(values=c("#FF9900", "#00b0f0", "#6633FF", "#04c728"), labels=c("control", "spring", "summer", "unburnt")) +
  scale_x_discrete(limits=c("UB","CO","SP", "SU"), breaks=c("UB","CO","SP", "SU"), labels=c("unburnt", "control", "spring", "summer"))+
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(colour="black", face="italic", size=12)) + 
  theme(axis.text.y = element_text(colour="black", size=12)) +
  labs(x= "Fire history", y=expression(Organic~matter~"("~italic("%")~")")) +
  theme(aspect.ratio=1) +
  theme(legend.position="none") 

## Supplementary plots
### Carbon:Nitrogen
lm.cn<-lm(C.N~fire, data=abiotic_data_2) #linear model
em.cn<-emmeans(lm.cn, "fire") #contrasts between fire histories
pairs(em.cn)
pred_cn <- ggpredict(lm.cn, terms="fire") #predictions
names(pred_cn) <- c("fire", "predicted", "std.error", "conf.low", "conf.high")

plot_cn <- 
  ggplot(pred_cn, aes(x = fire, y = predicted, colour=fire, group=fire)) + #
  geom_point(position=position_dodge(0.6), size=3) +
  geom_point(data = abiotic_data_2, aes(x = fire, y = C.N , fill = fire), position = position_dodge(width = 0.6), alpha = 0.5, size=2) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, width = 0.3), position=position_dodge(0.6), linewidth=1) +
  scale_color_manual(values=c("#FF9900", "#00b0f0", "#6633FF", "#04c728"), labels=c("control", "spring", "summer", "unburnt")) +
  scale_x_discrete(limits=c("UB","CO","SP", "SU"), breaks=c("UB","CO","SP", "SU"), labels=c("unburnt", "control", "spring", "summer"))+
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(colour="black", face="italic", size=12)) + 
  theme(axis.text.y = element_text(colour="black", size=12)) +
  labs(x= "Fire history", y="Carbon:Nitrogen ratio") +
  theme(aspect.ratio=1) +
  theme(legend.position="none") 

### Potassium
lm.k<-lm(K.BS~fire, data=abiotic_data_2) #linear model
em.k<-emmeans(lm.k, "fire") #contrasts between fire histories
pairs(em.k)
pred_k <- ggpredict(lm.k, terms="fire") #predictions
names(pred_k) <- c("fire", "predicted", "std.error", "conf.low", "conf.high")

plot_k<-
  ggplot(pred_k, aes(x = fire, y = predicted, colour=fire, group=fire)) + #
  geom_point(position=position_dodge(0.6), size=3) +
  geom_point(data = abiotic_data_2, aes(x = fire, y = K.BS , fill = fire), position = position_dodge(width = 0.6), alpha = 0.5, size=2) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, width = 0.3), position=position_dodge(0.6), linewidth=1) +
  scale_color_manual(values=c("#FF9900", "#00b0f0", "#6633FF", "#04c728"), labels=c("control", "spring", "summer", "unburnt")) +
  scale_x_discrete(limits=c("UB","CO","SP", "SU"), breaks=c("UB","CO","SP", "SU"), labels=c("unburnt", "control", "spring", "summer"))+
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(colour="black", face="italic", size=12)) + 
  theme(axis.text.y = element_text(colour="black", size=12)) +
  labs(x= "Fire history", y=expression(Potassium~"("~italic("% of CEC")~")")) +
  theme(aspect.ratio=1) +
  theme(legend.position="none") 

### Calcium
lm.ca<-lm(Ca.BS~fire, data=abiotic_data_2) #linear model
em.ca<-emmeans(lm.ca, "fire") #contrasts between fire histories
pairs(em.ca)
pred_ca <- ggpredict(lm.ca, terms="fire") #predictions
names(pred_ca) <- c("fire", "predicted", "std.error", "conf.low", "conf.high")

plot_ca<-
  ggplot(pred_ca, aes(x = fire, y = predicted, colour=fire, group=fire)) + #
  geom_point(position=position_dodge(0.6), size=3) +
  geom_point(data = abiotic_data_2, aes(x = fire, y = CA.BS , fill = fire), position = position_dodge(width = 0.6), alpha = 0.5, size=2) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, width = 0.3), position=position_dodge(0.6), linewidth=1) +
  scale_color_manual(values=c("#FF9900", "#00b0f0", "#6633FF", "#04c728"), labels=c("control", "spring", "summer", "unburnt")) +
  scale_x_discrete(limits=c("UB","CO","SP", "SU"), breaks=c("UB","CO","SP", "SU"), labels=c("unburnt", "control", "spring", "summer"))+
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(colour="black", face="italic", size=12)) + 
  theme(axis.text.y = element_text(colour="black", size=12)) +
  labs(x= "Fire history", y=expression(Calcium~"("~italic("% of CEC")~")")) +
  theme(aspect.ratio=1) +
  theme(legend.position="none") 

### Magnesium
lm.mg<-lm(Mg.BS~fire, data=abiotic_data_2) #linear model
em.mg<-emmeans(lm.mg, "fire") #contrasts between fire histories
pairs(em.mg)
pred_mg <- ggpredict(lm.mg, terms="fire") #predictions
names(pred_mg) <- c("fire", "predicted", "std.error", "conf.low", "conf.high")

plot_mg<-
  ggplot(pred_mg, aes(x = fire, y = predicted, colour=fire, group=fire)) + #
  geom_point(position=position_dodge(0.6), size=3) +
  geom_point(data = abiotic_data_2, aes(x = fire, y = Mg.BS , fill = fire), position = position_dodge(width = 0.6), alpha = 0.5, size=2) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, width = 0.3), position=position_dodge(0.6), linewidth=1) +
  scale_color_manual(values=c("#FF9900", "#00b0f0", "#6633FF", "#04c728"), labels=c("control", "spring", "summer", "unburnt")) +
  scale_x_discrete(limits=c("UB","CO","SP", "SU"), breaks=c("UB","CO","SP", "SU"), labels=c("unburnt", "control", "spring", "summer"))+
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(colour="black", face="italic", size=12)) + 
  theme(axis.text.y = element_text(colour="black", size=12)) +
  labs(x= "Fire history", y=expression(Magnesium~"("~italic("% of CEC")~")")) +
  theme(aspect.ratio=1) +
  theme(legend.position="none") 

### Sodium
lm.na<-lm(Na.BS~fire, data=abiotic_data_2) #linear model
em.na<-emmeans(lm.na, "fire") #contrasts between fire histories
pairs(em.na)
pred_na <- ggpredict(lm.na, terms="fire") #predictions
names(pred_na) <- c("fire", "predicted", "std.error", "conf.low", "conf.high")

plot_na<-
  ggplot(pred_na, aes(x = fire, y = predicted, colour=fire, group=fire)) + #
  geom_point(position=position_dodge(0.6), size=3) +
  geom_point(data = abiotic_data_2, aes(x = fire, y = Na.BS , fill = fire), position = position_dodge(width = 0.6), alpha = 0.5, size=2) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, width = 0.3), position=position_dodge(0.6), linewidth=1) +
  scale_color_manual(values=c("#FF9900", "#00b0f0", "#6633FF", "#04c728"), labels=c("control", "spring", "summer", "unburnt")) +
  scale_x_discrete(limits=c("UB","CO","SP", "SU"), breaks=c("UB","CO","SP", "SU"), labels=c("unburnt", "control", "spring", "summer"))+
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(colour="black", face="italic", size=12)) + 
  theme(axis.text.y = element_text(colour="black", size=12)) +
  labs(x= "Fire history", y=expression(Sodium~"("~italic("% of CEC")~")")) +
  theme(aspect.ratio=1) +
  theme(legend.position="none") 

