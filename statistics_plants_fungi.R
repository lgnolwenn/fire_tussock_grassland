# ANALYSES OF PLANT AND FUNGAL DATASETS
## Basic community statistics on plant and fungal data from Deep Stream grassland.
## Note: paths need to be changed to match your device.

# LOAD PACKAGES

library(dplyr)
library(ggplot2)
library(vegan)
library(ggordiplots)
library(ggpubr)
library(lme4)
library(emmeans)
library(forcats)
library(ggeffects)


# LOAD DATA

fungi_it<-read.csv("C:/Users/nolwe/Documents/M2/internship/project/data/fungi_inter_tussock_final_selection.csv")
fungi_ca <- read.csv("C:/Users/nolwe/Documents/M2/internship/project/data/fungi_chionochloa_final_selection.csv")
plants_allt <- read.csv("C:/Users/nolwe/Documents/M2/internship/project/data/temporal_plant_data.csv")

# ARRANGING DATA AND CREATING COMMUNITY MATRICES

## Plants all times
colnames(plants_allt)[1]="subplot_id"
plants_allt=arrange(plants_allt, desc(subplot_id))
plants_allt$time<-substr(plants_allt$subplot_id, 1, 2)
plants_allt$time<-gsub('_', '', plants_allt$time)
plants_allt$treatment<-substr(plants_allt$subplot_id, 3, 5)
plants_allt$treatment<-gsub('_', '', plants_allt$treatment)
plants_allt$subplot<-substr(plants_allt$subplot_id, 6, 12)
plants_allt$subplot<-gsub('_', '', plants_allt$subplot)
plants_allt$subplot<-gsub('O', '_', plants_allt$subplot)
for (i in 1:nrow(plants_allt)){
  x=plants_allt[i,]
  if (x$treatment=="UB"){
    a=x$subplot
    b=substr(a, 1, 1)
    b=as.numeric(b)+9
    c=substr(a, 2, 3)
    new_a=paste0(b, c)
    plants_allt[i, "subplot"]<-new_a
  }}
plants_allt$plot<-substr(plants_allt$subplot, 1, 2)
plants_allt$plot<-gsub('_', '', plants_allt$plot)
plants_allt<-plants_allt[,2:49]

### Community matrix
plants_allt_sp=plants_allt %>% select( -"treatment", -"plot", -"subplot", -"time")
plants_allt_sp_a<-colnames(plants_allt_sp) #get the list of species
plants_allt_com_a<-plants_allt$subplot #get the list of subplots
m_plants_allt=matrix(, nrow=length(plants_allt_com_a), ncol=length(plants_allt_sp_a), 
                     dimnames=list(plants_allt_com_a, plants_allt_sp_a)) #create a community matrix
#fill the matrix with the percent of cover
for (i in 1:length(plants_allt_com_a)){#for each subplot
  for(j in 1:length(plants_allt_sp_a)){ #for each species
    m_plants_allt[i,j]= plants_allt_sp[i, j]/100
  }
}

## Fungi IT 
fungi_it<-fungi_it[,2:105]
colnames(fungi_it)[2]<-"subplot"
fungi_it$plot<-substr(fungi_it$subplot, 1, 2)
fungi_it$plot<-gsub('_', '', fungi_it$plot)
### Community matrix
fungi_it_sp=fungi_it %>% select( -"treatment", -"plot", -"subplot")
fungi_it_sp_a<-colnames(fungi_it_sp) #get the list of species
fungi_it_com_a<-fungi_it$subplot #get the list of subplots
m_fungi_it=matrix(, nrow=length(fungi_it_com_a), ncol=length(fungi_it_sp_a), 
                   dimnames=list(fungi_it_com_a, fungi_it_sp_a)) #create a community matrix
#fill the matrix with the percent of cover
for (i in 1:length(fungi_it_com_a)){#for each subplot
  for(j in 1:length(fungi_it_sp_a)){ #for each species
    m_fungi_it[i,j]= fungi_it_sp[i, j]
  }
}

## Fungi CA 
fungi_ca<-fungi_ca[,2:100]
colnames(fungi_ca)[2]<-"subplot"
fungi_ca$plot<-substr(fungi_ca$subplot, 1, 2)
fungi_ca$plot<-gsub('_', '', fungi_ca$plot)
### Community matrix
fungi_ca_sp=fungi_ca %>% select( -"treatment", -"plot", -"subplot")
fungi_ca_sp_a<-colnames(fungi_ca_sp) #get the list of species
fungi_ca_com_a<-fungi_ca$subplot #get the list of subplots
m_fungi_ca=matrix(, nrow=length(fungi_ca_com_a), ncol=length(fungi_ca_sp_a), 
                  dimnames=list(fungi_ca_com_a, fungi_ca_sp_a)) #create a community matrix
#fill the matrix with the percent of cover
for (i in 1:length(fungi_ca_com_a)){#for each subplot
  for(j in 1:length(fungi_ca_sp_a)){ #for each species
    m_fungi_ca[i,j]= fungi_ca_sp[i, j]
  }
}

## All fungi
asv_c<-colnames(fungi_ca)
asv_it<-colnames(fungi_it)
asv<-unique(c(asv_c, asv_it))
fungi<-data.frame()
for (i in asv){
  for (j in 1:nrow(fungi_ca)){
    if (i %in% asv_c){
      fungi[j, i]=fungi_ca[j,i]
    }
    else{
      fungi[j, i]=0
    }
  }
}

for (i in asv){
  for (j in 35:(34+nrow(fungi_it))){
    if (i %in% asv_it){
      fungi[j, i]=fungi_it[j-34,i]
    }
    else{
      fungi[j, i]=0
    }
  }
}

fungi$X<-c(rep("C", 34), rep("IT", 34))
### Community matrix
fungi_asv=fungi %>% select(-"plot", -"subplot", -"treatment", -"X")
fungi_sp_a<-colnames(fungi_asv) #get the list of species
fungi_com_a<-fungi$subplot #get the list of subplots
m_f=matrix(, nrow=length(fungi_com_a), ncol=length(fungi_sp_a), 
           dimnames=list(fungi_com_a, fungi_sp_a)) #create a community matrix
#fill the matrix with the percent of cover
for (i in 1:length(fungi_com_a)){#for each subplot
  for(j in 1:length(fungi_sp_a)){ #for each species
    m_f[i,j]= fungi_asv[i, j]
  }
}

# COMPARING INTERTUSSOCK (IT) AND CHIONOCHLOA-ADJACENT (CA) DATA SETS

## Get a similarity matrix based on Bray-Curtis distances 
md_f<-vegdist(m_f, "bray")

## Calculate the PCoA ordination
PCoA_f<-wcmdscale(md_f, eig=TRUE)
pcoa.scores.f1<-scores(PCoA_f, display="sites", choices=1)
pcoa.scores.f2<-scores(PCoA_f, display="sites", choices=2)
percent_exp_PCoA_f<- round((100* PCoA_f$eig / sum(PCoA_f$eig)), digits = 1)

## Create a data frame
df_ord<-data.frame(fungi$X, fungi$treatment, fungi$subplot, pcoa.scores.f1, pcoa.scores.f2)
colnames(df_ord)=c("soil.origin", "treatment", "subplot", "pcoa1", "pcoa2")

## Plot
plot_PCoA_f<-gg_ordiplot(PCoA_f, groups=fungi$X, pt.size = 3, choices=c(1,2))
plot_PCoA_f<- plot_PCoA_f$plot
plot_PCoA_f<-plot_PCoA_f + theme_classic() +
  xlab(glue("PCO1 ({percent_exp_PCoA_f[1]}%)"))+
  ylab(glue("PCO2 ({percent_exp_PCoA_f[2]}%)"))

## Welch t.test of the sites scores to soil origin
t.test(df_ord$pcoa1~df_ord$soil.origin)
### the test is non significant, IT and CA data sets are similar. We will use
### only IT data set in the following analyses.

# RICHNESS

## Fungi IT

richness_fungi_it=rowSums(fungi_it!=0)[rowSums(fungi_it!=0)!=0]#count the number of fungal taxa
fungi_it_rich= data.frame(subplot=fungi_it$subplot, plot= fungi_it$plot, 
                          treatment=fungi_it$treatment, richness=richness_fungi_it)
### Linear model
lmer_fungi_it_rich=lmer( richness~ treatment+(1|plot), data=fungi_it_rich)
summary(lmer_fungi_it_rich)
### Contrasts between fire histories
emm_fungi_it_rich<-emmeans(lmer_fungi_it_rich, "treatment")
pairs(emm_fungi_it_rich)
### Get means and standard deviations
mean_richness_fungi_it= fungi_it_rich %>%
  group_by(treatment) %>%
  summarise_at(vars(richness), list(mean = mean, sd=sd))
### Compute predictions
emm_fungi_it_rich_2<-emmeans(lmer_fungi_it_rich, specs = pairwise ~ treatment)
table_emms_fungi_rich <- as.data.frame(emm_fungi_it_rich_2$contrast) %>%
  mutate(organism = "Soil fungi", variable = "Richness")
### Get predicted values for plotting
pred_fungi_rich <- ggpredict(lmer_fungi_it_rich, terms = "treatment")
names(pred_fungi_rich) <- c("treatment", "predicted", "std.error", "conf.low", "conf.high")
### Plot
plot_fungi_rich <- 
  ggplot(pred_fungi_rich, aes(x = treatment, y = predicted, colour=treatment, group=treatment)) + #
  geom_point(position=position_dodge(0.6), size=3, aes(shape=treatment)) +
  geom_point(data = fungi_it_rich, aes(x = treatment, y = richness , 
                                       fill = treatment, shape=treatment), 
             position = position_dodge(width = 0.6), alpha = 0.5, size=2) +
  scale_shape_manual(values=c(15, 15, 15, 15)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, width = 0.3), 
                position=position_dodge(0.6), linewidth=1) +
  scale_color_manual(values=c("#FF9900", "#00b0f0", "#6633FF", "#04c728"), 
                     labels=c("control", "spring", "summer", "unburnt")) +
  scale_x_discrete(limits=c("UB","CO","SP", "SU"), breaks=c("UB","CO","SP", "SU"), 
                   labels=c("unburnt", "control", "spring", "summer"))+
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(colour="black", face="italic", size=12)) + 
  theme(axis.text.y = element_text(colour="black", size=12)) +
  ylab("Species richness") +
  xlab(" Fire history") +
  theme(aspect.ratio=1) +
  theme(legend.position="none") 

## Plants all times

richness_plants_allt=rowSums(plants_allt!=0)[rowSums(plants_allt!=0)!=0]#count the number of plant taxa
plants_allt_rich= data.frame(time= plants_allt$time, subplot=plants_allt$subplot, 
                             plot= plants_allt$plot, treatment=plants_allt$treatment, 
                             richness=richness_plants_allt)
plants_allt_rich$time=as.factor(plants_allt_rich$time)
plants_allt_rich$time=relevel(plants_allt_rich$time, "2")
plants_allt_rich$treatment=as.factor(plants_allt_rich$treatment)
plants_allt_rich$treatment=relevel(plants_allt_rich$treatment, "UB")
### Removing unburnt plots for the statistical analysis
plants_allt_rich_2=arrange(plants_allt_rich, desc(treatment))
plants_allt_rich_2=plants_allt_rich_2[1:81,]
plants_allt_rich_2$time=as.factor(plants_allt_rich_2$time)
plants_allt_rich_2$time=relevel(plants_allt_rich_2$time, "2")
plants_allt_rich_2$treatment=as.factor(plants_allt_rich_2$treatment)
plants_allt_rich_2$treatment=relevel(plants_allt_rich_2$treatment, "CO")
### Linear model
lmer_plants_allt_rich_2=lmer( richness~ treatment*time+(1|plot/subplot), data=plants_allt_rich_2)
summary(lmer_plants_allt_rich_2)
### Contrasts between fire histories
emm_plants_allt_rich_2_t_tr<-emmeans(lmer_plants_allt_rich_2,  ~ time | treatment)
emm_plants_allt_rich_2_tr_t<-emmeans(lmer_plants_allt_rich_2,  ~ treatment | time)
pairs(emm_plants_allt_rich_2_t_tr)
pairs(emm_plants_allt_rich_2_tr_t)
### Get means and standard deviations
 mean_richness_plants_allt= plants_allt_rich_2 %>%
   group_by( time) %>%
  summarise_at(vars(richness), list(mean_rich = mean, sd_rich=sd))
### Compute predicted values
 emm_plants_rich_2<-emmeans(lmer_plants_allt_rich_2, specs = pairwise ~ treatment:time)
 table_emms_plants_rich <- as.data.frame(emm_plants_rich_2$contrast) %>%
   mutate(organism = "Plants", variable = "Richness")
### Get predicted values for plotting
 pred_plants_rich <- ggpredict(lmer_plants_allt_rich_2, terms = c("treatment", "time"))
 names(pred_plants_rich) <- c("treatment", "predicted", "std.error", "conf.low", "conf.high", "time")
#### For unburnt plots, manually get the needed values
 ub_plants_rich<-plants_allt_rich%>% filter(treatment=="UB"& time=="2")
pred_ub<-as.numeric(mean(ub_plants_rich$richness))
std.error <- function(x) sd(x)/sqrt(length(x))
se_ub<-as.numeric(std.error(ub_plants_rich$richness))
margin_ub <- as.numeric(qt(0.975,df=9-1)*sd(ub_plants_rich$richness)/sqrt(9))
conf_l_ub<-as.numeric(pred_ub-margin_ub)
conf_h_ub<-as.numeric(pred_ub+margin_ub)
ub_predicted<-data.frame(treatment="UB", predicted=pred_ub, std.error=se_ub, 
                         conf.low=conf_l_ub, conf.high=conf_h_ub, time="2")
levels(pred_plants_rich$treatment)<-c(levels(pred_plants_rich$treatment), "UB")
#### Combine predicted values
pred_plants_rich <- rbind( pred_plants_rich, ub_predicted)
### Combine richness data sets
plants_allt_rich_3<-rbind(plants_allt_rich_2, ub_plants_rich)
### Plot
plot_plants_rich <- 
  ggplot(pred_plants_rich, aes(x = treatment, y = predicted, colour=treatment, group=time)) + #
  geom_point(position=position_dodge(0.6), size=3, aes(shape=time)) +
  geom_point(data = plants_allt_rich_3, 
             aes(x = treatment, y = richness , fill = treatment, shape=time), 
             position = position_dodge(width = 0.6), alpha = 0.5, size=2) +
  scale_shape_manual(values=c(15, 19, 17)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, width = 0.3), 
                position=position_dodge(0.6), linewidth=1) +
  scale_color_manual(values=c("#FF9900", "#00b0f0", "#6633FF", "#04c728"), 
                     labels=c("control", "spring", "summer", "unburnt")) +
  scale_x_discrete(limits=c("UB","CO","SP", "SU"), breaks=c("UB","CO","SP", "SU"), 
                   labels=c("unburnt", "control", "spring", "summer"))+
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(colour="black", face="italic", size=12)) + 
  theme(axis.text.y = element_text(colour="black", size=12)) +
  ylab("Species richness") +
  xlab(" Fire history") +
  theme(aspect.ratio=1) +
  theme(legend.position="none") 


# DIVERSITY

## Fungi IT

diversity_fungi_it=diversity(m_fungi_it,index = "simpson")
fungi_it_div= data.frame(subplot=fungi_it$subplot, plot= fungi_it$plot, 
                         treatment=fungi_it$treatment, diversity=diversity_fungi_it)
### Linear model
lmer_fungi_it_div=lmer( diversity~ treatment+(1|plot), data=fungi_it_div)
summary(lmer_fungi_it_div)
### Contrasts between fire histories
emm_fungi_it_div<-emmeans(lmer_fungi_it_div, "treatment")
pairs(emm_fungi_it_div)
### Get means and standard deviations
mean_diversity_fungi_it= fungi_it_div %>%
  group_by(treatment) %>%
  summarise_at(vars(diversity), list(name = mean))
### Compute predicted values
emm_fungi_it_div_2<-emmeans(lmer_fungi_it_div, specs = pairwise ~ treatment)
table_emms_fungi_div <- as.data.frame(emm_fungi_it_div_2$contrast) %>%
  mutate(organism = "Soil fungi", variable = "diversity")
### Get predicted values for plotting
pred_fungi_div <- ggpredict(lmer_fungi_it_div, terms = "treatment")
names(pred_fungi_div) <- c("treatment", "predicted", "std.error", "conf.low", "conf.high")
### Plot
plot_fungi_div <- 
  ggplot(pred_fungi_div, aes(x = treatment, y = predicted, colour=treatment, group=treatment)) + #
  geom_point(position=position_dodge(0.6), size=3, aes(shape=treatment)) +
  geom_point(data = fungi_it_div, aes(x = treatment, y = diversity , 
                                      fill = treatment, shape=treatment), 
             position = position_dodge(width = 0.6), alpha = 0.5, size=2) +
  scale_shape_manual(values=c(15, 15, 15, 15)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, width = 0.3), 
                position=position_dodge(0.6), linewidth=1) +
  scale_color_manual(values=c("#FF9900", "#00b0f0", "#6633FF", "#04c728"), 
                     labels=c("control", "spring", "summer", "unburnt")) +
  scale_x_discrete(limits=c("UB","CO","SP", "SU"), breaks=c("UB","CO","SP", "SU"), 
                   labels=c("unburnt", "control", "spring", "summer"))+
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(colour="black", face="italic", size=12)) + 
  theme(axis.text.y = element_text(colour="black", size=12)) +
  ylab("Simpson diversity") +
  xlab(" Fire history") +
  theme(aspect.ratio=1) +
  theme(legend.position="none") 

## Plants all times

diversity_plants_allt=diversity(m_plants_allt,index = "simpson")
plants_allt_div= data.frame(time= plants_allt$time, subplot=plants_allt$subplot, 
                            plot= plants_allt$plot, treatment=plants_allt$treatment, 
                            diversity=diversity_plants_allt)
plants_allt_div$time=as.factor(plants_allt_div$time)
plants_allt_div$time=relevel(plants_allt_div$time, "2")
plants_allt_div$treatment=as.factor(plants_allt_div$treatment)
plants_allt_div$treatment=relevel(plants_allt_div$treatment, "UB")
### Removing unburnt plots for the analyses
plants_allt_div_2=plants_allt_div[10:90,]
plants_allt_div_2$time=as.factor(plants_allt_div_2$time)
plants_allt_div_2$time=relevel(plants_allt_div_2$time, "2")
plants_allt_div_2$treatment=as.factor(plants_allt_div_2$treatment)
plants_allt_div_2$treatment=relevel(plants_allt_div_2$treatment, "CO")
### Linear model
lmer_plants_allt_div_2=lmer( diversity~ treatment*time+(1|plot/subplot), data=plants_allt_div_2)
summary(lmer_plants_allt_div_2)
emm_plants_allt_div_2_t_tr<-emmeans(lmer_plants_allt_div_2,  ~ time | treatment)
emm_plants_allt_div_2_tr_t<-emmeans(lmer_plants_allt_div_2,  ~ treatment | time)
pairs(emm_plants_allt_div_2_t_tr)
pairs(emm_plants_allt_div_2_tr_t)
dpa<-ggplot(plants_allt_div_2)+geom_boxplot(aes(x=treatment, y=diversity, fill=time))+theme_classic()
mean_diversity_plants_allt= plants_allt_div_2 %>%
  group_by(treatment, time) %>%
  summarise_at(vars(diversity), list(mean_div = mean, sd_div=sd))
### Compute predicted values
emm_plants_div_2<-emmeans(lmer_plants_allt_div_2, specs = pairwise ~ treatment:time)
table_emms_plants_div <- as.data.frame(emm_plants_div_2$contrast) %>%
  mutate(organism = "Plants", variable = "diversity")
### Get predicted values for plotting
pred_plants_div <- ggpredict(lmer_plants_allt_div_2, terms = c("treatment", "time"))
names(pred_plants_div) <- c("treatment", "predicted", "std.error", "conf.low", "conf.high", "time")
#### For unburnt plots, manually get the needed values
ub_plants_div<-plants_allt_div%>% filter(treatment=="UB"& time=="2")
pred_ub<-as.numeric(mean(ub_plants_div$diversity))
std.error <- function(x) sd(x)/sqrt(length(x))
se_ub<-as.numeric(std.error(ub_plants_div$diversity))
margin_ub <- as.numeric(qt(0.975,df=9-1)*sd(ub_plants_div$diversity)/sqrt(9))
conf_l_ub<-as.numeric(pred_ub-margin_ub)
conf_h_ub<-as.numeric(pred_ub+margin_ub)
ub_predicted<-data.frame(treatment="UB", predicted=pred_ub, std.error=se_ub, 
                         conf.low=conf_l_ub, conf.high=conf_h_ub, time="2")
levels(pred_plants_div$treatment)<-c(levels(pred_plants_div$treatment), "UB")
#### Combine both predictions
pred_plants_div <- rbind( pred_plants_div, ub_predicted)
### Combine both data sets
plants_allt_div_3<-rbind(plants_allt_div_2, ub_plants_div)
### Plot
plot_plants_div <- 
  ggplot(pred_plants_div, aes(x = treatment, y = predicted, colour=treatment, group=time)) +
  geom_point(position=position_dodge(0.6), size=3, aes(shape=time)) +
  geom_point(data = plants_allt_div_3, aes(x = treatment, y = diversity , 
                                           fill = treatment, shape=time), 
             position = position_dodge(width = 0.6), alpha = 0.5, size=2) +
  scale_shape_manual(values=c(15, 19, 17)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, width = 0.3), 
                position=position_dodge(0.6), linewidth=1) +
  scale_color_manual(values=c("#FF9900", "#00b0f0", "#6633FF", "#04c728"), 
                     labels=c("control", "spring", "summer", "unburnt")) +
  scale_x_discrete(limits=c("UB","CO","SP", "SU"), breaks=c("UB","CO","SP", "SU"), 
                   labels=c("unburnt", "control", "spring", "summer"))+
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(colour="black", face="italic", size=12)) + 
  theme(axis.text.y = element_text(colour="black", size=12)) +
  ylab("Simpson diversity") +
  xlab(" Fire history") +
  theme(aspect.ratio=1) +
  theme(legend.position="none") 


# PRINCIPAL COORDINATES ANALYSES

## Fungi IT

### Get a similarity matrix based on Bray-Curtis distances 
m_dist_fungi_it<-vegdist(m_fungi_it, "bray")
### Calculate the PCoA ordination
PCoA_fungi_it<-wcmdscale(m_dist_fungi_it, eig=TRUE)
fungi_it_pcoa1<-scores(PCoA_fungi_it, display="sites", choices=1)
fungi_it_pcoa2<-scores(PCoA_fungi_it, display="sites", choices=2)
### Create data frame
df_pcoa_fungi_it=data.frame(time_cat=rep("2", 34), fire_history=fungi_it$treatment, 
                            pcoa1=fungi_it_pcoa1, pcoa2=fungi_it_pcoa2)
colnames(df_pcoa_fungi_it)=c("time_cat", "fire_history", "pcoa1", "pcoa2")
df_pcoa_fungi_it=df_pcoa_fungi_it %>% mutate(time_fire=paste0(time_cat, "_", fire_history))
### Prepare data to plot
groups <- df_pcoa_fungi_it$time_fire # to compare fire histories
mod <- betadisper(m_dist_fungi_it, groups, type = "centroid") # centroids of each group
percent_exp_pcoa_f<- round((100* mod$eig / sum(mod$eig)), digits = 2) # Get the 
# percentage of variation explained by each axis
#### Create a data frame with centroids
df_pcoa_fungi_it_mod<-df_pcoa_fungi_it
vectors<-as.data.frame(mod$vectors)
df_pcoa_fungi_it_mod$pcoa1<-vectors$PCoA1
df_pcoa_fungi_it_mod$pcoa2<-vectors$PCoA2
df_pcoa_fungi_it_mod$distances<- mod$distances
dist_pcoa_fungi_it<-df_pcoa_fungi_it_mod %>%
  group_by(time_cat, fire_history, time_fire) %>%
  summarise(pcoa1 = mean(pcoa1), pcoa2 = mean(pcoa2), distances = mean(distances))
### Plot
pcoa_ellipse_fungi_it<- df_pcoa_fungi_it_mod %>%
  ggplot(aes(x = pcoa1, y = pcoa2)) +
  stat_ellipse(aes(fill = time_fire), alpha = 0.3, geom = "polygon", level = 0.75) + 
  geom_point(aes(color = time_fire, shape = time_cat), size = 3, alpha = 0.7) +
  geom_point(data = dist_pcoa_fungi_it, aes(x = pcoa1, y = pcoa2, color = time_fire, 
                                            shape = time_cat), size = 6) + 
  scale_shape_manual(values = c(15), labels = c("2 months")) +
  scale_color_manual(values = c("#FF9900", "#00b0f0", "#6633FF", "#04c728"), 
                     labels = c( "control", "spring", "summer", "unburnt")) +
  scale_fill_manual(values = c("#FF9900", "#00b0f0", "#6633FF","#04c728"), 
                    labels = c( "control", "spring", "summer", "unburnt"))+
  xlab(paste("PCoA 1 (", percent_exp_pcoa_f[1], "%)", sep = "")) + 
  ylab(paste("PCoA 2 (", percent_exp_pcoa_f[2], "%)", sep = ""))+
  coord_fixed()+
  theme_classic() + 
  theme(aspect.ratio = 1,
    panel.background = element_rect(colour='black'),
        axis.text.y=element_text(colour="black", size = 14),
        axis.text.x = element_text(colour="black", size = 14),
        axis.title.y=element_text(colour="black", size = 18, margin=margin(t=20)),
        axis.title.x = element_text(colour="black", size = 18, margin=margin(t=20)),
        strip.text = element_text(colour="black", size = 18),
        plot.subtitle = element_text(colour="black", size = 18, vjust = -6),
        legend.position = "none")




## Plants all time

#Get a similarity matrix based on Bray-Curtis distances 
m_dist_plants_allt<-vegdist(m_plants_allt, "bray")
# Calculate the PCoA ordination
PCoA_plants_allt<-wcmdscale(m_dist_plants_allt, eig=TRUE)
plants_allt_pcoa1<-scores(PCoA_plants_allt, display="sites", choices=1)
plants_allt_pcoa2<-scores(PCoA_plants_allt, display="sites", choices=2)
### Create data frame
df_pcoa_plants_allt=data.frame(time_cat=plants_allt$time, 
                               fire_history=plants_allt$treatment, 
                               pcoa1=plants_allt_pcoa1, pcoa2=plants_allt_pcoa2)
for(i in 1:nrow(df_pcoa_plants_allt)){
  df_i=df_pcoa_plants_allt[i,]
  if(df_i$fire_history=="UB"){
  df_pcoa_plants_allt[i,"time_cat"]="2_u"
}
}
colnames(df_pcoa_plants_allt)=c("time_cat", "fire_history", "pcoa1", "pcoa2")
df_pcoa_plants_allt=df_pcoa_plants_allt %>% mutate(time_fire=paste0(time_cat, "_", fire_history))
### Prepare data to plot
groups <- df_pcoa_plants_allt$time_fire # to compare fire histories
mod <- betadisper(m_dist_plants_allt, groups, type = "centroid") # centroids of each group
percent_exp_pcoa_p<- round((100* mod$eig / sum(mod$eig)), digits = 2)# Get the 
# percentage of variation explained by each axis
#### Create a data frame with centroids
df_pcoa_plants_allt_mod<-df_pcoa_plants_allt
vectors<-as.data.frame(mod$vectors)
df_pcoa_plants_allt_mod$pcoa1<-vectors$PCoA1
df_pcoa_plants_allt_mod$pcoa2<-vectors$PCoA2
df_pcoa_plants_allt_mod$distances<- mod$distances
dist_pcoa_plants_allt<-df_pcoa_plants_allt_mod %>%
  mutate(time_cat = fct_recode(time_cat, "2 months"= "2", "13 months" = "13","26 months" = "26", "2 months - unburnt" ="2_u")) %>%
  mutate(time_cat = fct_relevel(time_cat , "2 months", "13 months", "26 months", "2 months - unburnt"))   %>%
  group_by(time_cat, fire_history, time_fire) %>%
  summarise(pcoa1 = mean(pcoa1), pcoa2 = mean(pcoa2), distances = mean(distances))
### Plot
pcoa_ellipse_plants_allt<- df_pcoa_plants_allt_mod %>%
  mutate(time_cat = fct_recode(time_cat, "2 months"= "2", "13 months" = "13","26 months" = "26", "2 months - unburnt" ="2_u")) %>%
  mutate(time_cat = fct_relevel(time_cat , "2 months", "13 months", "26 months", "2 months - unburnt")) %>%
  ggplot(aes(x = pcoa1, y = pcoa2)) +
  stat_ellipse(aes(fill = time_fire), alpha = 0.3, geom = "polygon", level = 0.75) + 
  geom_point(aes(color = fire_history, shape = time_cat), size = 3, alpha = 0.7) +
  geom_point(data = dist_pcoa_plants_allt, aes(x = pcoa1, y = pcoa2, color = fire_history, shape = time_cat), size = 6) + 
  geom_path(data = dist_pcoa_plants_allt, aes(x = pcoa1, y = pcoa2, group = fire_history, color=fire_history), arrow=arrow(length=unit(0.7,"cm")), size= 1.5)+
  scale_shape_manual(values = c(15,19, 17, 15), labels = c("2 months", "2 months - unburnt")) +
  scale_color_manual(values = c("#FF9900", "#00b0f0", "#6633FF", "#04c728"), labels = c( "control", "spring", "summer", "unburnt")) +
  scale_fill_manual(values = c("#FF9900", "#00b0f0", "#6633FF","#FF9900", "#00b0f0", "#6633FF", "#04c728","#FF9900", "#00b0f0", "#6633FF"), labels = c( "control", "spring", "summer","control", "spring", "summer", "unburnt", "control", "spring", "summer"))+
  xlab(paste("PCoA 1 (", percent_exp_pcoa_p[1], "%)", sep = "")) + 
  ylab(paste("PCoA 2 (", percent_exp_pcoa_p[2], "%)", sep = ""))+
  coord_fixed()+
  theme_classic() + 
  theme(panel.background = element_rect(colour='black'),
        axis.text.y=element_text(colour="black", size = 14),
        axis.text.x = element_text(colour="black", size = 14),
        axis.title.y=element_text(colour="black", size = 18, margin=margin(t=20)),
        axis.title.x = element_text(colour="black", size = 18, margin=margin(t=20)),
        strip.text = element_text(colour="black", size = 18),
        plot.subtitle = element_text(colour="black", size = 18, vjust = -6),
        legend.position = "none")

# LEGEND

fire_history<-c("unburnt", "control", "spring", "summer", "control", "spring", 
                "summer", "control", "spring", "summer")
time<-c("t1 (2 months)","t1 (2 months)", "t1 (2 months)", "t1 (2 months)", 
        "t2 (13 months)","t2 (13 months)", "t2 (13 months)", "t3 (26 months)", 
        "t3 (26 months)", "t3 (26 months)")
col_code<-c("#04c728", "#FF9900", "#00b0f0", "#6633FF", "#FF9900", "#00b0f0", 
            "#6633FF", "#FF9900", "#00b0f0", "#6633FF")
shape_code<-c(22, 22, 22, 22, 21, 21, 21, 24, 24, 24)
pos.x<-c(1:10)
pos.y<-c(1:10)
df_legend<-data.frame(fire_history, col_code, time, shape_code, pos.x, pos.y)
p_text1<-ggplot(df_legend, aes(pos.x, pos.y)) + geom_point(aes(shape=time))+
  scale_shape_manual(values=c(15, 19, 17), name= "Time after fire")+
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  theme_classic()+ 
  theme(legend.text = element_text(colour="black", size=16),
        legend.title = element_text(colour="black", size=18, face="bold"), 
        legend.position = "bottom")

p_text2<-ggplot(df_legend, aes(pos.x, pos.y)) + geom_point(aes(color=fire_history))+
  scale_color_manual(values=c("#04c728", "#FF9900", "#00b0f0", "#6633FF"),
                     name="Fire history", label=c("unburnt", "control", "spring", "summer"))+
    guides(color = guide_legend(override.aes = list(size = 3))) +
    theme_classic()+ 
    theme(legend.text = element_text(colour="black", size=16),
          legend.title = element_text(colour="black", size=18, face="bold"), 
          legend.position = "bottom")
    
leg1<-get_legend(p_text1)
leg2<-get_legend(p_text2)

