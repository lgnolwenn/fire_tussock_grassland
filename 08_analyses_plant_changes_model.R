# ANALYSES OF HMSC MODEL OF THE EFFECTS OF FIRE HISTORY ON CHANGES OF THE PLANT COMMUNITY OF DEEP STREAM
## Code for analyzing the results of the HMSC models of changes of plant community.
## Note: paths need to be changed to match your device.

#LOAD PACKAGES

library(tidyverse)
library(Hmsc)
library(corrplot)
library(ggtext)
library(reshape2)
library(glue)
library(ggcorrplot)
library(ggpubr) 
library(gridGraphics)
library(grid)

# LOAD DATA

load("C:/Users/nolwe/Documents/M2/internship/project/model/models_plants_changes_3.RData")
m=temporal_models_plants[[3]] #change the index in function of the model you are
# interested with, only the 2-year-post-fire period was interesting for recovery.
plant_traits<-read.csv("C:/Users/nolwe/Documents/M2/internship/project/data/plant_traits.csv")

# CONVERGENCE AND MODEL FIT

## Get the posteriors
mpost = convertToCodaObject(m)
## Plot the posteriors distribution with the potential reduction factors of each parameter
par(mfrow=c(3,2))
ess.beta = effectiveSize(mpost$Beta)
psrf.beta = gelman.diag(mpost$Beta, multivariate=FALSE)$psrf
hist(ess.beta)
hist(psrf.beta)
ess.gamma = effectiveSize(mpost$Gamma)
psrf.gamma = gelman.diag(mpost$Gamma, multivariate=FALSE)$psrf
hist(ess.gamma)
hist(psrf.gamma)
sppairs = matrix(sample(x = 1:10^2, size = 25))
tmp = mpost$Omega[[1]]
for (chain in 1:length(tmp)){
  tmp[[chain]] = tmp[[chain]][,sppairs]
}
ess.omega = effectiveSize(tmp)
psrf.omega = gelman.diag(tmp, multivariate=FALSE)$psrf
hist(ess.omega)
hist(psrf.omega)

## Evaluate model fit
par(mfrow=c(1,1))
pred=computePredictedValues(m)
MF = evaluateModelFit(hM=m, predY=pred)
hist(MF$RMSE, xlim = c(0,1), main=paste0("Mean = ", round(mean(MF$RMSE),2)))
hist(MF$R2, xlim = c(0,1), main=paste0("Mean = ", round(mean(MF$R2),2)))

## Chains
par(mar = c(1, 1, 1, 1))
plot(mpost$Beta)

# VARIANCE PARTITIONING

VP = computeVariancePartitioning(m)
## we also used "pred" and "MF"

## Variance data frame
variance<-data.frame(species=m$spNames, total=rep(1, 37), explained=MF$R2)
vpval=VP$vals
vpval=as.data.frame(vpval)
vpval=t(vpval)
vpval=as.data.frame(vpval)
variance$fire=vpval$treatment # variance explained by fire (our environmental covariate)
variance$random=vpval$`Random: plot`# variance explained by random effect of plot

## Create 3 variance data frames: unexplained, fire and random
unexplained_variance<-data.frame(species=variance$species, partition=rep("unexplained", 37), value=variance$total-variance$explained)
unexplained_variance<-arrange(unexplained_variance, value)
species_order<-unexplained_variance$species
fire_variance<-data.frame(species=variance$species, partition=rep("fire", 37), value=variance$fire*variance$explained)
fire_variance<-fire_variance[match(species_order, fire_variance$species),]
random_variance<-data.frame(species=variance$species, partition=rep("random effect of plot", 37), value=variance$random*variance$explained)
random_variance<-random_variance[match(species_order, random_variance$species),]
variance_plot<-rbind(unexplained_variance, fire_variance, random_variance)
## Combine the 3 previous data frame in a data frame that will be plot
variance_plot<-rbind(unexplained_variance, fire_variance, random_variance)
### As the model is perform at the species level, explained variation is computed at the species level
### We need to average them and take the standard deviation for representativness of the model variance
### partitioning.
uvar=round(mean(unexplained_variance$value)*100, digits=1) #get mean value of unexplained variance
uvar_sd=round(sd(unexplained_variance$value)*100, digits=1) #get standard deviation of unexplained variance
fvar=round(mean(fire_variance$value)*100, digits=1) #get mean value of variance explained by fire
fvar_sd=round(sd(fire_variance$value)*100, digits=1)#get standard deviation of variance explained by fire
rvar=round(mean(random_variance$value)*100,  digits=1) #get mean value of variance explained by plot
rvar_sd=round(sd(random_variance$value)*100, digits=1)#get standard deviation of variance explained by plot

p_VP<-ggplot(variance_plot, aes(x=fct_inorder(species), y=value, fill=partition))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("#FF9900", "#339966", "#999999"), name  ="Variance partition",
                    breaks=c("fire", "random effect of plot", "unexplained"),
                    labels=c(glue("Fire ({fvar} ± {fvar_sd} %)"), 
                             glue("Random effect of plot ({rvar} ± {rvar_sd}%)"), 
                             glue("Unexplained ({uvar} ± {uvar_sd}%)")))
p_VP<-p_VP+theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=7), 
                 axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(), 
                 panel.background = element_blank(), 
                 axis.title = element_blank())
p_VP
### how much do traits explain out of the variation among the species in their responses to environmental covariates.
VP$R2T$Beta

# ENVIRONMENTAL NICHE

postBeta = getPostEstimate(m, parName="Beta")

##Get species and covariate names from model
spNames<- m$spNames
covNames<-m$covNames

## Pull the data from postBeta for plotting
mbeta = postBeta$mean
betaP = postBeta$support

## Choose the support level for Beta plot
supportLevel <- 0.95

## Format data as mean for all values with .95 support level
toPlot_b = mbeta
toPlot_b = toPlot_b * ((betaP > supportLevel) + (betaP < (1 - supportLevel)) > 0)

## Format the data as matrix and add column and row names
betaMat = matrix(toPlot_b, nrow = m$nc, ncol = ncol(m$Y))
colnames(betaMat)<- spNames
rownames(betaMat) <- covNames
betaMat<- as.data.frame(betaMat) %>%
  dplyr::slice(-1) #remove the intercept

## Update covariate names
rownames(betaMat) <-c('spring', 'summer')
betaMat<- as.data.frame(betaMat) 
betaMat<-betaMat[,!colSums(betaMat)==0]#remove unsignificant species
betaMat[betaMat == 0] <- NA # replace 0 by NA

## Reorder species, so displays highest values for unburnt on top
betaMat<-betaMat[, order(unname(unlist(as.list(betaMat[1,]))), na.last=FALSE)]

## Reformat for ggplot
betaMatmelt<-as.data.frame(melt(as.matrix(betaMat)))

## Color plants by traits
betaMatmelt$Var2<-gsub(".", " ", betaMatmelt$Var2, fixed=TRUE)
GrowthForm<-c()
BioStatus<-c()
for (i in 1:nrow(betaMatmelt)){
  sp=as.character(betaMatmelt[i,"Var2"])
  traits_sp=as.data.frame(plant_traits%>% filter(PreferredName==sp))
  BioStatus[i]<-traits_sp$BioStatus
  GrowthForm[i]<-traits_sp$GrowthForm
}
betaMatmelt$BioStatus=BioStatus
betaMatmelt$GrowthForm=GrowthForm
### Add a color code for each trait
col_code<-c()
for (i in 1:nrow(betaMatmelt)){
  row_i=betaMatmelt[i,]
    if (row_i$BioStatus=="Indigenous"& row_i$GrowthForm=="Forb"){col_code[i]="#00BB00"}
    else if (row_i$BioStatus=="Exotic"& row_i$GrowthForm=="Forb"){col_code[i]="#00BB00"
    betaMatmelt[i, "Var2"]= paste0(row_i$Var2, "*")}
    else if (row_i$BioStatus=="Indigenous"& row_i$GrowthForm=="Graminoid"){col_code[i]="#008600"}
    else if (row_i$BioStatus=="Exotic"& row_i$GrowthForm=="Graminoid"){col_code[i]="#008600"
    betaMatmelt[i, "Var2"]= paste0(row_i$Var2, "*")}
    else if (row_i$BioStatus=="Indigenous"& row_i$GrowthForm=="Woody"){col_code[i]="#005000"}
    else if (row_i$BioStatus=="Exotic"& row_i$GrowthForm=="Woody"){col_code[i]="#005000"
    betaMatmelt[i, "Var2"]= paste0(row_i$Var2, "*")}
}
betaMatmelt$ColTraits<-col_code
### Create a name variable to color the species in the plot
betaMatmelt=betaMatmelt%>% mutate(
  name = glue("<i style='color:{ColTraits}'>{Var2}</i>"))

## Set the colors of the plot
colors = colorRampPalette(c('#006666', '#E5FFFF', '#FFE5CB', '#662700'))
colorLevels<-3
cols_to_use= colors(colorLevels)

## Plot
p_time_env<-ggplot(betaMatmelt, aes(x = Var1, y = fct_inorder(name), fill = value)) +
  labs(x = "Fire history", y = "Species", fill = "") +
  geom_tile(color = 'gray60')+
  scale_fill_steps2(n.breaks = 40,
                    high = cols_to_use[3],
                    mid = 'white',
                    low = cols_to_use[1],
                    nice.breaks = TRUE,
                    na.value = "white",
                    midpoint=0,
                    limits=c(-2.0, 2.0)
  ) +
  theme(plot.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = NA),
        legend.margin = margin(l = 0.5, unit = 'cm'),
        legend.title = element_text(hjust = 0.08, size=50, face="bold"),
        legend.key.width = unit(2, 'cm'),
        legend.key.height = unit(6, 'cm'),
        legend.text = element_text(size = 45),
        axis.text.x = element_text(size=50),
        axis.text.y = element_markdown(size=50),
        axis.title.x = element_text(size=50, face="bold"),
        axis.title.y = element_text(size=50, face="bold"))+
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0))

# TRAITS CONTRIBUTION TO ENVIRONMENTAL NICHE

postGamma = getPostEstimate(m, parName="Gamma")

## Get traits and covariate names from model
trNames<- m$trNames
covNames<-m$covNames

## Pull the data from postGamma for plotting
mgamma = postGamma$mean
gammaP = postGamma$support

## Choose the support level for gamma plot
supportLevel <- 0.95

## Format data as mean for all values with .95 support level
toPlot_g = mgamma
toPlot_g = toPlot_g * ((gammaP > supportLevel) + (gammaP < (1 - supportLevel)) > 0)
###replace zeros with NA
toPlot_g[toPlot_g == 0] <- NA

## Format the data as matrix and add column and row names
gammaMat = matrix(toPlot_g, nrow = m$nc, ncol = m$nt)
colnames(gammaMat)<- trNames
rownames(gammaMat) <- covNames

## Update covariate names
rownames(gammaMat) <-c("intercept", 'spring', 'summer')
colnames(gammaMat)<-c("intercept", "exotic", "forb", "woody")

## Reformat for ggplot
gammaMatmelt<-as.data.frame(melt(as.matrix(gammaMat)))

###Set the colors
colors = colorRampPalette(c('#006666', '#E5FFFF', '#FFE5CB', '#662700'))
colorLevels<-3
cols_to_use= colors(colorLevels)

## Plot
p_time_tr<-ggplot(gammaMatmelt, aes(x = Var1, y = Var2, fill = value)) +
  labs(x = "Environmental Niche", y = "Traits", fill = "") +
  geom_tile(color = 'gray60')+
  scale_fill_steps2(n.breaks = 7,
                    high = cols_to_use[3],
                    mid = cols_to_use[2],
                    low = cols_to_use[1],
                    nice.breaks = TRUE,
                    na.value = "white",
  ) +
  theme(plot.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = NA),
        legend.margin = margin(l = 0.5, unit = 'cm'),
        legend.title = element_text(hjust = 0.08),
        legend.key.width = unit(0.5, 'cm'),
        legend.key.height = unit(1.25, 'cm'),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        #axis.text.y = element_markdown()
  )+
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0))



# SPECIES ASSOCIATIONS

OmegaCor = computeAssociations(m)
supportLevel = 0.95
toPlot_sa = ((OmegaCor[[1]]$support > supportLevel)
             + (OmegaCor[[1]]$support < (1-supportLevel))> 0)* OmegaCor[[1]]$mean
toPlot_sa=toPlot_sa[,!colSums(toPlot_sa)==1]
toPlot_sa=toPlot_sa[!rowSums(toPlot_sa)==0,]

## Color species by traits (plants) or phyla (fungi)
### Pull out significant species
plants_cor<-colnames(toPlot_sa)
### Plants
#### Add information on plant traits
BioStatus<-c()
GrowthForm<-c()
i=0
for (sp in plants_cor){
  sp<-gsub(".", " ", sp, fixed = TRUE)
  traits_sp=as.data.frame(plant_traits%>% filter(PreferredName==sp))
  i=i+1
  BioStatus[i]<-traits_sp$BioStatus
  GrowthForm[i]<-traits_sp$GrowthForm
}
df_plant_colors_cor<-data.frame(species=plants_cor, BioStatus=BioStatus, GrowthForm=GrowthForm, name=gsub(".", " ", plants_cor, fixed = TRUE))
#### Add a color code for each trait
col_code<-c()
for (i in 1:nrow(df_plant_colors_cor)){
  row_i=df_plant_colors_cor[i,]
  if (row_i$BioStatus=="Indigenous"& row_i$GrowthForm=="Forb"){col_code[i]="#00BB00"}
  else if (row_i$BioStatus=="Exotic"& row_i$GrowthForm=="Forb"){col_code[i]="#00BB00"
  df_plant_colors[i, "name"]= paste0(row_i$name, "*")}
  else if (row_i$BioStatus=="Indigenous"& row_i$GrowthForm=="Graminoid"){col_code[i]="#008600"}
  else if (row_i$BioStatus=="Exotic"& row_i$GrowthForm=="Graminoid"){col_code[i]="#008600"
  df_plant_colors[i, "name"]= paste0(row_i$name, "*")}
  else if (row_i$BioStatus=="Indigenous"& row_i$GrowthForm=="Woody"){col_code[i]="#005000"}
  else if (row_i$BioStatus=="Exotic"& row_i$GrowthForm=="Woody"){col_code[i]="#005000"
  df_plant_colors[i, "name"]= paste0(row_i$name, "*")}
}
df_plant_colors_cor$col_code<-col_code

### Prepare all names and colors to use on the plot
df_cor<-df_plant_colors_cor[,4:5]
df_cor=df_cor%>% mutate(
  nam_col= glue("<i style='color:{col_code}'>{name}</i>"))
colnames(toPlot_sa)<-df_cor$name
rownames(toPlot_sa)<-df_cor$name
color_cor<-df_cor$col_code
color_cor_use<-c(color_cor[2], color_cor[6], color_cor[3], color_cor[7],
                 color_cor[5], color_cor[1], color_cor[4])

corrplot(toPlot_sa, 
         method = "color", 
         type="upper", 
         order="hclust", 
         tl.col=color_cor_use, # color species names
         col = rev(COL2('RdBu', 20))) # color of the correlation scale
grid.echo()
p_species <- grid.grab() #save the plot

# LEGEND FOR PLOTS

cat<-c("Plant traits", "Forb", "Graminoid","Woody", "* Exotic plant")
col_code<-c("black", "#00BB00","#008600","#005000", "black")
df_legend<-data.frame(cat, col_code)
df_legend$x<-rep(1, 5)
df_legend$y<-c(15, 12, 9, 6, 3)
p_text<-ggplot(df_legend, aes(x, y)) + 
  geom_text(aes(label=cat), color=col_code, size=c(14,12,12,12,12))+
  xlim(c(0.98,1.02))+
  theme_void()

