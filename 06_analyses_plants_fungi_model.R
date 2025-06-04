# ANALYSES OF HMSC MODEL OF THE EFFECTS OF FIRE HISTORY ON PLANT AND SOIL FUNGAL COMMUNITIES OF DEEP STREAM
## Code for analysing the results of the HMSC model of both plant and fungal communities.
## Note: paths need to be changed to match your device.

# LOAD PACKAGES

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

load("C:/Users/nolwe/Documents/M2/internship/project/model/model_all_2.RData")
m=the_model
fungi_taxa<-read.csv("C:/Users/nolwe/Documents/M2/internship/project/data/selected_fungi_taxa.csv")
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

## Get posteriors for environmental niche
mpost_ESS_PSRF <- 
  data.frame(Species = colnames(m$Y),
             Parameter = rep(c("fire"), 138),
             ESS = as.data.frame(effectiveSize(mpost$Beta))[, 1],
             Point_Est = as.data.frame(gelman.diag(mpost$Beta, multivariate = FALSE)$psrf)[, 1],
             Upper_CI = as.data.frame(gelman.diag(mpost$Beta, multivariate = FALSE)$psrf)[, 2])
mpost_ESS_PSRF

## Chains
par(mar = c(1, 1, 1, 1))
plot(mpost$Beta)

# VARIANCE PARTITIONING

VP = computeVariancePartitioning(m)
## we also used "pred" and "MF"

## Variance data frame
variance<-data.frame(species=m$spNames, total=rep(1, 138), explained=MF$R2)
vpval=VP$vals
vpval=as.data.frame(vpval)
vpval=t(vpval) #transpose the matrix
vpval=as.data.frame(vpval)
variance$fire=vpval$treatment #variance explained by fire (environmental covariate)
variance$random=vpval$`Random: plot`#variance explained by plot (random level)

## Create 3 variance data frames: unexplained, fire and random
unexplained_variance<-data.frame(species=variance$species, partition=rep("unexplained", 138), value=variance$total-variance$explained)
unexplained_variance<-arrange(unexplained_variance, value)
species_order<-unexplained_variance$species # variable to ensure the 3 data frames have the same order
fire_variance<-data.frame(species=variance$species, partition=rep("fire", 138), value=variance$fire*variance$explained)
fire_variance<-fire_variance[match(species_order, fire_variance$species),]
random_variance<-data.frame(species=variance$species, partition=rep("random effect of plot", 138), value=variance$random*variance$explained)
random_variance<-random_variance[match(species_order, random_variance$species),]

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
### Species with fire explaining most of their variation
fire_variance[fire_variance$value>0.5, "species"]
mean(VP$R2T$Beta)
sd(VP$R2T$Beta)

# ENVIRONMENTAL NICHE

postBeta = getPostEstimate(m, parName="Beta")

## Get species and covariate names from model
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
rownames(betaMat) <-c('control', 'spring', 'summer')
betaMat<- as.data.frame(betaMat) 
betaMat<-betaMat[,!colSums(betaMat)==0] #remove non-significant species
betaMat[betaMat == 0] <- NA #replace 0 by NA

## Reorder species so displays highest unburnt values on top
betaMat<-betaMat[, order(unname(unlist(as.list(betaMat[1,]))), 
                         unname(unlist(as.list(betaMat[2,]))), 
                         unname(unlist(as.list(betaMat[3,]))))]

## Reformat for ggplot
betaMatmelt<-as.data.frame(melt(as.matrix(betaMat)))

##Color species by traits (plants) or phyla (fungi)
### Get the names of significant species 
significant_species_env<-as.character(unique(betaMatmelt$Var2))
fungi_env<-significant_species_env[grepl("ASV", significant_species_env)]
fungi_env_taxa<-fungi_taxa %>% filter(ASV %in% fungi_env) # get the taxa corresponding to the significant ASV
plants_env<-significant_species_env[! grepl("ASV", significant_species_env)]
### Plants
#### Add information on traits
BioStatus<-c()
GrowthForm<-c()
i=0
for (sp in plants_env){
  sp<-gsub(".", " ", sp, fixed = TRUE)
  traits_sp=as.data.frame(plant_traits%>% filter(PreferredName==sp))
  i=i+1
  BioStatus[i]<-traits_sp$BioStatus
  GrowthForm[i]<-traits_sp$GrowthForm
}
df_plant_colors<-data.frame(species=plants_env, BioStatus=BioStatus, GrowthForm=GrowthForm, name=gsub(".", " ", plants_env, fixed = TRUE))
#### Add a color code for each trait
col_code<-c()
for (i in 1:nrow(df_plant_colors)){
  row_i=df_plant_colors[i,]
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

df_plant_colors$col_code<-col_code
### Fungi
df_fungi_colors<-data.frame(asv= fungi_env_taxa$ASV, name= fungi_env_taxa$name, phylum=fungi_env_taxa$Phylum )
fcol_code<-c()
for (i in 1:nrow(df_fungi_colors)){
  row_i=df_fungi_colors[i,]
  if (is.na(row_i$phylum)){fcol_code[i]="#FFBBFF"}
  else if (row_i$phylum =="Ascomycota"){fcol_code[i]="#500050"}
  else if (row_i$phylum =="Basidiomycota"){fcol_code[i]="#860086"}
  else if (row_i$phylum == "Chytridiomycota"){fcol_code[i]="#BB00BB"}
  else if (row_i$phylum == "Entorrhizomycota"){fcol_code[i]="#F100F1"}
  else if (row_i$phylum=="Glomeromycota"){fcol_code[i]="#FF50FF"}
  else if (row_i$phylum=="Mortierellomycota"){fcol_code[i]="#FF86FF"}
}
df_fungi_colors$col_code<-fcol_code
### Create a name variable with the species name and its color
name_sp<-c()
colTraits<-c()
for (i in 1:nrow(betaMatmelt)){ #for each species in the beta matrix
  sp_i=as.character(betaMatmelt[i, "Var2"]) #get the name of the species
  if (grepl ("ASV", sp_i)){ #for fungi
    name=unlist(df_fungi_colors %>% filter(asv==sp_i) %>% select(name))
    col=unlist(df_fungi_colors %>% filter(asv==sp_i) %>% select(col_code))
  }
  else{ # for plants
    name=unlist(df_plant_colors %>% filter(species==sp_i) %>% select(name))
    col=unlist(df_plant_colors %>% filter(species==sp_i) %>% select(col_code))
  }
  name_sp<-c(name_sp,name)
  colTraits<-c(colTraits, col)
}

betaMatmelt$ColTraits<-colTraits
betaMatmelt$Species<-name_sp
betaMatmelt=betaMatmelt%>% mutate(
  name = glue("<i style='color:{ColTraits}'>{Species}</i>")) #add the name variable

## Divide the data frame to get a nicer plot with species with a positive response on the
## one hand, and species with a negative response on the other hand.
### Get values by treatment (CO, SP, SU)
co_values<-betaMatmelt %>% filter(Var1== "control") %>% select(value, Var2) %>% mutate(value=sign(value))
sp_values<-betaMatmelt %>% filter(Var1== "spring") %>% select(value, Var2) %>% mutate(value=sign(value))
su_values<-betaMatmelt %>% filter(Var1== "summer") %>% select(value, Var2) %>% mutate(value=sign(value))
### First differentiate for control
species_co_pos<-as.character(unlist(co_values%>% filter(value==1)%>% select(Var2)))
species_co_neg<-as.character(unlist(co_values%>% filter(value==-1)%>% select(Var2)))
species_co_na<-as.character(unlist(co_values%>% filter(is.na(value))%>% select(Var2)))
### Then for spring
species_sp_pos<-as.character(unlist(sp_values%>% filter(value==1)%>% select(Var2)))
species_sp_neg<-as.character(unlist(sp_values%>% filter(value==-1)%>% select(Var2)))
### Finally for summer
species_su_pos<-as.character(unlist(su_values%>% filter(value==1)%>% select(Var2)))
species_su_neg<-as.character(unlist(su_values%>% filter(value==-1)%>% select(Var2)))
### Divide beta data frame between positive, negative, NA for control positive for others
### and NA for control negative for others
betaMatmeltpos<-betaMatmelt%>% filter(Var2 %in% species_co_pos)
betaMatmeltneg<-betaMatmelt%>% filter(Var2 %in% species_co_neg)
betaMatmeltnapos<-betaMatmelt %>% filter (Var2 %in% species_co_na & (Var2 %in% species_sp_pos|Var2 %in% species_su_pos) )
betaMatmeltnaneg<-betaMatmelt %>% filter (Var2 %in% species_co_na & (Var2 %in% species_sp_neg|Var2 %in% species_su_neg) )
### Combine to get only two partitions of the data
betaneg=rbind(betaMatmeltneg, betaMatmeltnaneg)
betapos=rbind(betaMatmeltnapos,betaMatmeltpos)

## Choose the correlation colors
colors = colorRampPalette(c('#006666', '#E5FFFF', '#FFE5CB', '#662700'))
colorLevels<-3
cols_to_use= colors(colorLevels)

## Plot
p_envpos<-ggplot(betapos, aes(x = Var1, y = fct_inorder(name), fill = value)) +
  labs(x = "Fire history", y = "Species", fill = "") +
  geom_tile(color = 'gray60')+
  coord_fixed(ratio = 1/2)+
  scale_fill_steps2(n.breaks = 40,
                    high = cols_to_use[3],
                    mid = 'white',
                    low = cols_to_use[1],
                    nice.breaks = TRUE,
                    na.value = "white",
                    midpoint=0,
                    limits=c(-2, 2)
  ) +
  theme(plot.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = NA),
        legend.margin = margin(l = 1, unit = 'cm'),
        legend.title = element_text(hjust = 0.08, size = 60, face="bold" ),
        legend.key.width = unit(4, 'cm'),
        legend.key.height = unit(10, 'cm'),
        legend.text = element_text(size = 50),
        axis.text.x = element_text(size=60, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_markdown(size=60),
        axis.title.x = element_text(size=60, face="bold"),
        axis.title.y = element_text(size=60, face="bold"),
        plot.margin = margin(t = 0.5,  # Top margin
                             r = 0.1,  # Right margin
                             b = 0.1,  # Bottom margin
                             l = 1,  # Left margin
                             unit = "cm"))+
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0))

p_envneg<-ggplot(betaneg, aes(x = Var1, y = fct_inorder(name), fill = value)) +
  labs(x = "Fire history", y = "Species", fill = "") +
  geom_tile(color = 'gray60')+
  coord_fixed(ratio = 1/2)+
  scale_fill_steps2(n.breaks = 40,
                    high = cols_to_use[3],
                    mid = 'white',
                    low = cols_to_use[1],
                    nice.breaks = TRUE,
                    na.value = "white",
                    midpoint=0,
                    limits=c(-2, 2)
  ) +
  theme(plot.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = NA),
        legend.margin = margin(l = 1, unit = 'cm'),
        legend.title = element_text(hjust = 0.08, size = 60, face="bold" ),
        legend.key.width = unit(4, 'cm'),
        legend.key.height = unit(10, 'cm'),
        legend.text = element_text(size = 50),
        axis.text.x = element_text( size=60, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_markdown(size=60),
        axis.title.x = element_text(size=60, face="bold"),
        axis.title.y = element_text(size=60, face="bold"),
        plot.margin = margin(t = 0.1,  # Top margin
                             r = 1,  # Right margin
                             b = 0.1,  # Bottom margin
                             l = 0.1,  # Left margin
                             unit = "cm"))+
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0))

p_env2<-ggarrange(p_envneg, p_envpos, nrow=1, ncol=2, common.legend = TRUE, legend = "right")

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

## Replace zeros with NA
toPlot_g[toPlot_g == 0] <- NA

## Format the data as matrix and add column and row names
gammaMat = matrix(toPlot_g, nrow = m$nc, ncol = m$nt)
colnames(gammaMat)<- trNames
rownames(gammaMat) <- covNames
gammaMat<- as.data.frame(gammaMat) %>%
  dplyr::slice(-1) #remove intercept

## Update covariate names
rownames(gammaMat) <-c('control', 'spring', 'summer')
colnames(gammaMat)<-c('plants', 'fungi')

## Reformat for ggplot
gammaMatmelt<-as.data.frame(melt(as.matrix(gammaMat)))

## Set the colors
colors = colorRampPalette(c('#006666', '#E5FFFF', '#FFE5CB', '#662700'))
colorLevels<-3
cols_to_use= colors(colorLevels)

## Plot
p_traits<-ggplot(gammaMatmelt, aes(x = Var1, y = Var2, fill = value)) +
  labs(x = "Fire history", y = "Organisms", fill = "") +
  geom_tile(color = 'gray60')+
  geom_text(aes(label = round(value, digits = 3)), color = "gray30", size = 24) +
  coord_fixed()+
  scale_fill_steps2(n.breaks = 40,
                    high = cols_to_use[3],
                    mid = 'white',
                    low = cols_to_use[1],
                    nice.breaks = TRUE,
                    na.value = "white",
                    midpoint=0,
                    limits=c(-2, 2)
  )  +
  theme(plot.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = NA),
        legend.margin = margin(l = 1, unit = 'cm'),
        legend.title = element_text(hjust = 0.08, size = 50, face="bold" ),
        legend.key.width = unit(3, 'cm'),
        legend.key.height = unit(4, 'cm'),
        legend.text = element_text(size = 50),
        axis.text.x = element_text(size=60),
        axis.text.y = element_text(size=60),
        axis.title.x = element_text(size=60, face="bold"),
        axis.title.y = element_text(size=60, face="bold"),
        plot.margin = margin(t = 0.1,  # Top margin
                            r = 0.1,  # Right margin
                            b = 0.1,  # Bottom margin
                            l = 0.1,  # Left margin
                                           unit = "cm"),
        legend.position = "none"
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
significant_species_cor<-colnames(toPlot_sa)
fungi_cor<-significant_species_cor[grepl("ASV", significant_species_cor)]
fungi_cor_taxa<-fungi_taxa %>% filter(ASV %in% fungi_cor) # retrieve taxa of 
# significant fungi
plants_cor<-significant_species_cor[! grepl("ASV", significant_species_cor)]
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
### Fungi
df_fungi_colors_cor<-data.frame(asv= fungi_cor, name= fungi_cor_taxa$name, phylum=fungi_cor_taxa$Phylum )
fcol_code<-c()
for (i in 1:nrow(df_fungi_colors_cor)){
  row_i=df_fungi_colors_cor[i,]
  if (is.na(row_i$phylum)){fcol_code[i]="#FFBBFF"}
  else if (row_i$phylum =="Ascomycota"){fcol_code[i]="#500050"}
  else if (row_i$phylum =="Basidiomycota"){fcol_code[i]="#860086"}
  else if (row_i$phylum == "Chytridiomycota"){fcol_code[i]="#BB00BB"}
  else if (row_i$phylum == "Entorrhizomycota"){fcol_code[i]="#F100F1"}
  else if (row_i$phylum=="Glomeromycota"){fcol_code[i]="#FF50FF"}
  else if (row_i$phylum=="Mortierellomycota"){fcol_code[i]="#FF86FF"}
}
df_fungi_colors_cor$col_code<-fcol_code
### Prepare all names and colors to use on the plot
df_plant_cor<-df_plant_colors_cor[,4:5]
df_fungi_cor<-df_fungi_colors_cor[,c(2,4)]
df_cor<-as.data.frame(rbind(df_plant_cor, df_fungi_cor))
df_cor=df_cor%>% mutate(
  nam_col= glue("<i style='color:{col_code}'>{name}</i>"))
colnames(toPlot_sa)<-df_cor$name
rownames(toPlot_sa)<-df_cor$name
color_cor<-df_cor$col_code
color_cor_use<-c(color_cor[7], color_cor[5], color_cor[1], color_cor[4],
                 color_cor[3], color_cor[2], color_cor[6])

corrplot(toPlot_sa, 
         method = "color", 
         type="upper", 
         order="hclust", 
         tl.col=color_cor_use, # color species names
         col = rev(COL2('RdBu', 20))) # color of the correlation scale
grid.echo()
p_species <- grid.grab() #save the plot

# LEGEND FOR PLOTS

cat<-c("Plant traits", "Forb", "Graminoid","Woody", "* Exotic plant", 
       "Fungal phyla",  "Ascomycota", "Basidiomycota", "Chytridiomycota",
       "Entorrhizomycota", "Glomeromycota","Mortierellomycota", "unknown phylum")
col_code<-c("black", "#00BB00","#008600","#005000", "black", 
            "black","#500050","#860086","#BB00BB",
            "#F100F1","#FF50FF","#FF86FF","#FFBBFF")
df_legend<-data.frame(cat, col_code)
df_legend$x<-rep(1, 13)
df_legend$y<-c(30, 27, 25, 23, 21, 17, 14, 12, 10, 8, 6, 4, 2)
p_text<-ggplot(df_legend, aes(x, y)) + 
  geom_text(aes(label=cat), color=col_code, size=c(14,12,12,12,12,14,12,12,12,12,12,12,12))+
  xlim(c(0.98,1.02))+
    theme_void()
