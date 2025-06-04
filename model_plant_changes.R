# HMSC MODEL OF THE EFFECTS OF FIRE HISTORY ON CHANGES OF THE PLANT COMMUNITY OF DEEP STREAM
## Code for defining and running the HMSC models of changes of plant community.
## Note: paths need to be changed to match your device.

# LOAD PACKAGES

library(tidyverse)
library(Hmsc)

# LOAD DATA

change_1_2<-read.csv("/nfs/home/leguyano/change_1_2.csv")
change_1_3<-read.csv( "/nfs/home/leguyano/change_1_3.csv")
change_2_3<-read.csv( "/nfs/home/leguyano/change_2_3.csv")
plant_info<-read.csv("/nfs/home/leguyano/plant.csv")

# MATRIXES

## Y, community matrix

### Y_1_2 changes between 1st (2 months post-fire) and 2nd sampling time (13 months post-fire)
#### Rearrange the matrix
colnames(change_1_2)[1]="subplot_id"
change_1_2=arrange(change_1_2, subplot_id)
subplot_id=change_1_2$subplot_id
# #### Make sure there are no empty columns
change_1_2=Filter(function(x) length(unique(x))>1, change_1_2)
#### Transform into percentage and taxe the Z score for standardization
change_1_2_abundance=change_1_2[,2:length(change_1_2)]
change_1_2_abundance=change_1_2_abundance/100
change_1_2_abundance_mean_standard <- change_1_2_abundance %>% 
  mutate(across(everything(), ~ (. - mean(.)) / sd(.)))
#### Define the Y matrix and make sure that row and column names are correct
Y_1_2=as.matrix(change_1_2_abundance_mean_standard)
rownames(Y_1_2)=change_1_2$subplot_id

### Y_2_3 changes between 2nd (13 months post-fire) and 3rd sampling time (26 months post-fire)
#### Rearrange the matrix 
colnames(change_2_3)[1]="subplot_id"
change_2_3=arrange(change_2_3, subplot_id)
# #### Make sure there are no empty columns
change_2_3=Filter(function(x) length(unique(x))>1, change_2_3)
#### Transform into percentage and taxe the Z score for standardization
change_2_3_abundance=change_2_3[,2:length(change_2_3)]
change_2_3_abundance=change_2_3_abundance/100
change_2_3_abundance_mean_standard <- change_2_3_abundance %>% 
  mutate(across(everything(), ~ (. - mean(.)) / sd(.)))
#### Define the Y matrix and make sure that row and column names are correct
Y_2_3=as.matrix(change_2_3_abundance_mean_standard)
rownames(Y_2_3)=change_2_3$subplot_id

### Y_1_3 changes between 1st (2 months post-fire) and 3rd sampling time (26 months post-fire)
#### Rearrange the matrix
colnames(change_1_3)[1]="subplot_id"
change_1_3=arrange(change_1_3, subplot_id)
# #### Make sure there are no empty columns
change_1_3=Filter(function(x) length(unique(x))>1, change_1_3)
#### Transform into percentage and taxe the Z score for standardization
change_1_3_abundance=change_1_3[,2:length(change_1_3)]
change_1_3_abundance=change_1_3_abundance/100
change_1_3_abundance_mean_standard <- change_1_3_abundance %>% 
  mutate(across(everything(), ~ (. - mean(.)) / sd(.)))
#### Define the Y matrix and make sure that row and column names are correct
Y_1_3=as.matrix(change_1_3_abundance_mean_standard)
rownames(Y_1_3)=change_1_3$subplot_id

## X, environmental matrix
treatment<-c(rep("CO", 9), rep("SP", 9), rep("SU", 9))
XData=data.frame(as.factor(treatment)
                 #, as.factor(time)
)
rownames(XData)=rownames(Y_1_2) #make sure rownames of X match the rownames (i.e. sampling units) of all Y
colnames(XData)<-c("treatment"
                   #, "time"
)
XData$treatment <- relevel(XData$treatment, ref="CO") # set the "control" as the reference

## Tr, traits matrix

### Retrieve the traits of plants
traits=plant_info %>%
  group_by(PreferredName)%>%
  select(BioStatus, GrowthForm)
traits=unique(traits)
###Simplify traits
#### Bio Status
for(i in 1:nrow(traits)){
  if(grepl("Indigenous",traits[i,"BioStatus"])==TRUE)
  {traits[i,"BioStatus"]<-"Indigenous"}
}
#### Growth Form
for(i in 1:nrow(traits)){
  if(grepl("Woody",traits[i,"GrowthForm"])==TRUE)
  {traits[i,"GrowthForm"]<-"Woody"}
  else if(grepl("Shrub",traits[i,"GrowthForm"])==TRUE)
  {traits[i,"GrowthForm"]<-"Woody"}
  else if (traits[i,"GrowthForm"]=="Fern")
  {traits[i,"GrowthForm"]<-"Forb"}
}
### Tr_1_2
species_1_2=colnames(Y_1_2)
species_1_2<-gsub(".", " ", species_1_2, fixed=TRUE)
traits_1_2=traits[match(species_1_2, traits$PreferredName), ]   
### Build the data frame
Tr_1_2=data.frame(as.factor(traits_1_2$BioStatus), as.factor(traits_1_2$GrowthForm))
colnames(Tr_1_2)=c("BioStatus", "GrowthForm")
rownames(Tr_1_2)=colnames(Y_1_2) #make sure rownames of Tr_1_2 match the colnames (i.e. species) of Y_1_2
Tr_1_2$BioStatus <- relevel(Tr_1_2$BioStatus, ref="Indigenous")
Tr_1_2$GrowthForm <- relevel(Tr_1_2$GrowthForm, ref="Graminoid")

### Tr_2_3
species_2_3=colnames(Y_2_3)
species_2_3<-gsub(".", " ", species_2_3, fixed=TRUE)
traits_2_3=traits[match(species_2_3, traits$PreferredName), ]   
### Build the data frame
Tr_2_3=data.frame(as.factor(traits_2_3$BioStatus), as.factor(traits_2_3$GrowthForm))
colnames(Tr_2_3)=c("BioStatus", "GrowthForm")
rownames(Tr_2_3)=colnames(Y_2_3) #make sure rownames of Tr_2_3 match the colnames (i.e. species) of Y_2_3
Tr_2_3$BioStatus <- relevel(Tr_2_3$BioStatus, ref="Indigenous")
Tr_2_3$GrowthForm <- relevel(Tr_2_3$GrowthForm, ref="Graminoid")

### Tr_1_3
species_1_3=colnames(Y_1_3)
species_1_3<-gsub(".", " ", species_1_3, fixed=TRUE)
traits_1_3=traits[match(species_1_3, traits$PreferredName), ]   
### Build the data frame
Tr_1_3=data.frame(as.factor(traits_1_3$BioStatus), as.factor(traits_1_3$GrowthForm))
colnames(Tr_1_3)=c("BioStatus", "GrowthForm")
rownames(Tr_1_3)=colnames(Y_1_3) #make sure rownames of Tr_1_3 match the colnames (i.e. species) of Y_1_3
Tr_1_3$BioStatus <- relevel(Tr_1_3$BioStatus, ref="Indigenous")
Tr_1_3$GrowthForm <- relevel(Tr_1_3$GrowthForm, ref="Graminoid")

##Formulae
XFormula = ~treatment
TrFormula= ~BioStatus+GrowthForm #no interactions between traits

## Study design
### Plot  and sample levels
plot.id<-c()
sample.id<-c()
sbp=change_1_2$subplot_id
for (i in 1:length(change_1_2$subplot_id)){
  sbp_i=sbp[i]
  #plot
  plot.id[i]=as.numeric(substr(sbp_i, 4, 4))
  #sample
  s=substr(sbp_i,8,9)
  sample.id[i]=paste0(plot.id[i], "_", s)
} 
plot.id<-as.factor(plot.id)
sample.id<-as.factor(sample.id)

### Study design and random levels
studyDesign = data.frame(sample = sample.id, plot = plot.id)
rL = HmscRandomLevel(units = studyDesign$plot) #plot is the unique random level here

# MODELS

temporal_models_plants = list()
for (i in 1:3){
  Y = switch(i, Y_1_2, Y_2_3, Y_1_3)
  TrData = switch(i, Tr_1_2, Tr_2_3, Tr_1_3)
  m = Hmsc(Y = Y, XData = XData, XFormula = XFormula, TrData = TrData, 
           TrFormula=TrFormula, studyDesign = studyDesign, 
           ranLevels = list("plot"=rL),distr = "normal")
  temporal_models_plants[[i]] = m
}

thin = 100 # thinning of sampling
samples = 250 # number of posterior samples per chain
transient = ceiling(0.5*samples*thin) # number of first iterations used as burn in after which sampling starts
nChains = 4 # number of MCMC chains sampled

for (i in 1:3){
  temporal_models_plants[[i]] = sampleMcmc(temporal_models_plants[[i]], thin = thin,
                           samples = samples,transient = transient,
                           nChains = nChains)
}

save(temporal_models_plants, file="/nfs/home/leguyano/models_plants_changes_3.RData")
