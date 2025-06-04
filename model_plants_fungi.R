# HMSC MODEL OF THE EFFECTS OF FIRE HISTORY ON PLANT AND SOIL FUNGAL COMMUNITIES OF DEEP STREAM
## Code for defining and running the HMSC model of both plant and fungal communities.
## Note: paths need to be changed to match your device.

# LOAD PACKAGES

library(tidyverse)
library(Hmsc)

# LOAD DATA
fungi_data<-read.csv("/nfs/home/leguyano/fungi_inter_tussock_final_selection.csv")
plants<-read.csv("/nfs/home/leguyano/plant_data.csv")

# MATRIX

## Y, community matrix

### Preparing the data
#### Fungi
fungi_data<-fungi_data[,-1]
fungi_data$X<-paste0(fungi_data$treatment, "_", fungi_data$subplot_id)
fungi_data=arrange(fungi_data, desc(fungi_data$X)) # in order to get the same order for both fungi and plants
fungi<-fungi_data
fungi<-fungi[,-c(1:2, length(fungi_data))] #removing unnecessary columns
#### Plants
colnames(plants)[1]="subplot_id"
plants=arrange(plants, desc(subplot_id)) # in order to get the same order for both fungi and plants
##### Renaming subplots to match the name of subplots in fungal dataset
plants$subplot_id<-gsub("_O_", "_", plants$subplot_id)
x<-c()
y<-c()
sbp=plants$subplot_id
for (i in 1:length(plants$subplot_id)){
  sbp_i=sbp[i]
  t=substr(sbp_i, 3, 4)
  if(t=="UB"){x[i]=as.numeric(substr(sbp_i, 6, 6))+9}
  else{x[i]=as.numeric(substr(sbp_i, 6, 6))}
  s=substr(sbp_i,8,9)
  y[i]=paste0(t, "_", x[i], "_", s)
} 
plants$subplot_id<-y
plants<-plants %>% filter (plants$subplot_id %in% fungi_data$X) #remove subplots 
#that were removed by the processing of fungal sequences
plant_abundance=plants[,2:length(plants)] #keep only useful data
plant_abundance=plant_abundance/100 # transform values into percentages

### Combined and Z-score
pf_com<-cbind(plant_abundance, fungi)
abundance_mean_standard <- pf_com %>% 
  mutate(across(everything(), ~ (. - mean(.)) / sd(.))) #Z-score to standardize the data

### Define the Y matrix and make sure that row and column names are correct
Y=as.matrix(abundance_mean_standard)
rownames(Y)=fungi_data$X

## X, environmental matrix
treatment<-fungi_data$treatment
XData=data.frame(as.factor(treatment)) # categorical data need to be treated as a factor
rownames(XData)=rownames(Y) #make sure rownames of X match the rownames (i.e. sampling units) of Y
colnames(XData)<-c("treatment")
XData$treatment <- relevel(XData$treatment, ref="UB") #set "unburnt" as the 
# reference instead of alphabetic order

## Tr, traits matrix
organism<-c(rep("plant", length(plant_abundance)), rep("fungi", length(fungi)))
TrData<- data.frame(as.factor(organism)) # categorical data need to be treated as a factor
colnames(TrData)<-c("organism")
rownames(TrData)<-colnames(Y) #make sure rownames of Tr match the colnames (i.e. species) of Y
TrData$organism<-relevel(TrData$organism, ref="plant") #set "plant" as the 
# reference instead of alphabetic order

## Formulae
XFormula = ~treatment
TrFormula= ~organism

## Study design
### Plot level
plot.id<-substr(fungi_data$subplot_id, 1, 2)
plot.id<-gsub("_", "", plot.id)
plot.id<-as.factor(plot.id)
### Sampling units level = subplot
sample.id<-paste0(fungi_data$location, fungi_data$subplot_id)
sample.id<-as.factor(sample.id)
### Study design and random levels
studyDesign = data.frame(sample = sample.id, plot = plot.id)
rL = HmscRandomLevel(units = studyDesign$plot) #plot is the unique random level here


# MODEL

the_model = Hmsc(Y = Y, XData = XData, XFormula = XFormula, TrData = TrData, TrFormula = TrFormula,
         studyDesign = studyDesign, ranLevels =list("plot"=rL), distr = "normal")

## Running MCMC
thin = 100 # thinning of sampling
samples = 250 # number of posterior samples per chain
transient = ceiling(0.5*samples*thin) # number of first iterations used as burn in after which sampling starts
nChains = 4 # number of MCMC chains sampled

the_model = sampleMcmc(the_model, thin = thin,
               samples = samples, transient = transient,
               nChains = nChains)

save(the_model, file="/nfs/home/leguyano/model_all_2.Rdata")