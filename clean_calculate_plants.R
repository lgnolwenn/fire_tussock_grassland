# PLANT DATA FROM DEEP STREAM
## Code for cleaning, calculating percentages of cover, retrieving traits and
## getting changes over time of plant data from Deep Stream grassland.
## Note: paths need to be changed to match your device.

# LOAD PACKAGES
library(tidyverse)

# LOAD DATA
plant<-read.csv("C:/Users/nolwe/Documents/M2/internship/project/data/plant.csv")

# CLEANING DATA

## Removing unnecessary data for our analyses
### exclosure subplots are not interesting for our analyses
plant_data_noE<-subset(plant, !grepl("_E_",plant$subplot_id))
### remove two extra unburnt subplots
plant_data_noE_noub<-subset(plant_data_noE,
                            !grepl("UB_1_O_4",plant_data_noE$subplot_id)
                            & !grepl("UB_1_O_5",plant_data_noE$subplot_id))

plant_data<-plant_data_noE_noub


# CALCULATING A PERCENTAGE OF COVER OF EACH PLANT

# For each subplot, and for each species, we calculated the percentage of cover.
# For this, we took the mid-point of the covers score range in each quadrat and 
# calculated the mean of the four quadrats of each subplot.

csp<-c()
cs=plant_data$cover_score
for(i in 1:nrow(plant_data)){
  if(cs[i]==1){csp[i]=0.5} #When cover score=1, cover is below 1%
  else if(cs[i]==2){csp[i]=3.5} #When cover score=2, cover is between 2 and 5%, mid-point being 3.5%
  else if(cs[i]==3){csp[i]=14.5}#When cover score=3, cover is between 6 and 25%, mid-point being 14.5%
  else if(cs[i]==4){csp[i]=38}#When cover score=4, cover is between 26 and 50%, mid-point being 38%
  else if(cs[i]==5){csp[i]=63}#When cover score=5, cover is between 51 and 75%, mid-point being 63%
  else if(cs[i]==6){csp[i]=88}#When cover score=2, cover is between 76 and 100%, mid-point being 88%
  else{csp[i]=NA} #in case one is not well marked
}
plant_data$cover_score_percent<-csp #add the percentage in a new column of the data frame

## Getting the mean
mean_csp=plant_data %>%
  group_by(months.after.fire, subplot_id, PreferredName) %>% #for the same subplot and the same species
  summarise(cover_score_percent=mean(cover_score_percent)) #calculate the mean cover of the four quadrats

# CREATING A DATA FRAME WITH MEAN COVER OF QUADRATS FOR EACH SUBPLOT

plant_species<-unique(mean_csp$PreferredName)
sampling_units<-unique(mean_csp$subplot_id)
months<-unique(mean_csp$months.after.fire)
#create time/subplot variable
time_sampling_units<-paste0(rep(months, each=length(sampling_units)),"_",sampling_units)
time_sampling_units_correct<-time_sampling_units[!grepl("13_UB", time_sampling_units)&!grepl("26_UB", time_sampling_units)]
#new data frame
plant_df<-data.frame(matrix(0, nrow = length(time_sampling_units), ncol = length(plant_species)))
s=0
for(t in months){
  tcover=mean_csp%>% filter(months.after.fire==t)
  r=s
  for (su in sampling_units){
    r=r+1
    sutcover=tcover%>% filter(subplot_id==su)
    c=0
    for (ps in plant_species){
      c=c+1
      v=sutcover%>% filter(PreferredName==ps)
      if(dim(v)[1]==0){plant_df[r,c]=0}
      else{plant_df[r,c]=v$cover_score_percent}
    }
  }
  s=r
}
rownames(plant_df)<-time_sampling_units
colnames(plant_df)<-plant_species
plant_df<-plant_df[time_sampling_units_correct,]

# KEEPING ONLY FIRST SAMPLING, 2 MONTHS POST-FIRE

plant_df_t1=subset(plant_df,!grepl("13_",rownames(plant_df))&!grepl("26_",rownames(plant_df)))
plant_df_t1=plant_df_t1%>%select(-which(colSums(plant_df_t1)==0))
plant_species_t1=unique(colnames(plant_df_t1))

# PLANT TRAITS

plant_info<-plant
### Retrieve the traits of plants
traits=plant_info %>%
  group_by(PreferredName)%>%
  select(BioStatus, GrowthForm)
traits=unique(traits)
###Simplify traits
####Bio Status
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

# CHANGES
plants_temp<-plant_df
## Rearranging the matrix
plants_temp$subplot_id<-time_sampling_units_correct
plants_temp<-plants_temp[, c(45, 1:44)]
colnames(plants_temp)[1]="subplot_id"
plants_temp=arrange(plants_temp, desc(subplot_id))
## Removing unburnt
ub_index=grepl("UB", plants_temp$subplot_id) #get indexes of UB
plants_temp_nub=plants_temp[!ub_index,] #remove the UB
## Removing species only present in unburnt plots
df_sp_rm=as.data.frame(colSums(plants_temp_nub[,-1])) #calculate the sum of cover for each species
colnames(df_sp_rm)="colsum"
plants_temp_nub_sp=plants_temp_nub[,-1] #remove informative, but non-numeric column
plants_temp_nub_sp=plants_temp_nub_sp[,df_sp_rm$colsum!=0] #if the sum is 0 the species was only in unburnt plots
plants_temp=cbind(plants_temp_nub[,1],plants_temp_nub_sp)#new data frame without species that are not there
## Rearranging the matrix again
colnames(plants_temp)[1]="subplot_id"
plants_temp=arrange(plants_temp, desc(subplot_id))
## Getting each time
### time 1
time1<- plants_temp[grepl("2_", plants_temp$subplot_id),]
time1<-time1[!grepl("13_", time1$subplot_id),]
time1<-time1[!grepl("26_", time1$subplot_id),]
### time 2
time2<-plants_temp[grepl("13_", plants_temp$subplot_id),]
### time 3
time3<-plants_temp[grepl("26_", plants_temp$subplot_id),]
## Getting the subplots names
sampling_units<- plants_temp[28:54,"subplot_id"]
sampling_units<-gsub("26_","", sampling_units)
## Constructing change data frames
change_1_2<-time2[,-1]-time1[,-1] #changes between time 1 (2 months) and time 2 (13 months)
rownames(change_1_2)<-sampling_units
change_1_3<-time3[,-1]-time1[,-1] #changes between time 2 (13 months) and time 3 (26 months)
rownames(change_1_3)<-sampling_units
change_2_3<-time3[,-1]-time2[,-1] #changes between time 1 (2 months) and time 3 (26 months)
rownames(change_2_3)<-sampling_units

# SAVING
## Percentage of plant cover all times
write.csv(plant_df, file="C:/Users/nolwe/Documents/M2/internship/project/data/temporal_plant_data.csv")
## Percentage of plant cover t1 = 2 months post-fire
write.csv(plant_df_t1, file="C:/Users/nolwe/Documents/M2/internship/project/data/plant_data.csv")
## Traits data
write.csv(traits, file="C:/Users/nolwe/Documents/M2/internship/project/data/plant_traits.csv")
## Changes data frame
write.csv(change_1_2, "C:/Users/nolwe/Documents/M2/internship/project/data/change_1_2.csv")
write.csv(change_1_3, "C:/Users/nolwe/Documents/M2/internship/project/data/change_1_3.csv")
write.csv(change_2_3, "C:/Users/nolwe/Documents/M2/internship/project/data/change_2_3.csv")