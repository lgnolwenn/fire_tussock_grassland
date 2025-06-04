# FUNGAL DATA FROM DEEP STREAM
## Code for cleaning and selecting the fungal Amplicon Sequence Variants (ASV)
## from Deep Stream grassland. Once selected, retrieve the taxonomic identification 
## of each ASV of intertussock fungi (those used in the model).
## Note: paths need to be changed to match your device.

# LOAD PACKAGES
library(tidyverse)
library(GUniFrac)

# LOAD DATA
ASV_original <- read.csv("C:/Users/nolwe/Documents/M2/internship/project/data/ASV_table_raw.csv")
fungi_taxa<-read.csv("C:/Users/nolwe/Documents/M2/internship/project/data/ftaxa.csv")

# RAREFACTION

## renaming columns
col_names=colnames(ASV_original)
col_names=col_names[-length(col_names)]
colnames(ASV_original)=c("samples",col_names)

## Inspecting number of reads per sample. Read data from csv file and group by samples
read_count <- ASV_original 
read_count <- pivot_longer(read_count, !samples, names_to = "ASV", values_to = "Count")
read_count <- group_by(read_count,samples)
num_reads <- read_count %>% 
  summarise(total_reads = sum(Count))
reads_vector <- num_reads$total_reads

## Removing samples with <15,000 reads
ASV_filt <- ASV_original %>%
  pivot_longer(., !samples, names_to = "ASV", values_to = "Count") %>%
  group_by(samples) %>%
  mutate(Reads = sum(Count)) %>%
  pivot_wider(.,names_from = ASV, values_from = Count) %>%
  filter(Reads > 15000)

## rarefaction to 15,000 reads
ASV_to_rarefy <- ASV_filt
ASV_to_rarefy <- subset(ASV_to_rarefy, select=c(-Reads, -samples))
is.data.frame(ASV_to_rarefy)
ASV_rarefied <- GUniFrac::Rarefy(ASV_to_rarefy, 15000)
ASV_rarefied_df <- do.call(rbind.data.frame, ASV_rarefied)

# REORGANIZING THE DATASET

## Renaming with a shorter name the dataframe
fdr=ASV_rarefied_df

## Add informativ columns
### Subplot
fdr$samples<-ASV_filt$samples 
fdr$subplot_id<- substring(fdr$samples, 5) #Retrieve the subplot identification
# rename subplots: no "O" because E plots are not used so there is no need to
# differenciate them; give unique number to unburnt (UB) plots.
fdr$subplot_id<-gsub("-O-", "_", fdr$subplot_id)
fdr$subplot_id<-gsub("-UB1-", "10_", fdr$subplot_id)
fdr$subplot_id<-gsub("-UB2-", "11_", fdr$subplot_id)
fdr$subplot_id<-gsub("-UB3-", "12_", fdr$subplot_id)
### Location
fdr$location<-substr(fdr$samples, 4, 4) #retrieve the location of sample, either 
#inter-tussock (A) or Chinochloa adjacent (S)
### Treatment
#CO for control not burnt in 2001 but in 2019, SU for summer burnt in 2001 and 
#burnt in 2019, SP for spring burnt in 2001 and burnt in 2019, UB for unburnt
fdr$treatment<-c()
for (i in 1:length(fdr$subplot_id)){
  sbp=fdr$subplot_id
  sbp_i=sbp[i]
  if (substr(sbp_i, 1, 2) %in% c(10,11,12)){fdr$treatment[i]="UB"}
  else if (substr(sbp_i, 1, 1) %in% c(1,6,7)){fdr$treatment[i]="CO"}
  else if (substr(sbp_i, 1, 1) %in% c(3,4,8)){fdr$treatment[i]="SP"}
  else if (substr(sbp_i, 1, 1) %in% c(2,5,9)){fdr$treatment[i]="SU"}
  else{fdr$treatment[i]=NA}
}
### Reorder
lfdr=length(fdr)
fdr<-fdr[, c(lfdr, lfdr-1, lfdr-2, lfdr-3, 1:(lfdr-4))] #Put identification columns first
fdr<-as.data.frame(fdr)

# CREATING TWO DATA FILES ONE FOR EACH LOCATION (Intertussock/Chionochloa-adjacent)

fdr_a<-fdr %>% filter(fdr$location=="A") #intertussock dataframe
fdr_s<-fdr %>% filter(fdr$location=="S") #chionochloa-adjacent dataframe

## Removing data not relevant anymore
###remove the useless location columns
fdr_a<-fdr_a%>% select(-location)
fdr_s<-fdr_s%>% select(-location)
### Removing ASV that were discarded by rarefaction
#### Intertussock
lfdra=length(fdr_a)
fdr_a_asv=fdr_a[,4:(lfdra-1)]
number_of_obs_asv=colSums(fdr_a_asv)
empty_asv=which(number_of_obs_asv==0) #ASV with no reads
fdr_a_asv=fdr_a_asv%>%select(-empty_asv)
#### Chionochloa-adjacent
lfdrs=length(fdr_s)
fdr_s_asv=fdr_s[,4:(lfdrs-1)]
number_of_obs_asv=colSums(fdr_s_asv)
empty_asv=which(number_of_obs_asv==0) #ASV with no reads
fdr_s_asv=fdr_s_asv%>%select(-empty_asv)

## Going from number of reads to proportions
read=15000 #number of reads
### Intertussock
fungi_prop_a <- fdr_a_asv %>% 
  mutate(across(everything(), ~ (as.numeric(.)) / read))
fdr_a_prop<-cbind(fdr_a[,c(1,2)], fungi_prop_a)
### Chionochloa-adjacent
fungi_prop_s <- fdr_s_asv %>% 
  mutate(across(everything(), ~ (as.numeric(.)) / read))
fdr_s_prop<-cbind(fdr_s[,c(1,2)], fungi_prop_s)

# FILTERING
fungia<-fdr_a_prop
fungis<-fdr_s_prop
### Intertussock
asv_prop_a<-colnames(fungia)[-c(1:3)]
fungi_prop_asv_a = fungia %>%
  select(asv_prop_a)
### Chionochloa-adjacent
asv_prop_s<-colnames(fungis)[-c(1:3)]
fungi_prop_asv_s = fungis %>%
  select(asv_prop_s)

## SAMPLES CONDITIONS
### Intertussock
#Number of samples in which the ASV is present
spl_without_asv_a<-34-colSums(fungi_prop_asv_a==0)
#Get the ASV to be present in 5 or more samples
names_s1_a=names(spl_without_asv_a[spl_without_asv_a>4])
# Get the ASv to be present in 32 or less samples
names_s2_a=names(spl_without_asv_a[spl_without_asv_a>31])
# Combine both
names_s_a=names_s1_a[!names_s1_a %in% names_s2_a]
### Chionochloa-adjacent
#Number of samples in which the ASv is present
spl_without_asv_s<-34-colSums(fungi_prop_asv_s==0)
#Get the ASV to be present in 5 or more samples
names_s1_s=names(spl_without_asv_s[spl_without_asv_s>4])
# Get the ASv to be present in 32 or less samples
names_s2_s=names(spl_without_asv_s[spl_without_asv_s>31])
# Combine both
names_s_s=names_s1_s[!names_s1_s %in% names_s2_s]

## READS CONDITIONS
###Intertussock
#The ASV scores at least 1% of the sample's reads for at least one sample
prop=150/15000
reads_names<-c()
for (j in 1:nrow(fungi_prop_asv_a)){
  for (k in 1:ncol(fungi_prop_asv_a)){
    asvjk=fungi_prop_asv_a[j,k]
    if(asvjk>prop){
      reads_names<-c(reads_names, colnames(fungi_prop_asv_a)[k])
    }
  }
}
names_r1_a<-unique(reads_names)
# The ASV scores at least 0.1% of total reads across all samples
prop2=15*34/15000
asv_sum_a=colSums(fungi_prop_asv_a)
names_r2_a=names(asv_sum_a[asv_sum_a>prop2])
# Combine both reads'conditions
names_r_a= names_r1_a[names_r1_a %in% names_r2_a]
### Chionochloa-adjacent
#The ASV scores at least 1% of the sample's reads for at least one sample
prop=150/15000
reads_names<-c()
for (j in 1:nrow(fungi_prop_asv_s)){
  for (k in 1:ncol(fungi_prop_asv_s)){
    asvjk=fungi_prop_asv_s[j,k]
    if(asvjk>prop){
      reads_names<-c(reads_names, colnames(fungi_prop_asv_s)[k])
    }
  }
}
names_r1_s<-unique(reads_names)
# The ASV scores at least 0.1% of total reads across all samples
prop2=15*34/15000
asv_sum_s=colSums(fungi_prop_asv_s)
names_r2_s=names(asv_sum_s[asv_sum_s>prop2])
# Combine both reads'conditions
names_r_s= names_r1_s[names_r1_s %in% names_r2_s]

## SAMPLES AND READS CONDITIONS
names_sr_a=names_s_a[names_s_a %in% names_r_a]
names_sr_s=names_s_s[names_s_s %in% names_r_s]

# SELECTING
fungi_a<-fungia %>% select (1, 2, names_sr_a)
fungi_s<-fungis %>% select (1, 2, names_sr_s)

# ATTRIBUTING TAXA TO ASV FOR INTERTUSSOCK
fungi_data<-fungi_a
##Curation of data
fungi_taxa<-fungi_taxa[, c(1:7, 15)]
asv_sel<-colnames(fungi_data)
asv_sel<-asv_sel[-(1:3)]
fungi_taxa_sel<-fungi_taxa %>% filter (ASV %in% asv_sel)
## Last level of phylogeny known
fungi_taxa_sel$name<-paste0(fungi_taxa_sel$Genus, " ", fungi_taxa_sel$Species)
for (i in 1:nrow(fungi_taxa_sel)){
  df=fungi_taxa_sel[i,]
  if(df$name=="NA NA"){
    if (!is.na(df$Family)){df$name= paste0(df$Family, " NA") }
    else if (!is.na(df$Order)){df$name= paste0(df$Order, " NA") }
    else if (!is.na(df$Class)){df$name= paste0(df$Class, " NA") }
    else if (!is.na(df$Phylum)){df$name= paste0(df$Phylum, " NA") }
    else {df$name= "Fungi NA"}
  }
  fungi_taxa_sel[i, "name"]=df$name
}
for (i in 1:nrow(fungi_taxa_sel)){
  df=fungi_taxa_sel[i,]
  if(grepl("NA", df$name)){
    num=substr(df$ASV, 5, 8)
    num_name=paste0("_", num)
    df$name<-gsub(" NA", num_name, df$name)
  }
  else{
    num=substr(df$ASV, 5, 8)
    df$name<-paste0(df$name, "_", num)
  }
  fungi_taxa_sel[i, "name"]=df$name
}

# SAVING
#all asv
write.csv(fdr_a_prop, file="C:/Users/nolwe/Documents/M2/internship/project/data/fungi_inter_tussock_data_prop.csv")
write.csv(fdr_s_prop, file="C:/Users/nolwe/Documents/M2/internship/project/data/fungi_chionochloa_data_prop.csv")
#selected asv
write.csv(fungi_a, file="C:/Users/nolwe/Documents/M2/internship/project/data/fungi_inter_tussock_final_selection.csv")
write.csv(fungi_s, file="C:/Users/nolwe/Documents/M2/internship/project/data/fungi_chionochloa_final_selection.csv")
# taxa of selected asv from intertussock
write.csv(fungi_taxa_sel, file="C:/Users/nolwe/Documents/M2/internship/project/data/selected_fungi_taxa.csv")
