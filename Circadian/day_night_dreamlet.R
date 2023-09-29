.libPaths(c("/sc/arion/projects/psychAD/Madeline/Rlib_4_3","/sc/arion/projects/CommonMind/kiran/RLib_4_3",.libPaths()))
.libPaths(c("/sc/arion/projects/CommonMind/kiran/RLib_4_3",.libPaths()))
library(zellkonverter)
library(basilisk)
library(dreamlet)
library(crumblr)
library(foreach)
library(doParallel)
suppressPackageStartupMessages({
library(SingleCellExperiment)
library(dreamlet)
library(zenith)
library(DelayedArray)
library(GSEABase)
library(scater)
library(tidyverse)
library(ggplot2)
library(cowplot) 
library(org.Hs.eg.db)
library(R.utils)
})
library(data.table)

args = commandArgs(trailingOnly=TRUE)

cell_groups=args[1] #subclass or subtype
celltype=args[2] #Astro or OPCs. 
age_group=args[3] # adulthood or adulthood or old


file="/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/h5ad_final/AGING_2023-04-01_16_14.h5ad"
sce = readH5AD(file, use_hdf5=TRUE, verbose=TRUE)
assayNames(sce)[1] = "counts"
sce1=subset(sce, ,subclass==celltype)

pb = aggregateToPseudoBulk(sce1,
      assay = "counts",
      cluster_id = "subclass",
      sample_id = "SubID")

colData(pb)$SubID=rownames(colData(pb))


### creating age groups

df=colData(pb)[,c("Age","Sex","Source","SubID")]
df$groups="Postnatal"
df$groups[df$Age>=1 & df$Age<12]="early_childhood"
df$groups[df$Age>=12 & df$Age<20]="late_childhood"
df$groups[df$Age>=20 & df$Age<40]="early_adulthood"
df$groups[df$Age>=40 & df$Age<60]="late_adulthood"
df$groups[df$Age>=60]="old"

#1.1 DxN of Individual Groups 
#Group 4_Early Adulthood(20-39):   Day - n = 18    Night - n = 24
#Group 5_Late Adulthood(40-59):    Day - n = 49    Night - n = 16
#Group 6_Old Age (>=60):            Day - n = 51    Night - n = 26

#1.2 DxN of Combined Groups 
#Group 4-5_Adulthood(20-59):       Day - n = 67    Night - n = 40
#Group 4-6_Adulthood+Old Age(>=20)  Day - n = 118   Night - n = 66

contrast_age_groups_list=list(
"early_adulthood"=df$SubID[df$Age>=20 & df$Age<40],
"late_adulthood"=df$SubID[df$Age>=40 & df$Age<60],
"old" = df$SubID[df$Age>=60],
"adulthood"=df$SubID[df$Age>=20 & df$Age<60],
"adulthood_old"=df$SubID[df$Age>=20])

keep_subids_age_group=contrast_age_groups_list[[age_group]]

### print number of samples in each group

lapply(contrast_age_groups_list,length)

### preparing data for covariate correction

df = as.data.frame(colData(pb))
df$SubID=rownames(df)

metadata_aggrs=as.data.frame(metadata(pb)$aggr_means)
metadata_df=merge(metadata_aggrs,df,by="SubID")
rownames(metadata_df)=metadata_df$SubID
metadata_df_flt=metadata_df[match(colnames(pb),metadata_df$SubID),]
identical(colnames(pb),as.character(metadata_df_flt$SubID))
metadata_df_flt$log_n_counts=log(metadata_df_flt$n_counts)
metadata_df_flt$scale_PMI=scale(metadata_df_flt$PMI)

### adding day night data for dreamlet analysis

metadata_daynight=read.csv("/sc/arion/projects/psychAD/aging/sleep_patterns/metadata_aging_TOD_DayNight.csv")
metadata_daynight_merged=merge(metadata_daynight,metadata_df_flt,by="SubID")
metadata_daynight_merged=metadata_daynight_merged[metadata_daynight_merged$Day_Night!=0,]

### subset samples for a given age group contrast

metadata_daynight_merged_subset=metadata_daynight_merged[metadata_daynight_merged$SubID %in% keep_subids_age_group,]

### subsetting pseudobulk for a given age group 

keep_subids=colnames(pb)[match(metadata_daynight_merged_subset$SubID,colnames(pb))]
pb_subset=subset(pb, ,SubID %in% keep_subids)
metadata_daynight_merged_ordered = metadata_daynight_merged_subset[match(colnames(pb_subset),as.character(metadata_daynight_merged_subset$SubID)),]
identical(as.character(metadata_daynight_merged_ordered$SubID),colnames(pb_subset))
rownames(metadata_daynight_merged_ordered)=metadata_daynight_merged_ordered$SubID
metadata_daynight_merged_ordered$cat_Day_Night=factor(metadata_daynight_merged_ordered$Day_Night,levels=c(1,2),labels=c("Day","Night"))
colData(pb_subset)=DataFrame(metadata_daynight_merged_ordered)
colData(pb_subset)$scale_age=scale(colData(pb_subset)$Age)

### making the day night contrast as factor
table(metadata_daynight_merged_ordered$cat_Day_Night)
colData(pb_subset)$scale_age=scale(colData(pb_subset)$Age) #[MRS note] think this is an error - same as line 107

#[MRS] I think this is what we need to add to add day night contrast as a factor
colData(pb_subset)$cat_Day_NightDay = factor(pb_subset$cat_Day_Night, level=c("Day"))
colData(pb_subset)$cat_Day_NightNight = factor(pb_subset$cat_Day_Night, level=c("Night"))

## to see the impact of Day night on gene expression

model1 <- as.formula(" ~ 0  + cat_Day_Night + (1 | Source) + (1 | prep) + (1 | pool) + scale(PMI) + log_n_counts + (1|Sex)")
contrasts <- c(Diff = "cat_Day_NightDay - cat_Day_NightNight")
res.proc_model1=processAssays(pb_subset,model1, min.count=5)
res.dl_model1=dreamlet(res.proc_model1,model1,contrasts = contrasts)
coefNames(res.dl_model1)

## to see the impact of Day night on gene expression after remove the age effect
model2 <- as.formula(" ~ 0  + cat_Day_Night + scale_age + (1 | Source) + (1 | prep) + (1 | pool) + scale(PMI) + log_n_counts + (1|Sex)")
contrasts <- c(Diff = "cat_Day_NightDay - cat_Day_NightNight")
res.proc_model2=processAssays(pb_subset,model2, min.count=5)
res.dl_model2=dreamlet(res.proc_model2,model2,contrasts = contrasts)
coefNames(res.dl_model2)

## to see the impact of Day night on gene expression as a function of age
model3 <- as.formula(" ~ 0  + cat_Day_Night*scale_age + (1 | Source) + (1 | prep) + (1 | pool) + scale(PMI) + log_n_counts + (1|Sex)")
res.proc_model3=processAssays(pb_subset,model3, min.count=5)
res.dl_model3=dreamlet(res.proc_model3,model3)
coefNames(res.dl_model3)

save(res.dl_model1,res.dl_model2,res.dl_model3,model1,model2,model3,file=paste0("/sc/arion/projects/psychAD/aging/sleep_patterns/analysis/dreamlet/",cell_groups,"_day_night_",celltype,"_",age_group,".RDATA"))
