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

celltype=args[1]
ii=as.numeric(args[2])


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
#Group 4_Early Adulthood(21-40):   Day - n = 18    Night - n = 23
#Group 5_Late Adulthood(41-60):    Day - n = 47    Night - n = 16
#Group 6_Old Age (>60):            Day - n = 51    Night - n = 25

#1.2 DxN of Combined Groups 
#Group 4-5_Adulthood(21-60):       Day - n = 65    Night - n = 39
#Group 4-6_Adulthood+Old Age(>20)  Day - n = 116   Night - n = 66

contrast_age_groups_list=list(
"early_adulthood"=df$SubID[df$Age>=20 & df$Age<40],
"late_adulthood"=df$SubID[df$Age>=40 & df$Age<60],
"old" = df$SubID[df$Age>=60],
"adulthood"=df$SubID[df$Age>=20 & df$Age<60],
"adulthood_old"=df$SubID[df$Age>=20])

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

metadata_daynight_merged_subset=metadata_daynight_merged[metadata_daynight_merged$SubID %in% contrast_age_groups_list[[ii]],]

### subsetting pseudobulk for a given age group 

keep_subids=colnames(pb)[match(metadata_daynight_merged_subset$SubID,colnames(pb))]
pb_subset=subset(pb, ,SubID %in% keep_subids)

metadata_daynight_merged_ordered = metadata_daynight_merged_subset[match(colnames(pb_subset),as.character(metadata_daynight_merged_subset$SubID)),]
identical(as.character(metadata_daynight_merged_ordered$SubID),colnames(pb_subset))
rownames(metadata_daynight_merged_ordered)=metadata_daynight_merged_ordered$SubID
colData(pb_subset)=DataFrame(metadata_daynight_merged_ordered)

### making the day night contrast as factor
metadata_daynight_merged_ordered$cat_Day_Night=factor(metadata_daynight_merged_ordered$Day_Night,levels=c(1,2),labels=c("Day","Night"))
table(metadata_daynight_merged_ordered$cat_Day_Night)


if (uniqueN(colData(pb_subset)$Source)>1) {

form <- as.formula(" ~ 0  + cat_Day_Night + (1 | Source) + (1 | prep) + (1 | pool) + scale(PMI) + log_n_counts + (1|Sex)")

} else {

form <- as.formula(" ~ 0 + cat_Day_Night + (1 | prep) + (1 | pool) + scale(PMI) + log_n_counts + (1|Sex)")

}

res.proc=processAssays(pb_subset,form, min.count=5)

contrasts <- c(Diff = "cat_Day_NightDay - cat_Day_NightNight")
res.dl=dreamlet(res.proc,form),contrasts = contrasts)
coefNames(res.dl)

#results <- topTable(res.dl[[1]], coef = "Diff")

save(res.dl,form,metadata_daynight_merged_ordered,file=paste0("/sc/arion/projects/psychAD/aging/sleep_patterns/dreamlet/Subid_day_night_",names(contrast_age_groups_list)[ii],".RDATA"))






