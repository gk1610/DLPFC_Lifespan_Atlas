.libPaths(c("/sc/arion/projects/psychAD/aging/kiran/RLib_4_3.2",.libPaths()))
library(zellkonverter)
library(dreamlet)
library(variancePartition)
library(HDF5Array)
library(dplyr)
source("/sc/arion/work/girdhk01/scripts/chipseq_files.R")
source('/sc/arion/projects/roussp01a/pengfei/hicchip/scripts/hic_helper.R')
source("/sc/arion/work/girdhk01/scripts/myscripts/CRD/CMC_SV_help.R")
source("/sc/arion/work/girdhk01/scripts/myscripts/CRD/TRH_help_functions.R")
source("/sc/arion/work/girdhk01/scripts/myscripts/CRD/HiC_pengfei.R")
source("/sc/arion/work/girdhk01/scripts/myscripts/chiq_test.R")

args = commandArgs(trailingOnly=TRUE)
celltypes=args[1]

average_model_list=list("poly_with_two_df_log_Age"= "~ 0 + poly(log2(Age+1),df=2)")


load("/sc/arion/projects/psychAD/aging/kiran/analysis/lifespan/subclass/best_model/lifespan_coefs_dream_results_best_model.RDATA")
load(paste0("/sc/arion/projects/psychAD/aging/kiran/analysis/residuals/residuals_subid_subclass_all.RDATA"))
metadata$SubID=rownames(metadata)

input_gene=lifespan_dream_coefs_df$ID[lifespan_dream_coefs_df$celltype == celltypes]
aging_dream_trajectories=list()
aging_dream_trajectories_scaled=list()
ct = 1

print(length(input_gene))

for (kkk in (1:length(input_gene))) {

coefs=lifespan_dream_coefs_df[lifespan_dream_coefs_df$ID %in% input_gene[kkk] & lifespan_dream_coefs_df$celltype == celltypes,c("coef_1","coef_2","Intercept")]

res_mat_subset=res_mat[[which(names(res_mat) %in% celltypes)]]
res_mat_subset_melted=reshape2::melt(res_mat_subset[rownames(res_mat_subset) %in% input_gene[kkk],])
res_mat_subset_melted$SubID=rownames(res_mat_subset_melted)
res_mat_subset_melted$ID=input_gene[kkk]
data_expression_Age=merge(res_mat_subset_melted,metadata[,c("SubID","Age")],by="SubID")

df=data_expression_Age
x_axis=model.matrix(as.formula(average_model_list[["poly_with_two_df_log_Age"]]),data=df)
df$fitted_dream_line=x_axis %*% unlist(coefs)[grep("Intercept", names(unlist(coefs)),invert = TRUE)]

fitted_dream_line=as.vector(x_axis %*% unlist(coefs)[grep("Intercept", names(unlist(coefs)),invert = TRUE)])
names(fitted_dream_line)=df$SubID

aging_dream_trajectories[[kkk]]=fitted_dream_line
aging_dream_trajectories_scaled[[kkk]]=fitted_dream_line-fitted_dream_line[which.min(df$Age)]
ct = ct + 1

}

names(aging_dream_trajectories)=input_gene
names(aging_dream_trajectories_scaled)=input_gene

aging_dream_trajectories_mat=do.call(rbind,aging_dream_trajectories)
aging_dream_trajectories_scaled_mat=do.call(rbind,aging_dream_trajectories_scaled)
res_mat_subset=res_mat[[which(names(res_mat) %in% celltypes)]]
res_mat=res_mat_subset

load(paste0("/sc/arion/projects/psychAD/aging/kiran/analysis/residuals/residuals_subid_subclass_all.RDATA"))
metadata$SubID=rownames(metadata)

save(metadata,res_mat,aging_dream_trajectories_scaled_mat,aging_dream_trajectories_mat,file=paste0("/sc/arion/projects/psychAD/aging/kiran/analysis/lifespan/subclass/shinyapp/average_trajectory_",celltypes,".RDATA"))
