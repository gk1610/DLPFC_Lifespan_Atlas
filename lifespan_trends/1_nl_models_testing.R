devtools::install_github("DiseaseNeurogenomics/crumblr")

suppressPackageStartupMessages({
library(zellkonverter)
library(basilisk)
library(dreamlet)
library(crumblr)
library(Seurat)
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
library(RColorBrewer)
library(aplot)
library(ggtree)
library(mgcv)
library(splines)
library(knitr)
library(rmarkdown)
library(gridExtra)
library(aplot)
library(ggtree)})

subclass_order=c("EN_L6_CT", "EN_L5_6_NP", "EN_L6B", "EN_L3_5_IT_1", "EN_L3_5_IT_2", "EN_L3_5_IT_3", "EN_L2_3_IT","EN_L6_IT_1",
    "EN_L6_IT_2","IN_LAMP5_RELN", "IN_LAMP5_LHX6","IN_ADARB2","IN_VIP","IN_PVALB_CHC","IN_PVALB", "IN_SST",
    "Astro","Oligo","OPC","Micro","Adaptive", "PVM","SMC","VLMC", "Endo","PC")

good_plot=function(g1,sz,rot){
g1=g1+theme_classic()+theme(axis.title= element_text(color="black", size=sz),axis.text = element_text(color="black", size=sz),axis.text.x=element_text(angle=rot,hjust=1,vjust=1))
g1=g1+theme(legend.position="right")+theme(strip.text=element_text(face="bold",color="black"),strip.background=element_rect(fill="#eeeeee",color="#eeeeee"))
g1+theme(axis.title= element_text(color="black", size=sz),axis.text = element_text(color="black", size=sz))+theme(aspect.ratio=1)+theme(panel.border = element_rect(size = 1, color = "black", fill = NA))
}

lifespan_dir="lifespan_trends" 
dir.create(lifespan_dir)

data_dir="lifespan_trends/all_models" 
dir.create(data_dir)

syn62064718=synGet(entity="syn62064718",downloadLocation=data_dir)
pb=readRDS(paste0(data_dir,"/lifespan_pseudobulk.rds"))
colData(pb)$SubID=rownames(colData(pb))

## linear and non-linear age trends using following models

model_list=list("linear"="~ Age + Sex + Source + scale(PMI) + log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes",
"log_Age"="~ log2(Age+1) + Sex + Source + scale(PMI) + log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes",

"poly_with_two_df"= "~ poly(Age,df=2) + Source + Sex + scale(PMI)+ log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes",
"poly_with_two_df_log_Age"= "~ poly(log2(Age+1),df=2) + Source + Sex + scale(PMI) + log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes",

"poly_with_three_df"= "~ poly(Age,df=3) + Source + Sex + scale(PMI)+ log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes",
"poly_with_three_df_log_Age"= "~ poly(log2(Age+1),df=3) + Source + Sex + scale(PMI) + log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes",

"poly_with_four_df"= "~ poly(Age,df=4) + Source + Sex + scale(PMI)+ log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes",
"poly_with_four_df_log_Age"= "~ poly(log2(Age+1),df=4) + Source + Sex + scale(PMI) + log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes",

"poly_with_five_df"= "~ poly(Age,df=5) + Source + Sex + scale(PMI)+ log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes",
"poly_with_five_df_log_Age"= "~ poly(log2(Age+1),df=5) + Source + Sex + scale(PMI) + log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes",

"poly_with_six_df"= "~ poly(Age,df=6) + Source + Sex + scale(PMI)+ log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes",
"poly_with_six_df_log_Age"= "~ poly(log2(Age+1),df=6) + Source + Sex + scale(PMI) + log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes")


### dreamlet calculation (use bsub or qsub for each model)

for (kk in (1:length(model_list))) {
form=as.formula(model_list[[kk]])
res.proc = processAssays(pb,form, assays=grep("EN_L5_ET",assayNames(pb),value=TRUE,invert=TRUE))
res.dl = dreamlet(res.proc, form)
save(res.dl,form,file=paste0(data_dir,names(model_list)[kk],"_nl.RDATA"))
}


### After getting dreamlet object from all models, we extract BIC for each gene for 26 subclasses

files=list.files(data_dir,pattern="*_nl.RDATA",full.names=TRUE)
BIC_list=list()
for ( ii in (1:length(files))) {
load(files[ii])
df=lapply(1:length(res.dl),function(i) data.frame("BIC"=res.dl[[i]]$BIC,"celltype"=names(res.dl)[i],"genes"=names(res.dl[[i]]$BIC)))
BIC_list[[ii]]=do.call(rbind,df)
BIC_list[[ii]]$model_type=gsub("_nl.RDATA","",basename(files[ii]))
}

save(BIC_list,file=paste0(data_dir,"/BIC_all_models.RDATA"))


  
