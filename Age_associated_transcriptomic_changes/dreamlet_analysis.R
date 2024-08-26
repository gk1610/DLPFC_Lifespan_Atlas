devtools::install_github("DiseaseNeurogenomics/dreamlet")

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

groups_list=c("Developmental","Young_Adulthood","Middle_Adulthood","Late_Adulthood")
colors_groups_list=c("#F9CFA1","#EE9B00","#CA6702","#8C510A")
names(colors_groups_list)=groups_list

good_plot=function(g1,sz,rot){
g1=g1+theme_classic()+theme(axis.title= element_text(color="black", size=sz),axis.text = element_text(color="black", size=sz),axis.text.x=element_text(angle=rot,hjust=1,vjust=1))
g1=g1+theme(legend.position="right")+theme(strip.text=element_text(face="bold",color="black"),strip.background=element_rect(fill="#eeeeee",color="#eeeeee"))
g1+theme(axis.title= element_text(color="black", size=sz),axis.text = element_text(color="black", size=sz))+theme(aspect.ratio=1)+theme(panel.border = element_rect(size = 1, color = "black", fill = NA))
}

## code starts here 
data_dir="age_associated_transcriptomic_changes" 
dir.create(data_dir)
syn62064718=synGet(entity="syn62064718",downloadLocation=data_dir)
pb=readRDS(paste0(data_dir,"/lifespan_pseudobulk.rds"))
colData(pb)$SubID=rownames(colData(pb))

## age groups assignment
colData(pb)$groups="NA"
colData(pb)$groups[colData(pb)$Age<20]="Developmental"
colData(pb)$groups[colData(pb)$Age>=20 & colData(pb)$Age<40]="Young_Adulthood"
colData(pb)$groups[colData(pb)$Age>=40 & colData(pb)$Age<60]="Middle_Adulthood"
colData(pb)$groups[colData(pb)$Age>=60]="Late_Adulthood"
colData(pb)$scaled_age=scale(colData(pb)$Age)

### estimate coefficient of age from dreamlet analysis and runs mashr analysis using the summarys stats from dreamlet for four age groups

for (kk in (1:length(groups_list))) {
  
pb_subset_final=pb[,colData(pb)$groups==groups_list[kk]]

if (input_groups == "Developmental") {

form=" ~ scaled_age + Sex + scale(PMI) + log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes"
form=as.formula(form)

} else {

form="~ scaled_age + Source + Sex + scale(PMI) + log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes"
form=as.formula(form)

}

res.proc = processAssays(pb_subset_final,as.formula(form),assays=grep("EN_L5_ET",assayNames(pb_subset_final),value=TRUE,invert=TRUE))

## estimates coefficient of all variables using dreamlet per subclass 
res.dl = dreamlet(res.proc, form)

## mashr object for shared genes analysis
res_mash=run_mash(res.dl, coef="scaled_age") 

save(res_mash,res.dl,form,file=paste0("data_dir/",groups_list[kk],".RDATA"))

}


