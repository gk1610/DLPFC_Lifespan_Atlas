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

## code starts here

data_dir="Circadian/residuals" 
dir.create(data_dir)

syn62064718=synGet(entity="syn62064718",downloadLocation=data_dir)
pb=readRDS(paste0(data_dir,"/lifespan_pseudobulk.rds"))
colData(pb)$SubID=rownames(colData(pb))

### keep only young, middle and late adulthood samples as they have tod information available 
pb_subset=pb[,colData(pb)$Age>=20]
colData(pb_subset)$SubID=rownames(colData(pb_subset))

## optimal covariates model
form="~ scaled(Age) + Source + Sex + scale(PMI) + log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes"
res.proc = processAssays(pb_subset,as.formula(form))
res.dl = dreamlet(res.proc, as.formula(form))
res_mat = residuals(res.dl,res.proc)

save(res_mat,form,file=paste0(data_dir,"/residuals_subclass.RDATA"))
