suppressPackageStartupMessages({
library(variancePartition)
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

## code starts here 

data_dir="VariancePartition_gene_expression" 
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

### here varpart starts

for (input_groups in c("Developmental","Young_Adulthood","Middle_Adulthood","Late_Adulthood")) {

pb_subset_final=pb[,colData(pb)$groups==input_groups]

if (input_groups == "Developmental") {
form ="~ scaled_age + Sex + scale(PMI) + log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes"
} else {
form ="~ scaled_age + Source + Sex + scale(PMI) + log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes"
}

form=as.formula(form)
res.proc = processAssays(pb_subset_final,as.formula(form),assays=grep("EN_L5_ET",assayNames(pb_subset_final),value=TRUE,invert=TRUE))
vp.list =  fitVarPart(res.proc,as.formula(form))
save(vp.list,form,file=paste0(data_dir,"/varpart_age_",input_groups,".RDATA"))

}
