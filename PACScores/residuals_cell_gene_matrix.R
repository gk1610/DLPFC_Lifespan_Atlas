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
library(Seurat)
})
library(data.table)

args = commandArgs(trailingOnly=TRUE)

cell_groups=args[1]
celltype=args[2]

file="/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/h5ad_final/AGING_2023-06-09_01_45.h5ad"
sce = readH5AD(file, use_hdf5=TRUE, verbose=TRUE)
assayNames(sce)[1] = "counts"
sce1=subset(sce, ,subclass==celltype & Source == "H")

colData(sce1)$scaled_PMI=scale(colData(sce1)$PMI)
colData(sce1)$prep_num=as.numeric(as.factor(colData(sce1)$prep))
colData(sce1)$pool_num=as.numeric(as.factor(colData(sce1)$pool))
colData(sce1)$Sex_num=as.numeric(as.factor(colData(sce1)$Sex))
colData(sce1)$SubID_num=as.numeric(as.factor(colData(sce1)$SubID))

sce_subset <- SingleCellExperiment(assays = list(counts = as.matrix(sce1@assays@data$counts)), colData=colData(sce1))
sce_seurat <- as.Seurat(sce_subset, counts = "counts",data=NULL)
sce_seurat_sc <- SCTransform(sce_seurat, vars.to.regress = c("scaled_PMI","prep_num","pool_num","Sex_num","SubID_num"),residual.features=NULL,return.only.var.genes=F,verbose = TRUE)

save(sce_seurat_sc,file=paste0("/sc/arion/projects/psychAD/aging/PACscores/residuals/",celltype,"_sct.RDATA"))
