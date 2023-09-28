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

args = commandArgs(trailingOnly=TRUE)
celltype=args[1]
brain_bank=args[2]

################################################################################################################################
file="/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/h5ad_final/AGING_2023-04-01_16_14.h5ad"
sce = readH5AD(file, use_hdf5=TRUE, verbose=TRUE)
assayNames(sce)[1] = "counts"

sce1=subset(sce, ,subclass==celltype)

pb = aggregateToPseudoBulk(sce1,
      assay = "counts",
      cluster_id = "subclass",
      sample_id = "Channel")


df = as.data.frame(colData(pb))
df$Channel=rownames(df)
metadata_aggrs=as.data.frame(metadata(pb)$aggr_means)
metadata_df=merge(metadata_aggrs,df,by="Channel")
rownames(metadata_df)=metadata_df$Channel
metadata_df_flt=metadata_df[match(colnames(pb),metadata_df$Channel),]
identical(colnames(pb),as.character(metadata_df_flt$Channel))
metadata_df_flt$log_n_counts=log(metadata_df_flt$n_counts)
metadata_df_flt$scale_PMI=scale(metadata_df_flt$PMI)
colData(pb)=DataFrame(metadata_df_flt)
colData(pb)$new_poolID=paste0(sapply(sapply(as.character(colData(pb)$poolID),strsplit,"-"),`[`,1),"-",sapply(sapply(as.character(colData(pb)$poolID),strsplit,"-"),`[`,2))


if (brain_bank!="all") {

form = as.formula("~ (1|SubID) + (1|new_poolID) + (1|Sex) + scale_PMI + log_n_counts + (1|prep)")
pb_subset=pb[,colData(pb)$Source==brain_bank]

} else {

form = as.formula("~ (1|Source) + (1|SubID) + (1|new_poolID) + (1|Sex) + scale_PMI + log_n_counts + (1|prep)")
pb_subset=pb

}

res.proc = processAssays(pb_subset,form, min.count=5,assays=celltype)
res.dl = dreamlet(res.proc, form)
res_mat = residuals(res.dl[[1]],res.proc[[1]], type="pearson")
res_mat_non_pearson = residuals(res.dl[[1]],res.proc[[1]])

metadata_all=as.data.frame(colData(res.proc))
metadata_celltype=metadata[match(colnames(res_mat),metadata$Channel),]
identical(as.character(metadata_celltype$Channel),colnames(res_mat))

save(metadata_all,metadata_celltype,res_mat,res_mat_non_pearson,form,file=paste0("/sc/arion/projects/CommonMind/aging/analysis/dream/residuals/subclass/channel_",celltype,"_",brain_bank,"_samples_pearson_and_non_pearson_residuals.RDATA"))
