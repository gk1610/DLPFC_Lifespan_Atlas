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

form = as.formula("~ (1|Source) + (1|SubID) + (1|new_poolID) + (1|Sex) + scale_PMI + log_n_counts + (1|prep)")
pb_subset=pb[,colData(pb)$Age>=20,]]

res.proc = processAssays(pb_subset,form, min.count=5,assays=celltype)
res.dl = dreamlet(res.proc, form)
res_mat = residuals(res.dl[[1]],res.proc[[1]])

metadata_all=as.data.frame(colData(res.proc))
metadata_celltype=metadata_all[match(colnames(res_mat),metadata_all$Channel),]
identical(as.character(metadata_celltype$Channel),colnames(res_mat))
print(dim(metadata_celltype))

save(metadata_celltype,res_mat,form,file=paste0("/sc/arion/projects/psychAD/aging/Circadian/analysis/residuals/channel_",celltype,"_adulthood_old.RDATA"))
