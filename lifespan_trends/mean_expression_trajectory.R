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

best_model_list=list("poly_with_two_df_log_Age"= "~ 0 + poly(log2(Age+1),df=2)")

load(paste0("lifespan_trends/best_model/lifespan_coefs_dream_results_best_model.RDATA"))
lifespan_dream_results_df$celltype_ID=paste0(lifespan_dream_results_df$celltype,"_",lifespan_dream_results_df$ID)

clusters_df=read.csv("lifespan_trends/best_model/clustered_data_nclust10.csv)
clusters_df$cluster=as.numeric(as.character(clusters_df$cluster))
clusters_df$ID=sub(".*_", "", clusters_df$celltype_ID)

clusters_lifespan_dream_coefs_df=merge(lifespan_dream_results_df,clusters_df[,c("celltype_ID","cluster")],by='celltype_ID')
clusters_lifespan_dream_coefs_df=clusters_lifespan_dream_coefs_df[clusters_lifespan_dream_coefs_df$adj.P.Val<.05,]


load(paste0("lifespan_trends/residuals/residuals_subclass.RDATA"))
metadata$SubID=rownames(metadata)

subclass_order=c("EN_L6_CT", "EN_L5_6_NP", "EN_L6B", "EN_L3_5_IT_1", "EN_L3_5_IT_2", "EN_L3_5_IT_3", "EN_L2_3_IT","EN_L6_IT_1",
    "EN_L6_IT_2","IN_LAMP5_RELN", "IN_LAMP5_LHX6","IN_ADARB2","IN_PVALB_CHC", "IN_PVALB", "IN_SST", "IN_VIP",
    "Astro","Oligo","OPC","Micro","Adaptive", "PVM","SMC","VLMC", "Endo","PC")

average_model_list=list("poly_with_two_df_log_Age"= "~ 0 + poly(log2(Age+1),df=2)")

aging_average_expression=list()
ct = 1

file_name=gsub('.csv','.pdf',input_file)

pdf(file_name,width=15)

ggplot(clusters_lifespan_dream_coefs_df,aes(coef_1,coef_2,color=factor(cluster)))+geom_point()+theme_bw()

for (mm in (1:max(clusters_lifespan_dream_coefs_df$cluster))) {

aging_data_plot=list()
celltypes=unique(clusters_lifespan_dream_coefs_df$celltype[clusters_lifespan_dream_coefs_df$cluster==mm])

for (kk in (1:length(celltypes))) {

print(mm)

input_gene=clusters_lifespan_dream_coefs_df$ID[clusters_lifespan_dream_coefs_df$cluster==mm & clusters_lifespan_dream_coefs_df$celltype == celltypes[kk]]

if (length(input_gene) >1) {
coefs=clusters_lifespan_dream_coefs_df[clusters_lifespan_dream_coefs_df$cluster==mm & clusters_lifespan_dream_coefs_df$celltype == celltypes[kk],c("coef_1","coef_2")]
coefs_mean=apply(coefs,2,mean)
print(length(input_gene))

res_mat_subset=res_mat[[which(names(res_mat) %in% celltypes[kk])]]
res_mat_subset_melted=melt(res_mat_subset[rownames(res_mat_subset) %in% input_gene,])
colnames(res_mat_subset_melted)[1:3]=c("ID","SubID","value")
data_expression_Age=merge(res_mat_subset_melted,metadata[,c("SubID","Age")],by="SubID")

df=data_expression_Age %>% dplyr::group_by(SubID)%>% dplyr::summarize(mean_expr=mean(value),Age=mean(Age))
fit=lm(paste0('mean_expr',average_model_list[["poly_with_two_df_log_Age"]]),data=df)
df$fitted_line=predict(fit)
x_axis=model.matrix(as.formula(average_model_list[["poly_with_two_df_log_Age"]]),data=df)
df$fitted_dream_line=x_axis %*% unlist(coefs_mean)[grep("PMI|n_genes|PC|Intercept", names(unlist(coefs_mean)),invert = TRUE)]
df$fitted_dream_line_scaled=df$fitted_dream_line-df$fitted_dream_line[which.min(df$Age)]
df=as.data.frame(df)
df$celltype=celltypes[kk]
  df$cluster=mm

aging_average_expression[[ct]]=df
aging_data_plot[[kk]]=df
ct = ct + 1

} else {

next

}

}
aging_data_plot_df=as.data.frame(do.call(rbind,aging_data_plot))

g1=ggplot(aging_data_plot_df,aes(Age,mean_expr))+geom_point()+geom_line(data=aging_data_plot_df,aes(Age,fitted_dream_line,color=celltype))
print(good_plot(g1,12,0)+ggtitle(paste0("cluster_",mm))+scale_color_manual(values=subclass_color_map))

g1=ggplot(aging_data_plot_df,aes(Age,mean_expr))+geom_line(data=aging_data_plot_df,aes(Age,fitted_dream_line,color=celltype))
print(good_plot(g1,12,0)+ggtitle(paste0("cluster_",mm))+scale_color_manual(values=subclass_color_map))+ylab("Predicted expression per celltype using mean of clustered coefficients of genes")

}

aging_average_expression_df=do.call(rbind,aging_average_expression)
aging_average_expression_df_summary=aging_average_expression_df %>% dplyr::group_by(cluster,SubID) %>% dplyr::summarize(mean_fitted_dream_line=mean(fitted_dream_line_scaled),Age=mean(Age))
aging_average_expression_df_summary$cluster=factor(aging_average_expression_df_summary$cluster,levels=cluster_order)

print(ggplot(aging_average_expression_df_summary,aes(Age,mean_fitted_dream_line,color=factor(cluster),group=cluster))+geom_line()+geom_point(color="black")+facet_wrap2(~cluster,ncol=5,strip=strip)+theme_bw()+theme(legend.position="none"))

dev.off()


