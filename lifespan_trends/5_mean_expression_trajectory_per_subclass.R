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

good_plot=function(g1,sz,rot){
g1=g1+theme_classic()+theme(axis.title= element_text(color="black", size=sz),axis.text = element_text(color="black", size=sz),axis.text.x=element_text(angle=rot,hjust=1,vjust=1))
g1=g1+theme(legend.position="right")+theme(strip.text=element_text(face="bold",color="black"),strip.background=element_rect(fill="#eeeeee",color="#eeeeee"))
g1+theme(axis.title= element_text(color="black", size=sz),axis.text = element_text(color="black", size=sz))+theme(aspect.ratio=1)+theme(panel.border = element_rect(size = 1, color = "black", fill = NA))
}

best_model_list=list("poly_with_two_df_log_Age"= "~ 0 + poly(log2(Age+1),df=2)")

### this is the list of coefs from best model 
load(paste0("lifespan_trends/best_model/lifespan_coefs_dream_results_best_model.RDATA"))
lifespan_dream_results_df$celltype_ID=paste0(lifespan_dream_results_df$celltype,"_",lifespan_dream_results_df$ID)

### this is the file which has clusters assigned to the coefs from subclass_genes saved as celltype_ID
clusters_df=read.csv("lifespan_trends/best_model/clustered_data_nclust10.csv")
clusters_df$cluster=as.numeric(as.character(clusters_df$cluster))+1
clusters_df$ID=sub(".*_", "", clusters_df$celltype_ID)

### here dreamlet coefs from best model are merged with cluster ids and filtered for which we have FDR<.05
clusters_lifespan_dream_coefs_df=merge(lifespan_dream_results_df,clusters_df[,c("celltype_ID","cluster")],by='celltype_ID')
clusters_lifespan_dream_coefs_df=clusters_lifespan_dream_coefs_df[clusters_lifespan_dream_coefs_df$adj.P.Val<.05,]

### now we take residuals which has all the covars regressed out to fit best model and obtain lifespan trajectory per cluster per subclass

load(paste0("lifespan_trends/residuals/residuals_subclass.RDATA"))
metadata$SubID=rownames(metadata)

subclass_order=c("EN_L6_CT", "EN_L5_6_NP", "EN_L6B", "EN_L3_5_IT_1", "EN_L3_5_IT_2", "EN_L3_5_IT_3", "EN_L2_3_IT","EN_L6_IT_1",
    "EN_L6_IT_2","IN_LAMP5_RELN", "IN_LAMP5_LHX6","IN_ADARB2","IN_PVALB_CHC", "IN_PVALB", "IN_SST", "IN_VIP",
    "Astro","Oligo","OPC","Micro","Adaptive", "PVM","SMC","VLMC", "Endo","PC")

model_list=list("poly_with_two_df_log_Age"= "~ 0 + poly(log2(Age+1),df=2)")

aging_average_expression=list()
ct = 1

for (mm in (1:max(clusters_lifespan_dream_coefs_df$cluster))) {

aging_data_plot=list()
celltypes=unique(clusters_lifespan_dream_coefs_df$celltype[clusters_lifespan_dream_coefs_df$cluster==mm])

for (kk in (1:length(celltypes))) {

input_genes=clusters_lifespan_dream_coefs_df$ID[clusters_lifespan_dream_coefs_df$cluster==mm & clusters_lifespan_dream_coefs_df$celltype == celltypes[kk]]

## first we get residuals of celltypes[kk] and subset it for input_gene
    
res_mat_subset=res_mat[[which(names(res_mat) %in% celltypes[kk])]]
res_mat_subset_melted=melt(res_mat_subset[rownames(res_mat_subset) %in% input_genes,])
colnames(res_mat_subset_melted)[1:3]=c("ID","SubID","value")
data_expression_Age=merge(res_mat_subset_melted,metadata[,c("SubID","Age")],by="SubID")

if (length(input_gene) >1) {
coefs=clusters_lifespan_dream_coefs_df[clusters_lifespan_dream_coefs_df$cluster==mm & clusters_lifespan_dream_coefs_df$celltype == celltypes[kk],c("coef_1","coef_2")]

### get mean of all coefs across a cluster and a celltype
coefs_mean=apply(coefs,2,mean)

## get mean residualized expression of all genes per donor
df=data_expression_Age %>% dplyr::group_by(SubID)%>% dplyr::summarize(mean_expr=mean(value),Age=mean(Age))

## fit the mean residualized expr to the nonlinear age trend from polynomial model
    
fit=lm(paste0('mean_expr',model_list[["poly_with_two_df_log_Age"]]),data=df)
df$fitted_line=predict(fit)
x_axis=model.matrix(as.formula(model_list[["poly_with_two_df_log_Age"]]),data=df)
df$fitted_dream_line=x_axis %*% unlist(coefs_mean)[grep("Age", names(unlist(coefs_mean)),invert = TRUE)]
df$fitted_dream_line_scaled=df$fitted_dream_line-df$fitted_dream_line[which.min(df$Age)]
df=as.data.frame(df)
df$celltype=celltypes[kk]
df$cluster=mm

aging_average_expression[[ct]]=df
aging_average_expression_celltype[[kk]]=df
ct = ct + 1

} else {

next

}

}
aging_average_expression_celltype_df=as.data.frame(do.call(rbind,aging_average_expression_celltype))

### cell type trajectory per cluster (Supplementary Figure 6b)
    
g1=ggplot(aging_average_expression_celltype_df,aes(Age,mean_expr))+geom_point()+geom_line(data=aging_average_expression_celltype_df,aes(Age,fitted_dream_line,color=celltype))
print(good_plot(g1,12,0)+ggtitle(paste0("cluster_",mm))

}


aging_average_expression_df=do.call(rbind,aging_average_expression)
aging_average_expression_df_summary=aging_average_expression_df %>% dplyr::group_by(cluster,SubID) %>% dplyr::summarize(mean_fitted_dream_line=mean(fitted_dream_line_scaled),Age=mean(Age))

### Fig. 2a

pdf(paste0(data_dir"/lifespan_trends.pdf"))

print(ggplot(aging_average_expression_df_summary,aes(Age,mean_fitted_dream_line,color=factor(cluster),group=cluster))+geom_line()+geom_point(color="black")+facet_wrap(~cluster,ncol=5)+theme_bw()+theme(legend.position="none"))

dev.off()


