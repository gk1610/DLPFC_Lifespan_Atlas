devtools::install_github("DiseaseNeurogenomics/crumblr")
devtools::install_github("DiseaseNeurogenomics/variancePartition")

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

data_dir="nuclei_composition" 
dir.create(data_dir)
syn62064718=synGet(entity="syn62064718",downloadLocation=data_dir)
pb=readRDS(paste0(data_dir,"/lifespan_pseudobulk.rds"))
colData(pb)$SubID=rownames(colData(pb))

## age groups assignment
pb_subset=pb
colData(pb_subset)$groups="NA"
colData(pb_subset)$groups[colData(pb_subset)$Age<1]="Developmental"
colData(pb_subset)$groups[colData(pb_subset)$Age>=1 & colData(pb_subset)$Age<12]="Developmental"
colData(pb_subset)$groups[colData(pb_subset)$Age>=12 & colData(pb_subset)$Age<20]="Developmental"
colData(pb_subset)$groups[colData(pb_subset)$Age>=20 & colData(pb_subset)$Age<40]="Young_Adulthood"
colData(pb_subset)$groups[colData(pb_subset)$Age>=40 & colData(pb_subset)$Age<60]="Middle_Adulthood"
colData(pb_subset)$groups[colData(pb_subset)$Age>=60]="Late_Adulthood"
colData(pb_subset)$scaled_age=scale(colData(pb_subset)$Age)

## remove the subclass which has nuclei less than < 500
cell_counts=cellCounts(pb_subset)
cell_counts=cell_counts[,colnames(cell_counts)!="EN_L5_ET"]

## nuclei composition estimation starts here
cobj = crumblr(cell_counts)

## covariate model
bestModel="~ Sex + Source + scale(PMI)"
form = as.formula(bestModel)

## regressing out covariates
fit = dream(cobj, form, colData(pb_subset))
residuals_mat=fit$residuals

cobj_df=as.data.frame(t(as.matrix(residuals_mat)))
cobj_df$SubID=rownames(cobj_df)
cobj_df=cobj_df[,grep("EN_L5_ET",colnames(cobj_df),invert=TRUE)]

metadata_df=as.data.frame(colData(pb_subset))
metadata_df$SubID=rownames(metadata_df)

cobj_metadata=merge(cobj_df,metadata_df[,c("Age","SubID")],by="SubID")
cobj_metadata=as.data.frame(cobj_metadata)
cobj_metadata_groups=add_groups_metadata(cobj_metadata)
cobj_metadata_groups_melted=melt(cobj_metadata_groups,id=c("Age","SubID","groups"))

## estimation of BIC for linear and log age trends

BIC_list=list()
ct = 1
for (cell_types in subclass_order) {

cobj_metadata$value=cobj_metadata[,cell_types]
df_order=cobj_metadata[order(cobj_metadata$Age),]
df_order=as.data.frame(df_order)
BIC_list[[ct]]=data.frame(celltype=cell_types)
BIC_list[[ct]]$log_BIC_model=BIC(lm(value~log2(Age+1),data=df_order))
BIC_list[[ct]]$linear_BIC_model=BIC(lm(value~(Age),data=df_order))
ct = ct + 1

}


BIC_list1=do.call(rbind,BIC_list)
BIC_list1$diff_BIC=BIC_list1$log_BIC_model-BIC_list1$linear_BIC_model
BIC_list_plot=melt(BIC_list1,id="celltype")
BIC_list_plot$celltype=factor(BIC_list_plot$celltype,levels=subclass_order)


save(BIC_list_plot,file=paste0(data_dir,"lifespan_trajectories_all_celltypes_BIC_models.RDATA"))

## plot Fig.S2a

pdf(paste0(data_dir,"/lifespan_trajectories_all_celltypes_BIC_models.pdf"))

g1=ggplot(subset(BIC_list_plot,variable=="diff_BIC"),aes(x=celltype, y=value,fill=factor(celltype))) + geom_bar(stat="identity",position="dodge2") + theme_bw()
g1+theme(legend.position="top")+xlab("")+ylab("BIC(log-linear)")

dev.off()


## plot Fig.S2b

### kmeans clustering of coefs of age from log model

crumblr_coefs_list=list()
crumblr_fitted_line_list=list()
ct = 1

for (cell_types in subclass_order) {

cobj_metadata$value=cobj_metadata[,cell_types]
df_order=cobj_metadata[order(cobj_metadata$Age),]
df_order$log_model=fitted(lm(value~log2(Age+1),data=df_order))
df_order$samples=df_order$SubID

crumblr_fitted_line_list[[ct]]=df_order[,c(cell_types,"samples","value","Age","log_model")]
colnames(crumblr_fitted_line_list[[ct]])[1]="cell_composition"
crumblr_fitted_line_list[[ct]]$celltype=cell_types
crumblr_coefs_list[[ct]]=data.frame(celltype=cell_types)
crumblr_coefs_list[[ct]]$log_estimate=summary(lm(value~log2(Age+1),data=df_order))[[4]][2,1]
crumblr_coefs_list[[ct]]$log_pvalue=summary(lm(value~log2(Age+1),data=df_order))[[4]][2,4]

ct = ct + 1
}


crumblr_coefs_list1=do.call(rbind,crumblr_coefs_list)
rownames(crumblr_coefs_list1)=crumblr_coefs_list1$celltype

crumblr_fitted_line_list1=do.call(rbind,crumblr_fitted_line_list)
rownames(crumblr_fitted_line_list1)=NULL

crumblr_coefs_list1$adj_pvalue=p.adjust(crumblr_coefs_list1$log_pvalue,method="fdr")

k <- 2  # number of clusters
kmeans_result <- kmeans(crumblr_coefs_list1$log_estimate, centers = k)
crumblr_coefs_list1$cluster=kmeans_result$cluster

pdf(paste0(data_dir,"/kmeans_results_lifespan_trajectories_all_celltypes.pdf"))

g1=ggplot(crumblr_coefs_list1,aes(factor(cluster),log_estimate,color=factor(celltype)))+geom_point()
g1+geom_hline(yintercept=0,linetype="dashed")+theme(legend.position="top")

dev.off()

crumblr_fitted_line_list_cluster=merge(crumblr_fitted_line_list1,crumblr_coefs_list1,by="celltype")
crumblr_fitted_line_list_cluster$celltype=factor(crumblr_fitted_line_list_cluster$celltype,level=subclass_order)

### make plot of whole lifespan trajectories of nuclei composition 

## plot Extended Data Figure 2a

pdf(paste0(data_dir,"/lifespan_trajectories_all_celltypes.pdf"))

g1=ggplot(subset(crumblr_fitted_line_list_cluster,adj_pvalue<.05)) + geom_line(aes(x=Age, y=log_model,group=celltype,color=celltype))+
                 facet_wrap(~cluster,ncol=3)+scale_color_manual(values=subclass_color_map)+theme_bw()+theme(legend.position="top")

g1=good_plot(g1,10,0)+theme(legend.position="top")+xlab("Age")+ylab("Residualized_counts")+geom_vline(xintercept=c(20,40,60),linetype="dashed")

print(g1)

dev.off()


## plot Extended Data Figure 2b

pdf(paste0(data_dir,"/lifespan_trajectories_all_celltypes_log_models.pdf"))

for (cell_types in subclass_order) {

cobj_metadata$value=cobj_metadata[,cell_types]
df_order=cobj_metadata[order(cobj_metadata$Age),]
df_order$log_model=fitted(lm(value~log2(Age+1),data=df_order))
df_order=as.data.frame(df_order)

g1=ggplot(df_order) +
                 geom_point(aes(x=Age, y=value),alpha=0.55, color="black") + geom_point(aes(x=Age, y=value),alpha=0.55, color="black") + geom_line(aes(x=Age, y=log_model),alpha=0.55, color="blue")+
                 theme_minimal() +ggtitle(paste0("log_fit_",cell_types))

g1=g1+theme(legend.position="top")+xlab("")+ylab("Residualized_counts")

print(g1)

}

dev.off()


#### quantification of coefficent of Age using crumblr limma based analysis as shown in Extended Data Figure 2c

bestModel=" ~ log2(Age + 1)"
form = as.formula(bestModel)
fit = dream(residuals_mat, form, metadata_df) ## we use here residuals after removing covariates 
fit = eBayes(fit) 
head(fit$coefficients)

df_table=topTable(fit, coef='log2(Age + 1)', number=Inf) %>%   
  select(logFC, AveExpr, t, P.Value, adj.P.Val)
head(df_table)

## multivariate hypothesis testing across all celltypes

# first we measure the correlation across all celltypes 
hcl = buildClusterTreeFromPB(pb_subset, assays = grep("EN_L5_ET",assayNames(pb),invert=TRUE,value=TRUE))
cobj$E=residuals_mat
res1 = treeTest(fit, cobj, hcl, coef="log2(Age + 1)")
fig.tree = plotTreeTest(res1) + theme(legend.position="left")+xlim(0, 15)
tab_scaled_age = topTable(fit, "log2(Age + 1)", number=Inf, sort.by="none")
tab_scaled_age$celltype =factor(rownames(tab_scaled_age), rev(subclass_order))
tab_scaled_age$se = with(tab_scaled_age, logFC/ t)
fig.tree$data$label=factor(fig.tree$data$label,levels=rev(subclass_order))


## plot Extended Data Figure 2c

pdf(paste0(data_dir,"/wls_with_res.pdf"))

fig.logFC = ggplot(tab_scaled_age, aes(celltype, logFC,color=celltype)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey50") + 
  geom_errorbar(aes(ymin=logFC - 1.96*se, ymax=logFC + 1.96*se), width=0) +
  geom_point(data=tab_scaled_age,aes(celltype, logFC,color=celltype,size=-log10(adj.P.Val))) +
  coord_flip() +
  theme_classic() +
  xlab('') +
  theme(aspect.ratio=3.65,axis.text.y = element_blank())+scale_color_manual(values=subclass_color_map)
 fig.logFC %>% insert_left(fig.tree)  #%>% insert_right(fig.vp)

dev.off()

write.csv(tab_scaled_age,file=paste0(data_dir,"wls_with_res.csv"))


#### age group specific analysis ####

# measure the variance explained by log2(Age+1) for each age group in nuclei composition

pb_subset_final=pb_subset[,colData(pb_subset)$groups=="Developmental"]
cobj = crumblr(cellCounts(pb_subset_final))
bestModel=" ~ scale(PMI) + Sex"
form = as.formula(bestModel)
fit = dream(cobj, form, colData(pb_subset_final))
residuals_mat=fit$residuals
bestModel=" ~ log2(Age+1)"
form=as.formula(bestModel)
res.vp1 = fitExtractVarPartModel(residuals_mat, form, colData(pb_subset_final))
res.vp1=res.vp1[grep("EN_L5_ET",rownames(res.vp1),invert=TRUE),]
res.vp1=as.data.frame(res.vp1)

pb_subset_final=pb_subset[,colData(pb_subset)$groups=="Young_Adulthood"]
colData(pb_subset_final)$scaled_age=scale(colData(pb_subset_final)$Age)
cobj = crumblr(cellCounts(pb_subset_final))
bestModel=" ~ scale(PMI) + Sex + Source"
form = as.formula(bestModel)
fit = dream(cobj, form, colData(pb_subset_final))
residuals_mat=fit$residuals
bestModel=" ~ log2(Age+1)"
form=as.formula(bestModel)
res.vp2 = fitExtractVarPartModel(residuals_mat, form, colData(pb_subset_final))
res.vp2=res.vp2[grep("EN_L5_ET",rownames(res.vp2),invert=TRUE),]
res.vp2=as.data.frame(res.vp2)

pb_subset_final=pb_subset[,colData(pb_subset)$groups=="Middle_Adulthood"]
colData(pb_subset_final)$scaled_age=scale(colData(pb_subset_final)$Age)
cobj = crumblr(cellCounts(pb_subset_final))
bestModel=" ~ scale(PMI) + Sex + Source"
form = as.formula(bestModel)
fit = dream(cobj, form, colData(pb_subset_final))
residuals_mat=fit$residuals
bestModel=" ~ log2(Age+1)"
form=as.formula(bestModel)
res.vp3 = fitExtractVarPartModel(residuals_mat, form, colData(pb_subset_final))
res.vp3=res.vp3[grep("EN_L5_ET",rownames(res.vp3),invert=TRUE),]
res.vp3=as.data.frame(res.vp3)

pb_subset_final=pb_subset[,colData(pb_subset)$groups=="Late_Adulthood"]
colData(pb_subset_final)$scaled_age=scale(colData(pb_subset_final)$Age)
cobj = crumblr(cellCounts(pb_subset_final))
bestModel=" ~ scale(PMI) + Sex + Source"
form = as.formula(bestModel)
fit = dream(cobj, form, colData(pb_subset_final))
residuals_mat=fit$residuals
bestModel=" ~ log2(Age+1)"
form=as.formula(bestModel)
res.vp4 = fitExtractVarPartModel(residuals_mat, form, colData(pb_subset_final))
res.vp4=res.vp4[grep("EN_L5_ET",rownames(res.vp4),invert=TRUE),]
res.vp4=as.data.frame(res.vp4)

res.vp1$celltype=rownames(res.vp1)
res.vp2$celltype=rownames(res.vp2)
res.vp3$celltype=rownames(res.vp3)
res.vp4$celltype=rownames(res.vp4)

res.vp1$groups="Developmental"
res.vp2$groups="Young_Adulthood"
res.vp3$groups="Middle_Adulthood"
res.vp4$groups="Late_Adulthood"

res.vp_all=cbind(res.vp1[,c("log2(Age + 1)","groups","celltype")],res.vp2[,c("log2(Age + 1)","groups","celltype")],res.vp3[,c("log2(Age + 1)","groups","celltype")],res.vp4[,c("log2(Age + 1)","groups","celltype")])
res.vp_all=res.vp_all[,grep("log",colnames(res.vp_all))]
colnames(res.vp_all)=c("Developmental","Young_Adulthood","Middle_Adulthood","Late_Adulthood")

## plot Extended Data Figure 2d

pdf(paste0(data_dir,"/group_specific_varpart_all_celltypes.pdf"))

plotVarPart(res.vp_all)+scale_fill_manual(values=colors_groups_list)

dev.off()

### here we quantify the coefficient of log2(Age+1) for each group using limma-based crumblr analysis

tab_scaled_age_list=list()
ct = 1

for (group in groups_list) {

pb_subset_final=pb_subset[,colData(pb_subset)$groups==group]
metadata_df=as.data.frame(colData(pb_subset_final))
cell_counts = cellCounts(pb_subset_final)
cell_counts=cell_counts[,colnames(cell_counts)!="EN_L5_ET"]
cobj = crumblr(cell_counts)

if (group =="Developmental"){
bestModel=" ~ log2(Age+1) + scale(PMI) + Sex"
} else {
bestModel=" ~ log2(Age+1) + scale(PMI) + Source + Sex"
}

form = as.formula(bestModel)
fit = dream(cobj, form, metadata_df)
fit = eBayes(fit) 
head(fit$coefficients)

hcl = buildClusterTreeFromPB(pb_subset_final, assays = grep("EN_L5_ET",assayNames(pb_subset_final),invert=TRUE,value=TRUE))

res1 = treeTest(fit, cobj, hcl, coef="log2(Age + 1)")
fig.tree = plotTreeTest(res1) + theme(legend.position="none")+xlim(0, 15)
tab_scaled_age = topTable(fit, "log2(Age + 1)", number=Inf, sort.by="none")
tab_scaled_age$celltype = factor(rownames(tab_scaled_age), rev(subclass_order))
tab_scaled_age$se = with(tab_scaled_age, logFC/ t)
tab_scaled_age$groups=group
tab_scaled_age_list[[ct]]=tab_scaled_age

write.csv(tab_scaled_age,file=paste0(data_dir,"/",group,"_crumblr_group_specific_results.csv"))

ct = ct + 1

}


tab_scaled_age_list1=do.call(rbind,tab_scaled_age_list)
tab_scaled_age_list1$groups=factor(tab_scaled_age_list1$groups,rev(groups_list))
tab_scaled_age_list1$celltype = factor(tab_scaled_age_list1$celltype, c(subclass_order))
tab_scaled_age_list1$mylabel=""
tab_scaled_age_list1$mylabel[tab_scaled_age_list1$P.Value<.05]="*"
tab_scaled_age_list1$mylabel[tab_scaled_age_list1$adj.P.Val<.05]="#"


## plot Extended Data Figure 2e

pdf(paste0(data_dir,"/group_specific_dream_all_celltypes.pdf"))

ggplot(tab_scaled_age_list1, aes(celltype, groups, fill=logFC, label=mylabel)) +
  geom_tile() +
  geom_text(vjust=1, hjust=0.5) +
  theme_classic() + 
  theme(aspect.ratio=0.3, axis.text.x = element_text(angle = 45,hjust=1,vjust=1)) +
  scale_fill_gradient2(low="blue", mid="white", high="red") 

dev.off()







