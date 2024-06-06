
##### Trajectories estimation starts here

setwd("/sc/arion/projects/CommonMind/aging/analysis/crumblr")

pb=readRDS("/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/pseudobulk/AGING_2024-02-01_22_23_PB_SubID_subclass.RDS")
colData(pb)$SubID=rownames(colData(pb))

#### keep final samples
samples_to_keep=read.table("/sc/arion/projects/CommonMind/aging/resources/AGING_2024-02-01_22_23_processAssays_SubID_subclass.txt")
pb_subset=pb[,colData(pb)$SubID %in% samples_to_keep$V1]

colData(pb_subset)$groups="NA"
colData(pb_subset)$groups[colData(pb_subset)$Age<1]="Childhood"
colData(pb_subset)$groups[colData(pb_subset)$Age>=1 & colData(pb_subset)$Age<12]="Childhood"
colData(pb_subset)$groups[colData(pb_subset)$Age>=12 & colData(pb_subset)$Age<20]="Childhood"
colData(pb_subset)$groups[colData(pb_subset)$Age>=20 & colData(pb_subset)$Age<40]="Young_Adulthood"
colData(pb_subset)$groups[colData(pb_subset)$Age>=40 & colData(pb_subset)$Age<60]="Middle_Adulthood"
colData(pb_subset)$groups[colData(pb_subset)$Age>=60]="Late_Adulthood"
colData(pb_subset)$scaled_age=scale(colData(pb_subset)$Age)

cell_counts=cellCounts(pb_subset)
cell_counts=cell_counts[,colnames(cell_counts)!="EN_L5_ET"]
cobj = crumblr(cell_counts)

### cell composition for dreamlet and crumblr analysis
bestModel="~ Sex + Source + scale(PMI)"
form = as.formula(bestModel)
fit = dream(cobj, form, colData(pb_subset))
residuals_mat=fit$residuals
cobj_df=as.data.frame(t(as.matrix(residuals_mat)))
cobj_df$SubID=rownames(cobj_df)
cobj_df=cobj_df[,grep("EN_L5_ET|EN_NF",colnames(cobj_df),invert=TRUE)]
metadata_df=as.data.frame(colData(pb_subset))
metadata_df$SubID=rownames(metadata_df)
cobj_metadata=merge(cobj_df,metadata_df[,c("Age","SubID")],by="SubID")
cobj_metadata=as.data.frame(cobj_metadata)
cobj_metadata_groups=add_groups_metadata(cobj_metadata)
cobj_metadata_groups_melted=melt(cobj_metadata_groups,id=c("Age","SubID","groups"))
cobj_metadata_groups_melted=add_class2_metadata(cobj_metadata_groups_melted,"variable")
table(cobj_metadata_groups_melted$class)


# test which function has best BIC

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
library(gridExtra)

save(BIC_list_plot,file="/sc/arion/projects/psychAD/aging/kiran/analysis/crumblr/lifespan/lifespan_trajectories_all_celltypes_BIC_models.RDATA")



#### Whole life span crumblr analysis

setwd("/sc/arion/projects/CommonMind/aging/analysis/crumblr")

pb=readRDS("/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/pseudobulk/AGING_2024-02-01_22_23_PB_SubID_subclass.RDS")
colData(pb)$SubID=rownames(colData(pb))

#### keep final samples
samples_to_keep=read.table("/sc/arion/projects/CommonMind/aging/resources/AGING_2024-02-01_22_23_processAssays_SubID_subclass.txt")
pb_subset=pb[,colData(pb)$SubID %in% samples_to_keep$V1]

colData(pb_subset)$groups="NA"
colData(pb_subset)$groups[colData(pb_subset)$Age<1]="Childhood"
colData(pb_subset)$groups[colData(pb_subset)$Age>=1 & colData(pb_subset)$Age<12]="Childhood"
colData(pb_subset)$groups[colData(pb_subset)$Age>=12 & colData(pb_subset)$Age<20]="Childhood"
colData(pb_subset)$groups[colData(pb_subset)$Age>=20 & colData(pb_subset)$Age<40]="Young_Adulthood"
colData(pb_subset)$groups[colData(pb_subset)$Age>=40 & colData(pb_subset)$Age<60]="Middle_Adulthood"
colData(pb_subset)$groups[colData(pb_subset)$Age>=60]="Late_Adulthood"
colData(pb_subset)$scaled_age=scale(colData(pb_subset)$Age)

metadata_df=as.data.frame(colData(pb_subset))
cell_counts=cellCounts(pb_subset)
cell_counts=cell_counts[,colnames(cell_counts)!="EN_L5_ET"]
cobj = crumblr(cell_counts)

bestModel=" ~ scale(PMI) + Sex + Source"
form = as.formula(bestModel)
fit = dream(cobj, form, metadata_df)
res_mat=residuals(fit)

bestModel=" ~ log2(Age + 1)"
form = as.formula(bestModel)
fit = dream(res_mat, form, metadata_df)
fit = eBayes(fit) 
head(fit$coefficients)

df_table=topTable(fit, coef='log2(Age + 1)', number=Inf) %>%   
  select(logFC, AveExpr, t, P.Value, adj.P.Val)
head(df_table)

hcl = buildClusterTreeFromPB(pb_subset, assays = grep("EN_L5_ET",assayNames(pb),invert=TRUE,value=TRUE))
cobj$E=res_mat
res1 = treeTest(fit, cobj, hcl, coef="log2(Age + 1)")
fig.tree = plotTreeTest(res1) + theme(legend.position="none")+xlim(0, 15)
tab_scaled_age = topTable(fit, "log2(Age + 1)", number=Inf, sort.by="none")
tab_scaled_age$celltype =factor(rownames(tab_scaled_age), rev(subclass_order))
tab_scaled_age$se = with(tab_scaled_age, logFC/ t)
fig.tree$data$label=factor(fig.tree$data$label,levels=rev(subclass_order))

pdf("/sc/arion/projects/psychAD/aging/kiran/analysis/crumblr/lifespan/wls_with_res.pdf")

fig.logFC = ggplot(tab_scaled_age, aes(celltype, logFC,color=celltype)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey50") + 
  geom_errorbar(aes(ymin=logFC - 1.96*se, ymax=logFC + 1.96*se), width=0) +
  geom_point(data=tab_scaled_age,aes(celltype, logFC,color=celltype,size=-log10(adj.P.Val))) +
  coord_flip() +
  theme_classic() +
  xlab('') +
  theme(aspect.ratio=3.65,axis.text.y = element_blank())+scale_color_manual(values=subclass_color_map)
 fig.logFC %>% insert_left(fig.tree)  #%>% insert_right(fig.vp)

plotCorrMatrix(cor(t(res_mat)))

dev.off()

write.csv(tab_scaled_age,file="/sc/arion/projects/psychAD/aging/kiran/analysis/crumblr/lifespan/wls_with_res.csv")


