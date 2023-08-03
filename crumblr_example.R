library(dreamlet)
library(crumblr)

## Example code for "no duplicates" matrices

pb=readRDS("/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/pseudobulk/AGING_2023-06-09_01_45_PB_SubID_subtype.RDS")
## this is your normalized cell counts matrices
cobj = crumblr(cellCounts(pb))
cobj_df=as.data.frame(t(as.matrix(cobj))) ### matrix
cobj_df$samples=rownames(cobj_df)

### this is the step to prepare covariates  whicha are Age, PMI, UMI counts, brain bank and sex
colData(pb)$scaled_Age=as.numeric(scale(colData(pb)$Age))
colData(pb)$scaled_PMI=as.numeric(scale(colData(pb)$PMI))
metadata(pb)$aggr_means$log_n_counts=as.numeric(log(metadata(pb)$aggr_means$n_counts))
metadata=as.data.frame(colData(pb))
metadata$samples=rownames(metadata)

df=merge(cobj_df,metadata[,c("Age","samples")],by="samples")
df=melt(df,id=c("samples","Age"))
df$celltypes=df$variable

### here we run cell composition test;  you can add TOD here 
## soucr is brain bank
bestModel=" ~ 0 + TOD + scaled_Age + (1 | Source) +  (1 | Sex) + scaled_PMI + log_n_counts"
form = as.formula(bestModel)
fit = dream(cobj, form, METADATA_aggr_metadata_df)
fit = eBayes(fit)
colnames(coef(fit))

library(ggtree)
library(aplot)

hc = buildClusterTreeFromPB(pb)
res = treeTest( fit, cobj, hc, coef="scaled_Age")
fig1 = plotTreeTest(res) + xlim(0, 15) + theme(legend.position="none") + ggtitle("scaled_Age")

## here you can extract the coefficient you want
tab = topTable(fit, coef="TOD", number=Inf)
tab$celltype = factor(rownames(tab), rev(get_taxa_name(fig1)))
tab$se = with(tab, logFC/t)
res=as_tibble(res)

fig2 = ggplot(tab, aes(celltype, logFC)) + 
  geom_hline(yintercept=0, linetype="dashed", color="grey", size=1) +
  geom_errorbar(aes(ymin = logFC - 1.96*se, ymax = logFC + 1.96*se), width=0) +
  geom_point(color="dodgerblue") +
  theme_classic() +
  coord_flip() +
  xlab('') + 
  ylab("Effect size") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# combine plots
fig2 %>% insert_left(fig1) 
tab$mylabel=""
tab$mylabel[tab$P.Value<.05]="*"
tab$mylabel[tab$adj.P.Val<.05]="#"

 fig2 = ggplot(tab, aes(celltype, logFC,label=mylabel)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey", size=1) +
  geom_errorbar(aes(ymin = logFC - 1.96*se, ymax = logFC + 1.96*se), width=0) +
  geom_point(color="dodgerblue") + geom_text()+
  theme_classic() +
  coord_flip() +
  xlab('') + 
  ylab("Effect size")



