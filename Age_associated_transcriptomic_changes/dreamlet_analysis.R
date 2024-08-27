devtools::install_github("DiseaseNeurogenomics/dreamlet")

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
subclass_color_map= c("#659EC7","#95B9C7","#6495ED","#79BAEC","#0020C2","#4863A0","#3BB9FF","#B4CFEC","#1589FF","#488AC7","#4E8975","#7BCCB5","#3B9C9C","#B2C248","#008080",
    "#728C00","#89C35C","#F75D59","#C11B17","#6F4E37","#E0B0FF","#4E387E","#B93B8F","#FFA62F","#C19A6B", "#ECE5B6","#FFF380","gray")
names(subclass_color_map)=c("EN_L2_3_IT","EN_L3_5_IT_1","EN_L3_5_IT_2","EN_L3_5_IT_3","EN_L5_6_NP","EN_L5_ET","EN_L6_CT","EN_L6_IT_1","EN_L6_IT_2",
                            "EN_L6B","IN_ADARB2","IN_LAMP5_LHX6","IN_LAMP5_RELN","IN_PVALB","IN_PVALB_CHC","IN_SST","IN_VIP","Micro","PVM","Adaptive","PC","SMC","VLMC","Endo","Astro","Oligo","OPC","NS")

groups_list=c("Developmental","Young_Adulthood","Middle_Adulthood","Late_Adulthood")
colors_groups_list=c("#F9CFA1","#EE9B00","#CA6702","#8C510A")
names(colors_groups_list)=groups_list

good_plot=function(g1,sz,rot){
g1=g1+theme_classic()+theme(axis.title= element_text(color="black", size=sz),axis.text = element_text(color="black", size=sz),axis.text.x=element_text(angle=rot,hjust=1,vjust=1))
g1=g1+theme(legend.position="right")+theme(strip.text=element_text(face="bold",color="black"),strip.background=element_rect(fill="#eeeeee",color="#eeeeee"))
g1+theme(axis.title= element_text(color="black", size=sz),axis.text = element_text(color="black", size=sz))+theme(aspect.ratio=1)+theme(panel.border = element_rect(size = 1, color = "black", fill = NA))
}

add_class2_metadata=function(metadata,col_name){
class2_order=c("EN","IN","Non_neurons","Adaptive_Endo_Mural")
metadata$class="Endo_Mural"
metadata$class[grep("^EN",metadata[,col_name])]="EN"
metadata$class[grep("^IN",metadata[,col_name])]="IN"
metadata$class[grep("Astro",metadata[,col_name])]="Non_neurons"
metadata$class[grep("OPC",metadata[,col_name])]="Non_neurons"
metadata$class[grep("Olig",metadata[,col_name])]="Non_neurons"
metadata$class[grep("Micro",metadata[,col_name])]="Non_neurons"
metadata$class[grep("Endo",metadata[,col_name])]="Adaptive_Endo_Mural"
metadata$class[grep("^PC|SMC|VLMC|Adaptive|PVM",metadata[,col_name])]="Adaptive_Endo_Mural"
metadata$class=factor(metadata$class,levels=(class2_order))
metadata

}
## code starts here 
data_dir="age_associated_transcriptomic_changes" 
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

### estimate coefficient of age from dreamlet analysis and runs mashr analysis using the summarys stats from dreamlet for four age groups

for (kk in (1:length(groups_list))) {
  
pb_subset_final=pb[,colData(pb)$groups==groups_list[kk]]

if (input_groups == "Developmental") {

form=" ~ scaled_age + Sex + scale(PMI) + log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes"
form=as.formula(form)

} else {

form="~ scaled_age + Source + Sex + scale(PMI) + log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes"
form=as.formula(form)

}

res.proc = processAssays(pb_subset_final,as.formula(form),assays=grep("EN_L5_ET",assayNames(pb_subset_final),value=TRUE,invert=TRUE))

## estimates coefficient of all variables using dreamlet per subclass 
res.dl = dreamlet(res.proc, form)

## mashr object for shared genes analysis
res_mash=run_mash(res.dl, coef="scaled_age") 

save(res_mash,res.dl,form,file=paste0(data_dir,"/",groups_list[kk],".RDATA"))

}

### collect dreamlet summary stats for each age group

setwd("age_associated_transcriptomic_changes")
files=list.files(pattern="*RDATA")

four_groups_specific=list()

for (ct in (1:length(groups_list))){
load(grep(groups_list[ct],files,value=TRUE))
four_groups_specific[[ct]]=topTable(res.dl, coef = "scaled_age",num=Inf)
four_groups_specific[[ct]]$age_groups=groups_list[ct]
}

four_groups_specific_summary_df=do.call(rbind,four_groups_specific)
save(four_groups_specific_summary_df,file=paste0(data_dir,"/four_groups_specific_dreamlet_summary.RDATA"))

# Supplementary Figure 3
four_groups_specific_summary_df$age_groups=factor(four_groups_specific_summary_df$age_groups,levels=groups_list)
four_groups_specific_summary_df$assay=factor(four_groups_specific_summary_df$assay,levels=subclass_order)
four_groups_specific_summary_df$plot="NS"
four_groups_specific_summary_df$plot[four_groups_specific_summary_df$adj.P.Val<.05]=as.character(four_groups_specific_summary_df$assay[four_groups_specific_summary_df$adj.P.Val<.05])
four_groups_specific_summary_df$type=1
four_groups_specific_summary_df$type[four_groups_specific_summary_df$adj.P.Val<.05]=2


png("all_groups_aging_effect_sizes.png", width = 4000, height = 4000, res = 300)
g1=ggplot(subset(four_groups_specific_summary_df),aes(x=factor(assay), y=logFC,size=factor(type),color=factor(plot)))+ geom_jitter()+scale_color_manual(values=subclass_color_map)+facet_wrap(~age_groups,ncol=1)+theme_bw()
good_plot(g1,12,45)+theme(legend.position="top")+theme(aspect.ratio=0.15)+scale_size_manual(values=c(0.15,2))+ylab("Age associated changes")+xlab("subclass")
dev.off()

# Extended Figure 3a

four_groups_cell_specific_summary=as.data.frame(four_groups_specific_summary_df %>% dplyr::group_by(celltype,age_groups) %>% dplyr::summarize(n_DGE_genes=uniqueN(ID[adj.P.Val<.05]),tot_genes=uniqueN(ID)))
four_groups_cell_specific_melted_summary=melt(four_groups_cell_specific_summary[,c("celltype","age_groups","n_DGE_genes","tot_genes")],id=c("celltype","age_groups","n_DGE_genes","tot_genes"))
four_groups_cell_specific_melted_summary$age_groups=factor(four_groups_cell_specific_melted_summary$age_groups,levels=groups_list)
four_groups_cell_specific_melted_summary$celltype=factor(four_groups_cell_specific_melted_summary$celltype,levels=subclass_order)
four_groups_cell_specific_melted_summary=add_class2_metadata(four_groups_cell_specific_melted_summary,"celltype")

df=four_groups_cell_specific_melted_summary %>% group_by(age_groups) %>% summarize(nDGE=sum(n_DGE_genes))
df$age_groups=factor(df$age_groups,levels=groups_list)
df=as.data.frame(df)

pdf(paste0(data_dir,"/all_groups_summary_bar_plot.pdf"))

df_plot %>%
  ggplot(aes(x=age_groups, y=logFC)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(data=df_plot,aes(color=factor(celltype)), size=2, alpha=0.9) +scale_color_manual(values=subclass_color_map)+
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle(gene) +
    xlab("")

dev.off()








