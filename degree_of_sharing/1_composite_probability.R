library(mashr)
library(dreamlet)

dir.create("degree_of_sharing")

groups_list=c("Developmental","Young_Adulthood","Middle_Adulthood","Late_Adulthood")

## grouped classes
glias <- c("Astro", "Micro","Oligo","OPC")
EN <- c("EN_L6_CT", "EN_L5_6_NP", "EN_L6B", "EN_L3_5_IT_1", "EN_L3_5_IT_2", "EN_L3_5_IT_3","EN_L2_3_IT","EN_L6_IT_1",
    "EN_L6_IT_2")
IN <- c("IN_LAMP5_RELN", "IN_LAMP5_LHX6","IN_ADARB2","IN_PVALB_CHC", "IN_PVALB", "IN_SST", "IN_VIP")
cell_class_list=list("EN"=EN,"IN"=IN,"glias"=glias)

### get mashr results from age_associated_transcriptomic_changes directory

setwd("age_associated_transcriptomic_changes")
files=list.files(pattern="*RDATA")

combined_prob=list()

for (group_id in (1:length(groups_list))){

load(grep(groups_list[group_id],files,value=TRUE))

prob_celltypes=list()

for (i in (1:length(cell_class_list))){
prob_celltypes[[i]] <- compositePosteriorTest(res_mash, include=cell_class_list[[i]],exclude=NULL,test="all")

# compositePosteriorTest function can be found here  https://github.com/GabrielHoffman/dreamlet/blob/devel/R/compositePosteriorTest.R 
# this function gives P_E, P_I and P_G
  
prob_celltypes[[i]]=as.data.frame(prob_celltypes[[i]])
colnames(prob_celltypes[[i]])="prob"
prob_celltypes[[i]]$genes=rownames(prob_celltypes[[i]])
prob_celltypes[[i]]$cell_class=names(cell_class_list)[i]

}

combined_prob[[group_id]] <- do.call(rbind, prob_celltypes)
combined_prob[[group_id]]$group=groups_list[[group_id]]
combined_prob[[group_id]]$group_cell_class_gene=paste0(combined_prob[[group_id]]$group,"_",combined_prob[[group_id]]$cell_class,"_",combined_prob[[group_id]]$genes)

}

## gather number of celltypes for each gene with age associated p value < .05 from dreamlet analysis 

load("age_associated_transcriptomic_changes/four_groups_specific_dreamlet_summary.RDATA")
genes_significance_cells_counts=list()
ct = 1

for (ii in (1:length(groups_list))){

for (i in (1:length(cell_class_list))){

temp_df=four_groups_specific_summary_df[four_groups_specific_summary_df$assay %in% cell_class_list[[i]] & four_groups_specific_summary_df$group==groups_list[ii],]
counts_cells_for_signficant_genes=temp_df %>% group_by(ID) %>% summarize(ncells=uniqueN(assay[P.Value<.05]),ncells_aDEG=uniqueN(assay[adj.P.Val<.05]))
genes_significance_cells_counts[[ct]]=counts_cells_for_signficant_genes
genes_significance_cells_counts[[ct]]$cell_class=names(cell_class_list)[i]
genes_significance_cells_counts[[ct]]$age_groups=groups_list[ii]
genes_significance_cells_counts[[ct]]$group_cell_class_gene=paste0(genes_significance_cells_counts[[ct]]$group,"_",genes_significance_cells_counts[[ct]]$cell_class,"_",genes_significance_cells_counts[[ct]]$ID)

ct = ct + 1

}

}
genes_significance_cells_counts_df=do.call(rbind,genes_significance_cells_counts)

combined_prob_age_groups_df <- do.call(rbind,combined_prob)
combined_prob_age_groups_ncells_df <- merge(combined_prob_age_groups_df, genes_significance_cells_counts_df[,c("group_cell_class_gene","ncells","ncells_aDEG","ncells_aDEG_per_class")], by = "group_cell_class_gene", all.x = TRUE)
combined_prob_age_groups_ncells_df$ncells[is.na(combined_prob_age_groups_ncells_df$ncells)] <- 0
combined_prob_age_groups_ncells_df$ncells_aDEG[is.na(combined_prob_age_groups_ncells_df$ncells_aDEG)] <- 0




