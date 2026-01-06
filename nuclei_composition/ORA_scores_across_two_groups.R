.libPaths(c("/sc/arion/projects/psychAD/aging/kiran/RLib_4_4.1",.libPaths()))
library(zellkonverter)
library(basilisk)
library(dreamlet)
library(crumblr)
library(Seurat)
source("/sc/arion/work/girdhk01/scripts/chipseq_files.R")
source('/sc/arion/projects/roussp01a/pengfei/hicchip/scripts/hic_helper.R')
source("/sc/arion/work/girdhk01/scripts/myscripts/CRD/CMC_SV_help.R")
source("/sc/arion/work/girdhk01/scripts/myscripts/CRD/TRH_help_functions.R")
source("/sc/arion/work/girdhk01/scripts/myscripts/CRD/HiC_pengfei.R")
source("/sc/arion/work/girdhk01/scripts/myscripts/chiq_test.R")
source("/sc/arion/projects/psychAD/aging/kiran/final_model/plotting_functions.R")
library(SingleCellExperiment)

args <- commandArgs(trailingOnly = TRUE)
id=as.numeric(args[1]) - 1

# get metadata here
pb=readRDS("/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/pseudobulk/AGING_2024-02-01_22_23_PB_SubID_subclass.RDS")
colData(pb)$SubID=rownames(colData(pb))
colData(pb)$groups="NA"
colData(pb)$groups[colData(pb)$Age<20]="Childhood"
colData(pb)$groups[colData(pb)$Age>=20 & colData(pb)$Age<40]="Adulthood"
colData(pb)$groups[colData(pb)$Age>=40 & colData(pb)$Age<60]="Adulthood"
colData(pb)$groups[colData(pb)$Age>=60]="Adulthood"
metadata=as.data.frame(colData(pb)[,c("SubID","Sex","PMI","Brain_bank","groups")])

### get nuclei metadata here

nuclei_metadata=fread('/sc/arion/projects/psychAD/aging/kiran/analysis/crumblr/ORA/Aging_nuclei_metadata.csv')
df=as.data.frame(nuclei_metadata)
nuclei_metadata_info=merge(df,metadata,by="SubID")
nuclei_metadata_info=add_class_metadata(nuclei_metadata_info,"subclass")
nuclei_metadata_info$class <- as.character(nuclei_metadata_info$class)


# extract all pathway scores here
pathway_scores=fread(paste0("/sc/arion/projects/psychAD/aging/kiran/analysis/crumblr/ORA/ORA_estimates_by_chunks/ORA_chunk_",id,".csv"))
pathway_scores=as.data.frame(pathway_scores)
rownames(pathway_scores)=pathway_scores$V1

pathway_scores_subset=pathway_scores[,-1]

nuclei_metadata_info=nuclei_metadata_info#[nuclei_metadata_info$class %in% c("EN", "IN"), ]
nuclei_metadata_info=nuclei_metadata_info[nuclei_metadata_info$subclass!="EN_L5_ET",]
keep_nuclei=intersect(rownames(pathway_scores_subset),nuclei_metadata_info$V1)

nuclei_metadata_subset=nuclei_metadata_info[match(keep_nuclei,nuclei_metadata_info$V1),]
pathway_scores_subset=pathway_scores_subset[match(keep_nuclei,rownames(pathway_scores_subset)),]

nuclei_metadata_subset=as.data.frame(nuclei_metadata_subset)


run_cluster_dream <- function(subset_meta, pathway_scores_subset, class_label,group_label) {
  
	subset_meta$scaled_Age <- scale(subset_meta$Age)
  rownames(subset_meta) <- subset_meta$V1
    
  path_scores <- pathway_scores_subset[match(subset_meta$V1, rownames(pathway_scores_subset)), ]
  path_scores_matrix <- t(path_scores)
      
  stopifnot(identical(colnames(path_scores_matrix), subset_meta$V1))
  
	if (group_label=="Childhood") {
      	formula <- ~ 0 + groups + Sex + scaled_Age + scale(PMI) + (1|SubID)
	} else {
      	formula <- ~ 0 + groups + Sex + scaled_Age + scale(PMI) + Brain_bank + (1|SubID)
	}


L <- makeContrastsDream(formula, subset_meta,
    contrasts = c(
    CH_YA = "groupsChildhood - groupsAdulthood"
  )
)
  

  fit <- dream(path_scores_matrix, formula, subset_meta, L)
  fit <- eBayes(fit)
  # Extract results for each contrast
  results_list <- list()
  
  for (contrast_name in colnames(L)) {
    res <- topTable(fit, coef = contrast_name, number = Inf)
    res$Contrast <- contrast_name
    res$Pathway <- rownames(res)
    res$Group <- group_label
    res$class <- class_label
    results_list[[contrast_name]] <- res
  }

  return(do.call(rbind, results_list))

}



class_groups <- unique(nuclei_metadata_subset$class)

results_list <- list()
ct <- 1

# List all subclasses
all_subclasses <- unique(nuclei_metadata_subset$subclass)

# ---- Subclass-wise analysis ----
for (subclass_label in all_subclasses) {
  subset_subclass <- nuclei_metadata_subset[nuclei_metadata_subset$subclass == subclass_label, ]
  if (nrow(subset_subclass) > 10) {
    res <- run_cluster_dream(subset_subclass, pathway_scores_subset, subclass_label, "All_groups")
    if (!is.null(res)) {
      res$ncells <- nrow(subset_subclass)
      results_list[[ct]] <- res
      ct <- ct + 1
    }
  }
}


# Combine all results
final_df <- do.call(rbind, results_list)


# Save the results
save(final_df, file = paste0("/sc/arion/projects/psychAD/aging/kiran/analysis/crumblr/ORA/dream_wls_across_two_groups_", id, ".RDATA"))




