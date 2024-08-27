
data_dir="degree_of_sharing" 
syn62064718=synGet(entity="syn62064718",downloadLocation=data_dir)
pb=readRDS(paste0(data_dir,"/lifespan_pseudobulk.rds"))
colData(pb)$SubID=rownames(colData(pb))

## age groups assignment
colData(pb)$groups="NA"
colData(pb)$groups[colData(pb)$Age<20]="Developmental"
colData(pb)$groups[colData(pb)$Age>=20 & colData(pb)$Age<40]="Young_Adulthood"
colData(pb)$groups[colData(pb)$Age>=40 & colData(pb)$Age<60]="Middle_Adulthood"
colData(pb)$groups[colData(pb)$Age>=60]="Late_Adulthood"


pb_subset=pb[,colData(pb)$groups=="Developmental"]

# Stack EN assays to estimate Tau scores 
pb.stack=stackAssays(pb, assays = grep("^EN",subclass_order,invert=FALSE,value=TRUE))
geneExpr <- assay(pb.stack, 1)
rownames(geneExpr) <- rownames(pb.stack)
geneExpr_keep=geneExpr[rownames(geneExpr) %in% PC_genes$gene_name,]

### convert counts to cpm
idx <- which(colSums(geneExpr_keep, useNames = TRUE) <5)
if (length(idx)>0){
geneExpr_keep <- geneExpr_keep[, -idx]
}

y <- DGEList(geneExpr_keep, remove.zeros = TRUE)
y <- calcNormFactors(y, method = "TMM")
geneExpr_cluster_cpm <- edgeR::cpm(y, log = FALSE)
geneExpr_cluster_cpm=as.data.frame(geneExpr_cluster_cpm)
geneExpr_cluster_cpm$genes=rownames(geneExpr_cluster_cpm)

geneExpr_cluster_cpm_by_donors <- geneExpr_cluster_cpm %>%
    tidyr::gather(celltype_ID, exp, colnames(select_if(.,is.numeric))) %>% mutate(celltypes=sub("_[^_]+$", "", celltype_ID)) %>%
    group_by(genes,celltypes) %>% summarize(median_exp=median(exp))

geneExpr_cluster_cpm_by_donors_mat=geneExpr_cluster_cpm_by_donors %>% spread(celltypes, median_exp)
geneExpr_cluster_cpm_by_donors_mat=as.data.frame(geneExpr_cluster_cpm_by_donors_mat)

#### This step is for finding cell specific scores

num_of_cells <- geneExpr_cluster_cpm_by_donors_mat %>% select_if(is.numeric) %>% names() %>% length()
result <- geneExpr_cluster_cpm_by_donors_mat %>%
    tidyr::gather(celltypes, exp, colnames(select_if(.,is.numeric))) %>%   # gather only numeric columns
    mutate_if(is.numeric, list(~ ifelse(. < 1, 1, .))) %>%  # if there's any expression less than 1, make it 1
    group_by(genes) %>%                             # calculations will be performed for a gene up to ungroup()
    mutate(lessthan1 = sum(exp <= 1),
           lessthan10 = sum(exp <= 10))%>%
    mutate(log2exp = log2(exp),                             # log transformation is needed for tau calculation
           sum_log = sum(log2exp),
           max_log = max(log2exp),
           tau = (num_of_cells - sum_log/max_log) / (num_of_cells -1))  %>%
    ungroup() %>%
    mutate(
      Status = case_when(
        lessthan1 == num_of_cells  ~ "Null expresssion",
        tau < 0.85                   ~ "Wide-spread expression",
        lessthan10 == num_of_cells ~ "Weak expresssion",
        TRUE                         ~ "Specific expression"
      )
    ) %>%     # distinct(EnsemblGeneID,Status) %>% count(Status) : status counts, see below
    mutate(tau = ifelse(Status == "Null expresssion", -1, tau)) %>%      # for null exp genes, assign tau as -1
    select(-(lessthan1:max_log), Tau.Score=tau)  %>%   # remove unnecessary columns before spread
    spread(celltypes, exp) %>%
    select(genes, Tau.Score, Status) %>%
    left_join(geneExpr_cluster_cpm_by_donors_mat, ., by="genes")  # had to use dot since new columns are wanted at the end

save(result,file=paste0(data_dir,"/tauscore_stacked_EN_data_developmental.RDATA"))





