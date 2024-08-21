# monocle3 conda environment
# Load libraries
library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(zellkonverter)
set.seed(222)
library(RColorBrewer)
library(viridis)
library(monocle3)
library(tidyselect)
library(grid)
library(mgcv)
library(colorspace)
library(ggrepel)
library(igraph)
library(pbapply)
library(devtools)
library(parallel)
library(evobiR)
library(tidyr)
library(cluster)
library(grDevices)
library(repr)
library(zoo)
library(ggnewscale)
library(VennDiagram)

setClass("cell_data_set_ext", contains = "cell_data_set", slots=c(graphs = "list", lineages="list", expression="list", expectation="list", pseudotime="list")) -> cell_data_set_ext

import_monocle <-function(cds){
  cds <- as(cds,"cell_data_set_ext")
  return(cds)
}

pseudotime_infer_func <- function(file_lin, file_out){
    d <- readH5AD(file_lin)
    d.seurat <- as.Seurat(d, counts = "raw-cts_pre-ds", data = "ds_norm_cts")

    # prepare cell_data_set object
    o0 <- CreateSeuratObject(counts = d.seurat[["originalexp"]]@counts, meta.data = d.seurat@meta.data)
    o0 <- SetAssayData(o0, slot = "data", new.data = d.seurat[["originalexp"]]@data)
    o0 <- SetAssayData(o0, slot = "scale.data", new.data = d.seurat[["originalexp"]]@scale.data)
    o0[["UMAP"]] <- CreateDimReducObject(embeddings = Embeddings(d.seurat, "X_umat"), key = "UMAP")
    cds <- SeuratWrappers::as.cell_data_set(o0)
    head(colData(cds))
    cds <- estimate_size_factors(cds)
    head(colData(cds))
    
    Idents(d.seurat) <- d.seurat@meta.data$subtype_uni
    s.clusters <- as.character(Idents(d.seurat))
    names(s.clusters) <- names(Idents(d.seurat))
    s.clusters <- s.clusters[colnames(cds)]
    
    cds@clusters$"UMAP"$"clusters" <- s.clusters
    cds@clusters$UMAP$partitions <- cds@clusters$UMAP$clusters
    cds@clusters$UMAP$partitions[cds@clusters$UMAP$partitions != "1"] <- "1"
    
    # construct trajectory graph
    graph_control <- setNames(list(3, 1, 1, FALSE, TRUE, FALSE, NULL, 10, 1e-5, 0.005, 0.01), c("euclidean_distance_ratio","geodesic_distance_ratio","minimal_branch_len","orthogonal_proj_tip","prune_graph","scale","rann.k","maxiter","eps","L1.gamma","L1.sigma"))
    cds <- learn_graph(cds, use_partition = F, learn_graph_control = graph_control)
    
    #import the monocle object to add new slots
    cds <- import_monocle(cds)
    
    #connect nodes for trajectories not joined in the original graph
    #plot nodes
    Y <- cds@principal_graph_aux[["UMAP"]]$dp_mst
    d <- as.data.frame(t(Y))
    
    rowData(cds)$gene_name <- rownames(cds)
    rowData(cds)$gene_short_name <- rowData(cds)$gene_name
    saveRDS(cds, file = paste0(file_out, "_beforeConstruction_cds.RDS"))
    saveRDS(d, file = paste0(file_out, "_beforeConstruction_d.RDS"))
}

args <- commandArgs(T)
file_lin <- args[1] # 
file_out <- args[2] # 

pseudotime_infer_func(file_lin, file_out)