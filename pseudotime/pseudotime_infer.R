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

source("0_color_settings.r")
setClass("cell_data_set_ext", contains = "cell_data_set", slots=c(graphs = "list", lineages="list", expression="list", expectation="list", pseudotime="list")) -> cell_data_set_ext

import_monocle <-function(cds){
  cds <- as(cds,"cell_data_set_ext")
  return(cds)
}

pseudotime_infer_func <- function(file_lin, file_out, fig_width, fig_height){
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

    # Visualization
    p_subclass_uni <- plot_cells(cds,
               reduction_method = 'UMAP',
               color_cells_by = 'subclass_uni',
               cell_size = 0.8,
               alpha = 0.6,
               label_branch_points = T,
               label_leaves = F,
               label_roots = F,
               group_label_size = 4,label_groups_by_cluster = F) + 
               scale_color_manual(values = cols_subclass_uni) + 
               theme(legend.text=element_text(size=8), plot.title = element_text(hjust = 0.5)) +
               ggtitle('subclass_uni')

    p_subtype_uni <- plot_cells(cds,
               reduction_method = 'UMAP',
               color_cells_by = 'subtype_uni',
               cell_size = 0.8,
               alpha = 0.6,
               label_branch_points = T,
               label_leaves = F,
               label_roots = F,
               group_label_size = 4,label_groups_by_cluster = F) + 
               # scale_color_manual(values = cols_subtype_uni) + 
               theme(legend.text=element_text(size=8), plot.title = element_text(hjust = 0.5)) +
               ggtitle('subtype_uni')

    p_sex <- plot_cells(cds,
               reduction_method = 'UMAP',
               color_cells_by = 'Sex',
               cell_size = 0.8,
               alpha = 0.6,
               label_branch_points = T,
               label_leaves = F,
               label_roots = F,
               group_label_size = 1,label_groups_by_cluster = F) + 
               theme(legend.position = 'top', legend.text=element_text(size=8), plot.title = element_text(hjust = 0.5)) +
               ggtitle('Sex') +
               scale_fill_manual(values = colorRampPalette(brewer.pal('Dark2', n = 8))(2))

    p_stage_id <- plot_cells(cds,
               reduction_method = 'UMAP',
               color_cells_by = 'stage_id',
               cell_size = 0.8,
               alpha = 0.6,
               label_branch_points = T,
               label_leaves = F,
               label_roots = F,
               group_label_size = 4,label_groups_by_cluster = F) + 
               scale_color_manual(values = cols_stage_id) + 
               theme(legend.text=element_text(size=8), plot.title = element_text(hjust = 0.5)) +
               ggtitle('stage_id')

    p_age <- plot_cells(cds,
               reduction_method = 'UMAP',
               color_cells_by = 'Age',
               cell_size = 0.8,
               alpha = 0.6,
               label_branch_points = T,
               label_leaves = F,
               label_roots = F,
               group_label_size = 1,label_groups_by_cluster = F) + 
               scale_color_viridis(direction = -1) +
               theme(legend.text=element_text(size=8), plot.title = element_text(hjust = 0.5)) +
               ggtitle('Age')

    pdf(paste0(file_out, '_visualization.pdf'), width = fig_width, height = fig_height)
        pushViewport(viewport(layout = grid.layout(2, 3)))
        vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
        print(p_subclass_uni, vp = vplayout(1, 1))
        print(p_subtype_uni, vp = vplayout(1, 2))
        print(p_sex, vp = vplayout(1, 3))
        print(p_stage_id, vp = vplayout(2, 1))
        print(p_age, vp = vplayout(2, 2))
    dev.off()
}

args <- commandArgs(T)
file_lin <- args[1] # 
file_out <- args[2] # 
fig_width <- as.numeric(args[3]) # 
fig_height <- as.numeric(args[4]) # 

pseudotime_infer_func(file_lin, file_out, fig_width, fig_height)