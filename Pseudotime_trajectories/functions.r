source("initial_settings.r")

branch_lineage <- c("L6B" = "Deep-non-IT", "L6_CT" = "Deep-non-IT", "L5_6_NP" = "Deep-non-IT", 
                    "L6_IT_1" = "Deep-IT", "L6_IT_2" = "Deep-IT", 
                    "L2_3_IT" = "Upper-IT", "L3_5_IT_1" = "Upper-IT", "L3_5_IT_2" = "Upper-IT", "L3_5_IT_3" = "Upper-IT", 
                    "SST" = "MGE", "PVALB" = "MGE", "PVALB_CHC" = "MGE", 
                    "VIP_TRPC6" = "CGE", "VIP_BCL11B" = "CGE", "ADARB2_RAB37" = "CGE", "ADARB2_COL12A1" = "CGE", "ADARB2_SYT10" = "CGE", "ADARB2_SV2C" = "CGE", "LAMP5_RELN" = "CGE", "LAMP5_LHX6" = "CGE")


regress_pt_along_age_trajectory <- function(cds, lineage){
    par(pty = "s")
    plot(NULL, xlab = "Age", ylab = "Pseudotime", main = paste0("Regression of pseudotime along age in ", lineage), 
         xlim = c(min(colData(cds)$numerical_age), max(colData(cds)$numerical_age)), 
         ylim = c(min(pseudotime(cds)), max(pseudotime(cds))), type="n")
    for(sub_lin in names(cds@lineages)){
        mod.smsp <- smooth.spline(x = cds@colData['numerical_age'][cds@lineages[[sub_lin]], ], 
                                  y = pseudotime(cds)[cds@lineages[[sub_lin]]], 
                                  spar = 1)
        lines(mod.smsp$x, mod.smsp$y, lwd = 3, col = cols_traj[sub_lin])
    }
    legend("bottomright", legend = names(cds@lineages), col = cols_traj[names(cds@lineages)], lwd = 3, inset = c(0.02, 0.02))
}

regress_pt_along_age_branch <- function(cds, lineage, branches){
    par(pty = "s")
    plot(NULL, xlab = "Age", ylab = "Pseudotime", main = paste0("Regression of pseudotime along age in ", lineage), 
         xlim = c(min(colData(cds)$numerical_age), max(colData(cds)$numerical_age)), 
         ylim = c(min(pseudotime(cds)), max(pseudotime(cds))), type="n")
    for(sub_branch in branches){
        age <- c()
        pseudo <- c()

        for(sub_lin in names(branch_lineage[branch_lineage == sub_branch])){
            age <- c(age, cds@colData['numerical_age'][cds@lineages[[sub_lin]], ])
            pseudo <- c(pseudo, pseudotime(cds)[cds@lineages[[sub_lin]]])
        }
        
        mod.smsp <- smooth.spline(x = age, 
                                  y = pseudo, 
                                  spar = 1)
        lines(mod.smsp$x, mod.smsp$y, lwd = 3, col = cols_traj[sub_branch])
    }
    legend("bottomright", legend = branches, col = cols_traj[branches], lwd = 3, inset = c(0.02, 0.02))
}


# Expression dynamics of traDEGs
scale_expression <- function(expr) {
    (expr - min(expr, na.rm = TRUE)) / 
    (max(expr, na.rm = TRUE) - min(expr, na.rm = TRUE))
}

adjust_fitted_matrix <- function(cds){
    fitted <- cds@expectation
    
    # Initialize a list to store the adjusted fitted matrices
    adjusted_fitted <- list()
    
    # Loop through each lineage
    for (lineage in names(fitted)) {
        fitted_matrix <- fitted[[lineage]]
        
        # If it's the first lineage, calculate the reference fitted for the first meta-cell
        if (lineage == names(fitted)[1]) {
            reference_fitted <- fitted_matrix[1, ]
        }
        
        # Calculate the difference for each gene
        adjustment <- reference_fitted - fitted_matrix[1, ]
        
        # Adjust the fitted matrix for the lineage
        adjusted_fitted_matrix <- sweep(fitted_matrix, 2, adjustment, '+')
        
        # Add the adjusted matrix to the list
        adjusted_fitted[[lineage]] <- adjusted_fitted_matrix
    }
    cds@expectation <- adjusted_fitted
    return(cds)
}

# Function to order rows within a cluster
order_rows_within_cluster <- function(cluster_data) {
    # Calculate row means (excluding the cluster column)
    row_means <- rowMeans(cluster_data[, -ncol(cluster_data)])
    # Order rows based on row means
    ordered_fitted_mat <- cluster_data[order(row_means), ]
    return(ordered_fitted_mat)
}

obtain_peak_times <- function(df, n_metacell){
    n_groups <- ncol(df) / n_metacell
    peak_times <- matrix(NA, nrow(df), n_groups)
    
    for (i in 1:n_groups) {
      cols <- ((i - 1) * 500 + 1):(i * 500)
    
      # Find the index of the max value in each row for these columns
      peak_times[, i] <- apply(df[, cols], 1, which.max)
    }
    
    # Optional: Adjust indices if necessary (e.g., to reflect actual times)
    
    # Calculate the average peak time for each row
    row_averages <- rowMeans(peak_times, na.rm = TRUE)
    return(row_averages)
}

# Source: https://gist.github.com/mathzero/a2070a24a6b418740c44a5c023f5c01e
save_pheatmap <- function(x, filename, width=12, height=12){
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    if(grepl(".png",filename)){
        png(filename, width=width, height=height, units = "in", res=300)
        grid::grid.newpage()
        grid::grid.draw(x$gtable)
        dev.off()
    }
    else if(grepl(".pdf",filename)){
        pdf(filename, width=width, height=height)
        grid::grid.newpage()
        grid::grid.draw(x$gtable)
        dev.off()
    }
    else{
      print("Filename did not contain '.png' or '.pdf'")
    }
}


cluster_and_func_enrich <- function(fitted_mat, n_tree, pseudo_range, lineages, branches, n_metacell, n_term, n_char, out_label, plotTextSize=11, width=12, height=12){
    fitted_mat <- t(apply(fitted_mat, 1, scale_expression))

    
    ###### kmeans clustering
    # 1st
    set.seed(222)
    ks_1st <- kmeans(fitted_mat, centers = n_tree, iter.max = 1000)

    
    # prepare for 2nd clustering
    if(length(pseudo_range) == 2){
        peak_times <- abs(apply(fitted_mat, 1, which.max) - n_metacell)
    }else{
        peak_times <- obtain_peak_times(fitted_mat, n_metacell)
    }
    
    # assign each gene to its cluster and add peak times
    gene_info <- data.frame(Cluster = ks_1st$cluster, PeakTime = peak_times)
    # calculate mean peak expression time for each cluster
    cluster_peak_times <- aggregate(PeakTime ~ Cluster, data = gene_info, mean)
    # order clusters based on mean peak expression times
    ordered_clusters <- cluster_peak_times[order(cluster_peak_times$PeakTime), "Cluster"]
    
    
    # 2nd: kmeans secondly according to specific order
    set.seed(222)
    centerOrder <- ks_1st$centers[ordered_clusters, ]
    ks_2nd <- kmeans(fitted_mat, centers = centerOrder, iter.max = 1000)
    # print(ks_2nd)
    annotation_row <- data.frame(Cluster = factor(ks_2nd$cluster))

    # Add cluster information to your fitted_mat  
    fitted_mat <- cbind(fitted_mat, data.frame(cluster = ks_2nd$cluster))

    # Apply the function to each cluster and combine the results
    ordered_fitted_mat <- do.call(rbind, by(fitted_mat, fitted_mat$cluster, order_rows_within_cluster))
    
    # Annotation for rows
    annotation_row <- data.frame(Cluster = factor(ordered_fitted_mat$cluster))
    rownames(annotation_row) <- rownames(ordered_fitted_mat)
    ordered_fitted_mat$cluster <- NULL

    # Annotation for columns
    annotation_col = data.frame(
                                CellType = factor(rep(lineages, each = n_metacell)), 
                                Branch = factor(branches))
    rownames(annotation_col) = colnames(ordered_fitted_mat)

    
    ###### Create the heatmap
    par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))
    cols_cluster <- c(brewer.pal(n = n_tree, name = "Paired"))
    names(cols_cluster) <- unique(annotation_row$Cluster)
    ann_colors <- list(CellType = cols_traj[lineages], Branch = cols_traj[unique(branches)], Cluster = cols_cluster)
    ph <- pheatmap(ordered_fitted_mat, 
             cluster_rows = F, 
             cluster_cols = F, 
             scale = "none", 
             annotation_row = annotation_row,
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             show_rownames = F, 
             show_colnames = F, 
             breaks = seq(0, 1, 0.01),
             color = cols_scale, 
             border_color = NA)
    save_pheatmap(ph, paste0(out_label, "_ks_", as.character(n_tree), "_heatmap.png"), width = width, height = height)

    ###### Plot average trend of each cluster
    # Set the layout for the plots - 1 row and n_tree columns
    pdf(paste0(out_label, "_ks_", as.character(n_tree), "_average_line.pdf"), width = 3 * n_tree, height = 3)
    # options(repr.plot.width = 3 * n_tree, repr.plot.height = 3, repr.plot.res = 300)
    par(mfrow = c(1, n_tree), mar = c(2, 1.8, 2, 1))
    
    for(sub_cluster in seq(1, n_tree, 1)){
        plot(0, type = "n", xlim = c(0, max(pseudo_range)), ylim = c(0, 1), xlab = "Pseudotime", ylab = "Relative expression", main = paste0("Cluster ", sub_cluster, " (", length(ks_2nd$cluster[ks_2nd$cluster == sub_cluster]), ")"))

        # Your lines plotting code
        if(length(pseudo_range) == 2){ 
            lines(seq(0, pseudo_range[1], length.out = n_metacell), 
                  rev(colMeans(fitted_mat[rownames(fitted_mat) %in% names(ks_2nd$cluster[ks_2nd$cluster == sub_cluster]), 1:n_metacell])), 
                  col = cols_traj[names(pseudo_range[1])], 
                  lwd = 10) 
        } else {
            lines(seq(0, pseudo_range[1], length.out = n_metacell), 
                  colMeans(fitted_mat[rownames(fitted_mat) %in% names(ks_2nd$cluster[ks_2nd$cluster == sub_cluster]), 1:n_metacell]), 
                  col = cols_traj[names(pseudo_range[1])], 
                  lwd = 10)  
        }
    
        if(length(pseudo_range) >= 2){
            for(i_branch in 2:length(pseudo_range)){
                lines(seq(0, pseudo_range[i_branch], length.out = n_metacell), 
                      colMeans(fitted_mat[rownames(fitted_mat) %in% names(ks_2nd$cluster[ks_2nd$cluster == sub_cluster]), n_metacell*(i_branch-1)+(1:n_metacell)]), 
                      col = cols_traj[names(pseudo_range)[i_branch]], 
                      lwd = 10)
            } 
        }
    }
    dev.off()

    # Reset default plotting layout
    par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))
    
    # Plot the legend in a new plot
    plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", main = "Legend")
    legend("topleft", legend = names(pseudo_range), col = cols_traj[names(pseudo_range)], lwd = 6, bty = "n")

    
    
    ###### Functional analysis
    options(repr.plot.width = 10, repr.plot.height = 7, repr.plot.res = 300)
    enriched_ggplot <- list()
    enriched_terms <- list()
    for(sub_cluster in seq(1, n_tree, 1)){
        GENES <- names(ks_2nd$cluster[ks_2nd$cluster == sub_cluster])
        write.table(unique(GENES), file = paste0(out_label, "_ks_", n_tree, "_cluster_", sub_cluster, ".csv"), sep = ",", row.names = F, col.names = F, quote = F)
        enriched <- enrichr(GENES, dbs)
        for(sub_enrich in names(enriched)){
            all_term <- enriched[[sub_enrich]]
            enriched_ggplot[[sub_enrich]] <- rbind(enriched_ggplot[[sub_enrich]], 
                                       data.frame(Cluster = rep(sub_cluster, n = nrow(all_term)), 
                                              Term = all_term$Term, 
                                              Count = sapply(strsplit(all_term$Genes, ";"), length), 
                                              P.value = all_term$P.value, 
                                              Adjusted.P.value = all_term$Adjusted.P.value,
                                              Neg.Log.Adj.P.Val = -log10(all_term$Adjusted.P.value)))
            top_term <- enriched[[sub_enrich]][order(enriched[[sub_enrich]]$Adjusted.P.value, decreasing = F)[1:n_term], ]
            top_term <- top_term[top_term$Adjusted.P.value < 0.05, ]
            if(nrow(top_term) > 0){
                enriched_terms[[sub_enrich]] <- unique(c(enriched_terms[[sub_enrich]], top_term$Term))               
            }
        }
    }
    
    for(sub_enrich in names(enriched_ggplot)){
        pdf(paste0(out_label, "_ks_", as.character(n_tree), "_func_enrich_", sub_enrich,".pdf"), width = 10, height = 7)
        enriched_data <- enriched_ggplot[[sub_enrich]]
        enriched_data$Cluster <- factor(enriched_data$Cluster, levels = seq(1, n_tree, 1))
        enriched_data <- enriched_data[enriched_data$Term %in% enriched_terms[[sub_enrich]], ]
        enriched_data$Term <- factor(enriched_data$Term, levels = enriched_terms[[sub_enrich]])
        # Add significance label
        enriched_data$plotLabel = ""
        # enriched_data[enriched_data$P.value < 0.05, "plotLabel"] <- "."
        enriched_data[enriched_data$Adjusted.P.value < 0.05, "plotLabel"] <- "#"
        sub_p <- ggplot(enriched_data, aes(x = Cluster, y = Term, size = Count, color = Neg.Log.Adj.P.Val)) +
            geom_point() +
            scale_size(range = c(5, 10)) +
            geom_text(aes(label = plotLabel), size = plotTextSize * 0.55, color = "black") +
            theme_minimal() +
            labs(title = sub_enrich, size = "Count", color = "Neg.Log.Adj.P.Val") +
            scale_color_gradient2(high="#D7191C", mid="white", low="#2C7BB6", name = "-log10(adj.p)")
        print(sub_p)
        dev.off()
    }
    return(list(ks = ks_2nd, fitted_mat = ordered_fitted_mat, enriched = enriched_ggplot))
}


# Plot LDsc / MAGMA heatmap
magma_heatmap <- function(df, plotCol="P", FDR_PER_TRAIT=F, myPadjustMethod="BH", markNominallySignificant=T,plotTextSize=11) {
    df = data.frame(df)
    df$sumstatName <- recode(df$sumstatName, 
                             "sz3"="SCZ",
                             "bip2"="BD",
                             "asd"="ASD",
                             "adhd_ipsych"="ADHD",
                             "mdd_ipsych"="MDD",
                             "ocd"="OCD",
                             "insomn2"="Insomnia",
                             "alcoholism_2019"="Alcoholism",
                             "tourette"="Tourettes",
                             "intel"="IQ",
                             "eduAttainment"="Education",
                             
                             "alzBellenguezNoApoe"="AD",
                             "ms"="MS",
                             "pd_without_23andMe"="PD",
                             "migraines_2021"="Migraines",
                             "als2021"="ALS",
                             "stroke"="Stroke",
                             "epilepsyFocal"="Epilepsy",
                                                        
                             "dm2"="T2D",
                             "ibd"="IBD",
                             "ra"="RArthritis",
                             "lipidCholTotal"="Cholesterol",
                             "obesity"="Obesity",
                             "uc"="UC")
    
    df$plotLabel = ""
    if(markNominallySignificant) {   df$plotLabel[(10^-df[,plotCol]) < 0.05] = "·" }
    df$plotLabel[p.adjust((10^-df[,plotCol]),method=myPadjustMethod) < 0.05]="#"
    
    if(FDR_PER_TRAIT) {
      for(i in unique(df$annoID)) {
        df$adj.P.value[df$annoID==i] = fdrtool(df[,..plotCol][df$annoID==i], statistic=c("pvalue"), cutoff.method=c("fndr"), plot = F)$qval
      }
    }
        
    zz = ggplot(df, aes(annoName, sumstatName, fill = tempScoreCol)) + geom_tile() + scale_y_discrete(expand = c(0, 0)) + scale_x_discrete(expand = c(0, 0)) + xlab("Cluster") + 
      ylab("Trait") + theme_classic(base_size = plotTextSize) + theme(axis.text.x = element_text(colour = "black", angle = -45, hjust = 0)) + 
      coord_fixed() + theme(legend.title = element_text(size = 10, face = "bold"))    
    print(zz + geom_tile(aes(fill = tempScoreCol)) + scale_fill_gradient2(high="#D7191C", mid="white", low="#2C7BB6", name = "-logP") + geom_text(aes(label = plotLabel), size = plotTextSize * 0.55))
}



compress_scdrs_traits <- function(cds, traits, start, cell_type, window = F, N = 500, cores = F){
    # Func: compress the scDRS for lineages of cds object
    cds_name <- deparse(substitute(cds))
    if(cores != F){
        cl <- makeCluster(cores)
        clusterEvalQ(cl, c(library(evobiR)))
    }

    # step0. Initialize
    scdrs_list <- list()
    for(sub_trait in traits){
        scdrs_list[[sub_trait]] <- data.frame(matrix(nrow = N, ncol = 0)) # column: lineage 1, 2, ...
    }
    
    for(sub_lin in names(cds@lineages)){
        # step1. Obtain object for specific lineage
        input <- paste0("get_lineage_object(",cds_name,", lineage = '", sub_lin, "', start = ", start, ")")
        cds_subset <- eval(parse(text=input))
        cells <- colnames(cds_subset)
    
        
        window <- ncol(cds_subset)/N
        step <- (ncol(cds_subset)-window)/N

        print(sub_lin)
        print(paste0("Window: ", window))
        print(paste0("Step: ", step))
        
   
        # step2. Extract the pseudotime information
        pt <- cds_subset@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
        pt <- pt[order(pt)]

        for(sub_trait in traits){
            # step3. Prepare for scDRS information
            scdrs <- as.data.frame(fread(paste0("files/scdrs/scDRS_sc_", sub_trait, ".csv")))
            rownames(scdrs) <- scdrs[[1]]
            scdrs <- scdrs[-1]
            scdrs <- scdrs[names(pt), "zscore"]
            
            # step4. Compress using sliding window
            scdrs.comp <- SlidingWindow("mean", scdrs, window, step)
            scdrs_list[[sub_trait]] <- cbind(scdrs_list[[sub_trait]], t(data.frame(matrix(scdrs.comp, nrow = 1, ncol = N), row.names = sub_lin)))
        }
    }
    saveRDS(scdrs_list, file = paste0("files/scdrs/scDRS_compressed_zscore_", cell_type,".RDS"))
    return(scdrs_list)
}


obtain_fitted_scdrs <- function(cds, file_scdrs, model = "scdrs ~ splines::ns(pseudotime, df=3)"){
    scdrs_list <- readRDS(file_scdrs)
    
    # step0. Initialize
    scdrs_fitted_list <- scdrs_list

    for(sub_trait in names(scdrs_list)){
        for(sub_lin in colnames(scdrs_list[[sub_trait]])){
            pt <- as.numeric(unlist(cds@pseudotime[[sub_lin]]))
            score <- scdrs_list[[sub_trait]][[sub_lin]]
            df_sel <- data.frame(pseudotime = pt, scdrs = score)
            d <- data.frame(pseudotime = seq(from = 0, to = max(pt), by = max(pt)/(length(pt) - 1)))

            fit <- speedglm(model, data = df_sel, acc = 1e-3, model = F, y = F)
            fitted <- as.numeric(stats::predict(fit, newdata = d, type = "response"))
            scdrs_fitted_list[[sub_trait]][[sub_lin]] <- fitted
        }
    }
    return(scdrs_fitted_list)
}

plot_heatmap_fitted_scdrs <- function(scdrs_list, lineages, branches, n_metacell, out_file, figure_width = 8, figure_height = 5, vmin = -1, vmax = 1, border_color = NA){
    df_scdrs <- data.frame(matrix(ncol=length(lineages)*n_metacell, nrow = 0))
    for(sub_trait in names(trait_info)){
        if(length(lineages) == 2){
            df_sub_trait <- data.frame(matrix(c(rev(scdrs_list[[sub_trait]][[1]]), scdrs_list[[sub_trait]][[2]]), nrow = 1))
        }else{
            df_sub_trait <- data.frame(matrix(as.numeric(unlist(scdrs_list[[sub_trait]])), nrow = 1))
        }
        rownames(df_sub_trait) <- trait_info[sub_trait]
        df_scdrs <- rbind(df_scdrs, df_sub_trait)
    }


    # Annotation for columns
    annotation_col = data.frame(
                                CellType = factor(rep(lineages, each = n_metacell)), 
                                Branch = factor(branches))
    rownames(annotation_col) = colnames(df_scdrs)
    ann_colors <- list(CellType = cols_traj[lineages], Branch = cols_traj[unique(branches)])


    ph <- pheatmap(df_scdrs, 
                   cluster_rows = F, 
                   cluster_cols = F, 
                   scale = "none", 
                   show_rownames = T, 
                   show_colnames = F, 
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   breaks = seq(vmin, vmax, length.out = 100),
                   # color = color_scale_heatmap, 
                   color = colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(100),
                   border_color = border_color)
    save_pheatmap(ph, out_file, width=figure_width, height=figure_height)
}

################## Transcriptional regulator inference ##################
identify_pt_tfs_ratio <- function(net, cds, genes, lineages, minsize = 5, ratio = 0.01){
    pt_genes <- c()
    for(lineage in lineages){
        # Express matrix at metacells
        expr <- t(cds@expression[[lineage]][, genes])
        
        # Infer TF activities using ulm (statistic, source, condition, score, p_value)
        acts <- tryCatch({
              decoupleR::run_ulm(
                mat = expr,
                net = net,
                .source = 'source',
                .target = 'target',
                .mor = 'mor',
                minsize = minsize
              )
            }, error = function(e) {
              message("run_ulm() failed: ", e$message)
              return(NULL)  # or return an empty tibble/data.frame if preferred
            })
        if (!is.null(acts)){
            acts_formatted <- acts %>%
                tidyr::pivot_wider(id_cols = 'source', 
                                   names_from = 'condition',
                                   values_from = 'p_value') %>%
                tibble::column_to_rownames('source')
            acts_formatted <- acts_formatted[, as.character(seq(1, 500, length.out = 500))]
            res <- rownames(acts_formatted[rowSums(acts_formatted < 0.05) >= ratio*500, ])
            
            if(length(res) > 0) pt_genes <- union(pt_genes, res)
        }
    }
    return(pt_genes)
}

# PPI network per module for both TF and gene
ppi_per_celltype_tf_gene <- function(df_celltype, celltype_sel, module_sel, tf_list) {  
  # Add is_TF column based on tf_list
  df_celltype <- df_celltype %>%
    mutate(is_TF = gene %in% tf_list)
  
  # Map genes to STRING identifiers
  mapped_genes <- string_db$map(df_celltype, "gene", removeUnmappedRows = TRUE)
  
  # Retrieve PPI edges
  ppi_df <- string_db$get_interactions(mapped_genes$STRING_id)
  
  # Subset to within-gene list interactions and remove duplicated edges
  ppi_subset <- ppi_df %>%
    filter(from %in% mapped_genes$STRING_id & to %in% mapped_genes$STRING_id) %>%
    rowwise() %>%
    mutate(pair = paste(sort(c(from, to)), collapse = "--")) %>%
    ungroup() %>%
    distinct(pair, .keep_all = TRUE) %>%
    select(-pair)
  
  # Build graph
  ppi_graph <- graph_from_data_frame(ppi_subset, directed = FALSE)
  
  # Map gene names
  gene_map <- mapped_genes %>% select(STRING_id, gene)
  V(ppi_graph)$gene <- gene_map$gene[match(V(ppi_graph)$name, gene_map$STRING_id)]
  
  # Remove unmapped vertices
  ppi_graph <- delete_vertices(ppi_graph, is.na(V(ppi_graph)$gene))
  
  # Add gene type (TF or not)
  V(ppi_graph)$is_TF <- V(ppi_graph)$gene %in% tf_list
  V(ppi_graph)$shape <- ifelse(V(ppi_graph)$is_TF, "square", "circle")
  
  # Degree-based styling
  deg <- degree(ppi_graph)
  V(ppi_graph)$degree <- deg

  # Label size: TFs = 6pt, others = 4pt
  V(ppi_graph)$label.cex <- ifelse(V(ppi_graph)$is_TF, 6/12, 4/12)

  # Label color:
  # - TFs: always black
  # - Genes: black if degree ≥ 5, else gray
  V(ppi_graph)$label.color <- ifelse(
    V(ppi_graph)$is_TF | deg >= 5, "black", "gray80"
  )

  # Font style:
  # - Bold if degree ≥ 5, else normal
  V(ppi_graph)$label.font <- ifelse(deg >= 5, 2, 1)
  V(ppi_graph)$label.family <- "sans"
    
  # ✨ Assign color by component
  comps <- components(ppi_graph)
  V(ppi_graph)$component <- comps$membership
  num_comps <- comps$no
  palette_fn <- colorRampPalette(brewer.pal(8, "Dark2"))
  comp_colors <- setNames(palette_fn(num_comps), as.character(1:num_comps))
  V(ppi_graph)$color <- comp_colors[as.character(V(ppi_graph)$component)]
  
  # Plot
  plot(ppi_graph,
       vertex.label = V(ppi_graph)$gene,
       vertex.label.cex = V(ppi_graph)$label.cex,
       vertex.label.font = V(ppi_graph)$label.font,
       vertex.label.family = V(ppi_graph)$label.family,
       vertex.label.color = V(ppi_graph)$label.color,
       vertex.color = V(ppi_graph)$color,
       vertex.shape = V(ppi_graph)$shape,
       edge.color = "gray40",
       main = paste0("PPI network (", celltype_sel, " : ", module_sel, ")"))

  # Legend
  legend("bottomright", legend = c("TF", "Gene"), pch = c(22, 21),
         pt.bg = "gray80", bty = "n", title = "Node type")

  return(ppi_graph)
}

# PPI network per module for both TF and gene (only visualize those sub-networks size no less than 5)
ppi_per_celltype_tf_gene_sel_5 <- function(df_celltype, celltype_sel, module_sel, tf_list) {
  # df_celltype: gene, extTF
  # Add is_TF column based on tf_list
  df_celltype <- df_celltype %>%
    mutate(is_TF = gene %in% tf_list)
  
  # Map genes to STRING identifiers
  mapped_genes <- string_db$map(df_celltype, "gene", removeUnmappedRows = TRUE)
  
  # Retrieve PPI edges
  ppi_df <- string_db$get_interactions(mapped_genes$STRING_id)
  
  # Subset to within-gene list interactions and remove duplicated edges
  ppi_subset <- ppi_df %>%
    filter(from %in% mapped_genes$STRING_id & to %in% mapped_genes$STRING_id) %>%
    rowwise() %>%
    mutate(pair = paste(sort(c(from, to)), collapse = "--")) %>%
    ungroup() %>%
    distinct(pair, .keep_all = TRUE) %>%
    select(-pair)
  
  # Build graph
  ppi_graph <- graph_from_data_frame(ppi_subset, directed = FALSE)
  
  # Map gene names
  gene_map <- mapped_genes %>% select(STRING_id, gene)
  V(ppi_graph)$gene <- gene_map$gene[match(V(ppi_graph)$name, gene_map$STRING_id)]
  
  # Remove unmapped vertices
  ppi_graph <- delete_vertices(ppi_graph, is.na(V(ppi_graph)$gene))
  
  # Add gene type (TF or not)
  V(ppi_graph)$is_TF <- V(ppi_graph)$gene %in% tf_list
  V(ppi_graph)$shape <- ifelse(V(ppi_graph)$is_TF, "square", "circle")
  
  # Degree-based styling
  deg <- degree(ppi_graph)
  V(ppi_graph)$degree <- deg

  # Label size: TFs = 6pt, others = 4pt
  V(ppi_graph)$label.cex <- ifelse(V(ppi_graph)$is_TF, 6/12, 4/12)

  # Label color:
  # - TFs: always black
  # - Genes: black if degree ≥ 5, else gray
  V(ppi_graph)$label.color <- ifelse(
    V(ppi_graph)$is_TF | deg >= 5, "black", "gray80"
  )

  # Font style:
  # - Bold if degree ≥ 5, else normal
  V(ppi_graph)$label.font <- ifelse(deg >= 5, 2, 1)
  V(ppi_graph)$label.family <- "sans"
    
  # ✨ Identify and filter by component size
  comps <- components(ppi_graph)
  V(ppi_graph)$component <- comps$membership
  
  # Only keep components with ≥ 5 nodes
  comp_sizes <- table(comps$membership)
  keep_comps <- as.integer(names(comp_sizes[comp_sizes >= 5]))
  ppi_graph <- induced_subgraph(ppi_graph, vids = V(ppi_graph)[component %in% keep_comps])
  
  # Recalculate component membership on the filtered graph
  comps <- components(ppi_graph)
  V(ppi_graph)$component <- comps$membership
  num_comps <- comps$no
    
  if(num_comps > 0){
      # Generate color palette for remaining components
      palette_fn <- colorRampPalette(brewer.pal(8, "Dark2"))
      comp_colors <- setNames(palette_fn(num_comps), as.character(1:num_comps))
      V(ppi_graph)$color <- comp_colors[as.character(V(ppi_graph)$component)]
    
      
      # Plot
      plot(ppi_graph,
           vertex.label = V(ppi_graph)$gene,
           vertex.label.cex = V(ppi_graph)$label.cex,
           vertex.label.font = V(ppi_graph)$label.font,
           vertex.label.family = V(ppi_graph)$label.family,
           vertex.label.color = V(ppi_graph)$label.color,
           vertex.color = V(ppi_graph)$color,
           vertex.shape = V(ppi_graph)$shape,
           edge.color = "gray40",
           main = paste0("PPI network (", celltype_sel, " : ", module_sel, ")"))
    
      # Legend
      legend("bottomright", legend = c("TF", "Gene"), pch = c(22, 21),
             pt.bg = "gray80", bty = "n", title = "Node type")
    
      return(ppi_graph)
  }
}

# Classify TF status in network
classify_tf_network_behavior <- function(ppi_graph, tf_list) {
  tf_info <- data.frame(TF = tf_list, status = "Not in network", stringsAsFactors = FALSE)
  
  if (is.null(ppi_graph) || vcount(ppi_graph) == 0) {
    return(tf_info)
  }
  
  # Get components and degree
  comps <- components(ppi_graph)
  deg <- degree(ppi_graph)
  V(ppi_graph)$component <- comps$membership
  V(ppi_graph)$degree <- deg
  
  for (i in seq_along(tf_list)) {
    tf <- tf_list[i]
    v <- V(ppi_graph)[V(ppi_graph)$gene == tf]
    if (length(v) == 1) {
      comp_size <- comps$csize[V(ppi_graph)[v]$component]
      tf_deg <- V(ppi_graph)[v]$degree
      tf_info$status[i] <- dplyr::case_when(
        comp_size >= 5 & tf_deg >= 3 ~ "Present, Comp≥5, Degree≥3",
        comp_size >= 5 & tf_deg < 3  ~ "Present, Comp≥5, Degree<3",
        comp_size < 5                ~ "Present, Comp<5",
        TRUE                         ~ "Present"
      )
    }
  }
  
  return(tf_info)
}


# Plot the presence of key TFs in networks
plotKeyTFsPresence <- function(modules_sel, n_header = 3){
  # step1. Identify key TFs
  keyTFs <- c()
  for(sub_module in modules_sel){
      # PPI networks not less than 5 nodes
      ppi_graph <- readRDS(paste0("data/inhouse/string-db/traDEG_modules/r_PPI_ext_TF_3_5_all_ratio_0.02_", sub_module, ".rds"))
    
      # ✨ Identify and filter by component size
      comps <- components(ppi_graph)
      V(ppi_graph)$component <- comps$membership
      
      # Only keep components with ≥ 5 nodes
      comp_sizes <- table(comps$membership)
      keep_comps <- as.integer(names(comp_sizes[comp_sizes >= 5]))
      ppi_graph <- induced_subgraph(ppi_graph, vids = V(ppi_graph)[component %in% keep_comps])
      
      # Focus on the top 3 TFs with degree ≥ 3
      if(!is.null(ppi_graph)){
          tf_high_deg_sorted <- tibble(gene = V(ppi_graph)$gene, 
                                       degree = V(ppi_graph)$degree, 
                                       is_TF = V(ppi_graph)$is_TF) %>%
          filter(is_TF, degree >= 3) %>%
          arrange(desc(degree)) %>% 
          head(n_header)
          
          keyTFs <- union(keyTFs, tf_high_deg_sorted$gene)
      }
  }

  # step2. Classify the TF presence in PPI networks
  tf_classifications <- lapply(modules_sel, function(mod){
      graphs_by_module <- readRDS(paste0("data/inhouse/string-db/traDEG_modules/r_PPI_ext_TF_3_5_all_ratio_0.02_", mod, ".rds"))
      classify_tf_network_behavior(graphs_by_module, keyTFs) %>%
          mutate(module = mod)
  })

  df_tf_behavior <- bind_rows(tf_classifications)
  df_tf_behavior$TF <- factor(df_tf_behavior$TF, levels = rev(keyTFs))
  df_tf_behavior$module <- factor(df_tf_behavior$module, levels = modules_sel)
  
  df_tf_behavior$status <- factor(
    df_tf_behavior$status,
    levels = c("Present, Comp≥5, Degree≥3",
               "Present, Comp≥5, Degree<3",
               "Present, Comp<5",
               "Not in network")
  )

  # step3. Prepare dot plot of TF presence in PPI networks
  dot_plot <- ggplot(df_tf_behavior, aes(x = module, y = TF, color = status)) +
    geom_point(aes(size = status), shape = 16) +
    scale_color_manual(values = status_colors) +
    scale_size_manual(values = c(5, 3, 2, 1)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "TF Presence and Network Context per Module",
         y = "Transcription Factor",
         x = "Module",
         color = "Network Status",
         size = "Network Status")
  
  # step4. Identify the source of TF
  # Clean and standardize names
  modules_clean <- lapply(modules, function(x) unique(trimws(as.character(x))))
  infTFs_list_clean <- lapply(infTFs_list, function(x) unique(trimws(as.character(x))))
  keyTFs_clean <- unique(trimws(as.character(keyTFs)))

  # Diagnostic function
  safe_membership_trace <- function(tf, mod, lst, source_name) {
    if (!mod %in% names(lst)) {
      message("⚠️ Module '", mod, "' not found in ", source_name)
      return(FALSE)
    }
    contents <- lst[[mod]]
    if (is.null(contents)) {
      message("⚠️ NULL module content in ", source_name, " for ", mod)
      return(FALSE)
    }
    if (!is.character(contents)) {
      message("⚠️ Non-character contents in ", source_name, " for ", mod, ": class=", class(contents))
      return(FALSE)
    }
    return(tf %in% contents)
  }

  # Build presence status safely
  tfPresense <- expand.grid(TF = keyTFs_clean, module = modules_sel, stringsAsFactors = FALSE) %>%
    mutate(
      in_a = mapply(function(tf, mod) safe_membership_trace(tf, mod, modules_clean, "modules"), TF, module),
      in_b = mapply(function(tf, mod) safe_membership_trace(tf, mod, infTFs_list_clean, "infTFs_list"), TF, module),
      status = case_when(
        in_a & in_b ~ "Both",
        in_a & !in_b ~ "Included only",
        !in_a & in_b ~ "Inferred only",
        TRUE ~ "None"
      )
    )
  df_keyTFs <- tfPresense[tfPresense$status != "None", ]
  return(list(dot_plot = dot_plot, df_keyTFs = df_keyTFs))
}

# Longify and min-max normalization
longify_norm <- function(df, value_name, tfs_of_interest, celltypes_ordered) {
  df[, c('transcription_factor', celltypes_ordered)]  %>%
    filter(transcription_factor %in% tfs_of_interest) %>%
    pivot_longer(-transcription_factor, names_to = "celltype", values_to = value_name) %>%
    mutate(celltype = factor(celltype, levels = celltypes_ordered)) %>%
    group_by(transcription_factor) %>%
    mutate(!!value_name := (get(value_name) - min(get(value_name), na.rm = TRUE)) /
                            (max(get(value_name), na.rm = TRUE) - min(get(value_name), na.rm = TRUE))) %>%
    ungroup()
}
                    
# Compute square tile corners
gen_square <- function(x, y) {
  tibble(x = c(x - 0.5, x - 0.5, x + 0.5, x + 0.5),
         y = c(y - 0.5, y + 0.5, y + 0.5, y - 0.5))
}
                    
################## Visium validation ##################
avg_PC1_per_spot <- function(genes){
    df <- data.frame(PC1 = c(), BayesSpace_harmony_09 = c(), Sample_ID = c())
    
    genes <- rowData(spe)$gene_search[
        rowData(spe)$gene_name %in% genes
    ]
    
    if(length(genes) > 0){
        for(sub_sample in unique(spe@int_metadata$imgData$sample_id)){
            # Perform PCA
            p09_sub <- vis_gene_pca(
                    spe = spe,
                    geneid = genes,
                    sampleid = sub_sample,
                    multi_gene_method = "pca",
                    point_size = 1.5
            )
            
            # Extract spatial metadata
            spatial_data <- as.data.frame(colData(spe))
            
            # Ensure that the sample ID matches
            spatial_data <- spatial_data[spatial_data$sample_id == sub_sample, ]
            
            spatial_data$PC1 <- p09_sub
            
            # Assuming "Sp09_cluster" is the column name for the Sp09 cluster values
            ggplot_df <- spatial_data %>%
                dplyr::select(PC1, BayesSpace_harmony_09)
            
            ggplot_df$Sample_ID <- sub_sample
            
            df <- rbind(df, ggplot_df)
        }
    }  
    
    df$BayesSpace_harmony_09 <- factor(df$BayesSpace_harmony_09, levels = c(1, 2, 3, 5, 8, 4, 7, 6, 9))    
    return(df)
}

square_heatmap <- function(A, B = NULL,
                           row_labels = NULL, col_labels = NULL,
                           limits = NULL,
                           size_range = c(0.2, 0.9),
                           fill_sequential = "Reds",
                           fill_diverging = "PRGn",
                           base_fill = "white",
                           base_border = "grey70",
                           square_border = "grey30") {
  
  nrow <- nrow(A)
  ncol <- ncol(A)
  
  if (is.null(B)) B <- matrix(1, nrow = nrow, ncol = ncol)
  if (is.null(row_labels)) row_labels <- paste0("Row", 1:nrow)
  if (is.null(col_labels)) col_labels <- paste0("Col", 1:ncol)
  
  # Build full data
  df <- expand.grid(y = 1:nrow, x = 1:ncol) %>%
    mutate(value = as.vector(unlist(A)),
           size = rescale(as.vector(unlist(B)), to = size_range)) %>%
    rowwise() %>%
    mutate(
      x_center = x - 0.5,
      y_center = y - 0.5,
      # Rescaled size square
      xmin = x_center - size / 2,
      xmax = x_center + size / 2,
      ymin = y_center - size / 2,
      ymax = y_center + size / 2,
      # Base square (always full size)
      base_xmin = x_center - 0.5,
      base_xmax = x_center + 0.5,
      base_ymin = y_center - 0.5,
      base_ymax = y_center + 0.5
    )

  # Determine fill scale
  val_range <- range(df$value, na.rm = TRUE)
  if (is.null(limits)) limits <- val_range
  fill_scale <- if (val_range[1] >= 0) {
    scale_fill_gradientn(
      colors = colorRampPalette(brewer.pal(9, fill_sequential))(100),
      limits = limits,
      oob = squish
    )
  } else {
    scale_fill_gradientn(
      # colors = rev(colorRampPalette(brewer.pal(11, fill_diverging))(100)),
      colors = colorRampPalette(brewer.pal(11, fill_diverging))(100),
      limits = limits,
      oob = squish
    )
  }

  # Plot
  ggplot(df) +
    # Base grid square
    geom_rect(aes(xmin = base_xmin, xmax = base_xmax, ymin = base_ymin, ymax = base_ymax),
              fill = base_fill, color = base_border, linewidth = 0.4) +
    
    # Value square
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = value),
              color = square_border) +
    
    fill_scale +
    coord_fixed() +
    scale_y_reverse(breaks = 1:nrow, labels = row_labels) +
    scale_x_continuous(breaks = 1:ncol, labels = col_labels) +
    theme_minimal(base_size = 14) +
    labs(x = NULL, y = NULL, fill = "Value") +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_text(hjust = 1)
    )
}

square_size_legend <- function(breaks, size_range = c(0.2, 0.9),
                               title = "Size", fill = "grey60", border = "black") {
  scaled_sizes <- rescale(breaks, to = size_range)
  df <- data.frame(
    value = breaks,
    size = scaled_sizes,
    x_center = 0,
    y_center = seq_along(breaks)
  ) %>%
    mutate(
      xmin = x_center - size / 2,
      xmax = x_center + size / 2,
      ymin = y_center - size / 2,
      ymax = y_center + size / 2
    )

  ggplot(df) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = fill, color = border) +
    scale_y_continuous(breaks = df$y_center, labels = df$value) +
    coord_fixed() +
    labs(y = title, x = NULL) +
    theme_void() +
    theme(
      axis.text.y = element_text(size = 10, hjust = 1),
      plot.margin = margin(5, 10, 5, 5)
    )
}