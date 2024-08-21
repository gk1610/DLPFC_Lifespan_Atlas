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

