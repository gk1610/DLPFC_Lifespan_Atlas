# Modify based on kriestein's Science paper (2023)
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

setClass("cell_data_set_ext", contains = "cell_data_set", slots=c(graphs = "list", lineages="list", expression="list", expectation="list", pseudotime="list")) -> cell_data_set_ext

import_monocle <-function(cds){
  cds <- as(cds,"cell_data_set_ext")
  return(cds)
}

path.distance <- function(path){
    # Func: computes the average Euclidean distance between successive points in a given path
    dists=c()
    for(i in 2:nrow(path)){
        x1 = path[i-1,1]
        y1 = path[i-1,2]
        x2 = path[i,1]
        y2 = path[i,2]
        d.x = x2 - x1
        d.y = y2 - y1
        dist = sqrt(d.x*d.x + d.y*d.y)
        dists = append(dist,dists)
    }
    return(mean(dists))
}

cell.selector_sub2 <- function(cell, coords, r){
    # Func: determines whether a given point (or cell) is within a specified radius r from a reference point 
    x2 = cell[1]
    y2 = cell[2]
    d.x = x2 - coords[1]
    d.y = y2 - coords[2]
    dist = sqrt(d.x*d.x + d.y*d.y)
    if(dist <= r){
        return(TRUE)
    }
    else{
        return(FALSE)
    }
}

#' @export
selector_sub <- function(node, cells, r){
    # Func: which cells (points in a 2D space) from a provided set are within a specified radius r from a reference node
    x1 = node[1]
    y1 = node[2]
    res = apply(cells, 1, cell.selector_sub2, coords = c(x1, y1), r = r, simplify = T)
    res = names(res[res == TRUE])
    return(res)
}

cell.selector <- function(path, cells, r, cl){
    # Func: determines which cells (points in a 2D space) from a provided set are within a specified radius r from any point in a given path.
    sel.cells = c()
    sel.cells = pbapply(path, 1, selector_sub, cells = cells, r = r, cl = cl, simplify = T)
    return(unique(unlist(sel.cells)))
}

make_graph <- function(sub.graph){
    # Func: constructs a graph object using the igraph package based on a sequence of nodes provided in sub.graph. 
    edges = names(sub.graph)
    start.edges = c()
    end.edges = c()
    for(i in 1:(length(edges)-1)){
        start.edges = append(start.edges, edges[i])
        end.edges = append(end.edges, edges[i+1])
    }
    d = cbind(start.edges, end.edges)
    g = graph_from_data_frame(d, directed = F)
    return(g)
}

included <- function(graph, include_nodes){
    # Func: checks if all nodes specified in include_nodes are present in the given graph.
    all(include_nodes %in% names(graph))
}

connect_nodes <- function(cds, node1, node2){
    # Func: add an edge between two specified nodes (node1 and node2) in the principal graph of a cds object (presumably from monocle or a related single-cell analysis package).
    graph.old = cds@principal_graph[["UMAP"]]
    graph.new <- add_edges(graph.old, c(node1, node2))
    cds@principal_graph[["UMAP"]] <- graph.new
    return(cds)
}

isolate_graph_sub <- function(cds, start, end, lineage, include_nodes = NULL){
    # Func: extract a subgraph (or path) from a principal graph of a cds object (likely from a single-cell analysis package like Monocle). 
    # This subgraph would be the shortest path between a start and end node, potentially filtered to include specific nodes.
    #get lineage graph
    reduction_method = "UMAP"
    graph = cds@principal_graph[[reduction_method]]
    #select cells that are 1) progenitor cells from the region of interest (MGE, CGE) or 2) lineage-committed cells
    sub.graph = all_simple_paths(graph, paste0("Y_", start), paste0("Y_", end))
    if(length(include_nodes) > 0){
          sub.graph = sub.graph[sapply(sub.graph, included, include_nodes = include_nodes)]
    }
    lengths = lengths(sub.graph)
    #get the shortest path
    n = which(lengths==min(lengths))[1]
    sub.graph = sub.graph[[n]]
}

isolate_graph <- function(cds, start, end, lineage, include_nodes = NULL){
    # Func: take a cds object, a starting node (start), an ending node (end), a lineage, and potentially some nodes that need to be included in the path (include_nodes). 
    # The function then modifies the cds object by isolating a graph between the start and end nodes and saves this subgraph back into the cds object.
    #get lineage graph
    cds_name = deparse(substitute(cds))
    sub.graph = isolate_graph_sub(cds, start, end, lineage, include_nodes = include_nodes)
    input = paste0(cds_name, "@graphs$", lineage, " <- make_graph(sub.graph)")
    eval(parse(text=input))
    eval(parse(text=paste0("return(", cds_name, ")")))
}

isolate_lineage_sub <- function(cds, lineage, sel_clusters = F, start_regions = F, starting_clusters = NULL, subset = FALSE, N = 5, cl = 1){
    # Func:  isolate specific cells from a lineage based on certain criteria and then return the selected cell names.
    cds_name = deparse(substitute(cds))
    input = paste0("sub.graph = ",cds_name,"@graphs$", lineage)
    eval(parse(text=input))
    nodes_UMAP = cds@principal_graph_aux[["UMAP"]]$dp_mst
    if(subset == F){
        nodes_UMAP.sub = as.data.frame(t(nodes_UMAP[,names(V(sub.graph))]))
    }
    else{
        g = principal_graph(cds)[["UMAP"]]
        dd = degree(g)
        names1 = names(dd[dd > 2 | dd == 1])
        names2 = names(dd[dd == 2])
        names2 = sample(names2, length(names2)/subset, replace = F)
        names = c(names1, names2)
        names = intersect(names(V(sub.graph)), names)
        nodes_UMAP.sub = as.data.frame(t(nodes_UMAP[,names]))
    }
    #select cells along the graph
    mean.dist = path.distance(nodes_UMAP.sub)
    r = mean.dist*N
    cells_UMAP = as.data.frame(reducedDims(cds)["UMAP"])
    cells_UMAP = cells_UMAP[,c("UMAP_1", "UMAP_2")]
    sel.cells = cell.selector(nodes_UMAP.sub, cells_UMAP, r, cl = cl)
    #only keep cells in the progenitor and lineage-specific clusters
    sel.cells1 = c()
    sel.cells2 = sel.cells
    if(starting_clusters != F){
        sel.cells1 = names(cds@"clusters"[["UMAP"]]$clusters[cds@"clusters"[["UMAP"]]$clusters %in% starting_clusters])
    }
    if(start_regions != F){
        sel.cells1 = sel.cells1[sel.cells1 %in% rownames(cds@colData[cds@colData$region %in% start_regions,])]
    }
    if(length(sel_clusters) > 0){
        sel.cells2 = names(cds@"clusters"[["UMAP"]]$clusters[cds@"clusters"[["UMAP"]]$clusters %in% sel_clusters])
    }
    cells = unique(c(sel.cells1, sel.cells2))
    sel.cells = sel.cells[sel.cells %in% cells]
    return(sel.cells)
}

isolate_lineage <- function(cds, lineage, sel_clusters = NULL, start_regions = F, starting_clusters = F, subset = FALSE, N = 5, cl = 1){
    # Func: isolating cells from a given lineage based on various criteria and then storing the selected cell names into the provided cds object under a specific slot
    cds_name = deparse(substitute(cds))
    sel.cells = isolate_lineage_sub(cds, lineage, sel_clusters = sel_clusters, start_regions = start_regions, starting_clusters = starting_clusters, subset = subset, N = N, cl = cl)
    input = paste0(cds_name, "@lineages$", lineage, " <- sel.cells")
    eval(parse(text=input))
    return(cds)
}

get_lineage_object <- function(cds, lineage = FALSE, start, N = FALSE){
    # Func: create a subset of the main cds object based on specific lineage criteria and other inputs. 
    cds_name = deparse(substitute(cds))
    if(lineage != FALSE){
        input = paste0("sub.graph = ",cds_name,"@graphs$", lineage)
        eval(parse(text=input))
        input = paste0("sel.cells = ",cds_name,"@lineages$", lineage)
        eval(parse(text=input))
    }
    else{
        sel.cells = colnames(cds)
    }
    sel.cells = sel.cells[sel.cells %in% colnames(cds)]
    nodes_UMAP = cds@principal_graph_aux[["UMAP"]]$dp_mst
    if(N != FALSE){
        if(N < length(sel.cells)){
            sel.cells = sample(sel.cells, N)
        }
    }
    #subset the moncole object
    cds_subset = cds[,sel.cells]
    #set the graph, node and cell UMAP coordinates
    if(lineage == FALSE){
        sub.graph = principal_graph(cds_subset)[["UMAP"]]
    }
    principal_graph(cds_subset)[["UMAP"]] <- sub.graph
    cds_subset@principal_graph_aux[["UMAP"]]$dp_mst <- nodes_UMAP[,names(V(sub.graph))]
    cds_subset@clusters[["UMAP"]]$partitions <- cds_subset@clusters[["UMAP"]]$partitions[colnames(cds_subset)]
    #recalculate closest vertex for the selected cells
    cells_UMAP = as.data.frame(reducedDims(cds_subset)["UMAP"])
    closest_vertex = apply(cells_UMAP[,c("UMAP_1", "UMAP_2")], 1, calculate_closest_vertex, nodes = as.matrix(nodes_UMAP[,names(V(sub.graph))]))
    closest_vertex = as.data.frame(closest_vertex)
    cds_subset@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex <- closest_vertex
    source("learn_graph.R")
    cds_subset <- project2MST(cds_subset, project_point_to_line_segment, F, T, "UMAP", nodes_UMAP[,names(V(sub.graph))])
    cds_subset <- order_cells(cds_subset, root_pr_nodes = c(paste0("Y_", as.character(start))))
    return(cds_subset)
}

calculate_closest_vertex <- function(cells, nodes){
    # Func: calculates the closest vertex (node) to a given point based on their coordinates.
    new.pos = as.numeric(cells)
    nearest.idx <- which.min(colSums((nodes - new.pos)^2))
    out = as.integer(gsub("Y_", "", names(nearest.idx)))
}

combine_lineages <- function(cds, start){
    # Func: integrate multiple lineage graphs from the main cds object into a single combined graph
    cds_name = deparse(substitute(cds))
    lineage = names(cds@lineages)[1]
    input = paste0(cds_name, "@graphs$", lineage)
    if(length(names(cds@lineages)) > 1){
        for(lineage in names(cds@lineages)[2:length(names(cds@lineages))]){
            input = paste0(input, ",", cds_name,"@graphs$", lineage)
        }
        input = paste0("igraph::union(", input, ")")
    }
    g = eval(parse(text=input))
    nodes_UMAP = cds@principal_graph_aux[["UMAP"]]$dp_mst
    principal_graph(cds)[["UMAP"]] <- g
    cds@principal_graph_aux[["UMAP"]]$dp_mst <- nodes_UMAP[,names(V(g))]
    cells_UMAP = as.data.frame(reducedDims(cds)["UMAP"])
    closest_vertex = apply(cells_UMAP[,c("UMAP_1", "UMAP_2")], 1, calculate_closest_vertex, nodes = as.matrix(nodes_UMAP[,names(V(g))]))
    closest_vertex = as.data.frame(closest_vertex)
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex <- closest_vertex
    source("learn_graph.R")
    cds <- project2MST(cds, project_point_to_line_segment, F, T, "UMAP", nodes_UMAP[,names(V(g))])
    cds <- order_cells(cds, root_pr_nodes = as.character(paste0("Y_",start)))
    return(cds)
}

compress_lineages <- function(cds, start, window = F, N = 500, cores = F){
    # Func:  iterates through each lineage in your cds object and performs some compression or processing operation on it.
    lineages = names(cds@lineages)
    for(lineage in lineages){
        print(lineage)
        cds = compress_lineage(cds, lineage = lineage, start = start, window = window, gene = FALSE, N = N, cores = cores)
        gc()
    }
    return(cds)
}

#' @export
compress_lineage <- function(cds, lineage, start, window = F, gene = FALSE, N = 500, cores = F, cells = FALSE){
    # Func: a helper function that focuses on compressing the expression of a particular lineage in your main data object (cds). 
    cds_name = deparse(substitute(cds))
    if(gene == FALSE){
        input = paste0("compress_expression(",cds_name,", lineage = '", lineage, "', start = ", start, ", window = ", window, ", gene = ", gene, ", N = ", N, ", cores = ", cores, ")")
    }
    else{
        input = paste0("compress_expression(",cds_name,", lineage = '", lineage, "', start = ", start, ", window = ", window, ", gene = '", gene, "', N = ", N, ", cores = ", cores, ")")
    }
    exp = eval(parse(text=input))
    input = paste0(cds_name, "@expression$", lineage, " <- exp$expression")
    eval(parse(text=input))
    input = paste0(cds_name, "@expectation$", lineage, " <- exp$expectation")
    eval(parse(text=input))
    input = paste0(cds_name, "@pseudotime$", lineage, " <- exp$pseudotime")
    eval(parse(text=input))
    eval(parse(text=paste0("return(",cds_name, ")")))
}

#' @export
compress_expression <- function(cds, lineage, start, window = F, gene = FALSE, N = 500, cores = F){
    # Func: a helper function that aims to process, compress, and fit models to expression data for a particular lineage of your main data object (cds).
    cds_name = deparse(substitute(cds))
    if(cores != F){
        cl <- makeCluster(cores)
        clusterEvalQ(cl, c(library(evobiR)))
    }
    input = paste0("get_lineage_object(",cds_name,", lineage = '", lineage, "', start = ", start, ")")
    cds_subset = eval(parse(text=input))
    family = stats::quasipoisson()
    model = "expression ~ splines::ns(pseudotime, df=3)"
    names(cds_subset) <- rowData(cds_subset)$gene_short_name
    exp = as_matrix(exprs(cds_subset))
    exp = t(exp) /  pData(cds_subset)[, 'Size_Factor']
    pt <- cds_subset@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
    pt = pt[order(pt)]
    exp = exp[names(pt),]
    if(window == FALSE){
        window = nrow(exp)/N
    }
    step = ((nrow(exp)-window)/N)
    #use sliding window to compress expression values and pseudotime
    print(paste0("Window: ", window))
    print(paste0("Step: ", step))
    pt.comp = SlidingWindow("mean", pt, window, step)
    max.pt = max(pt.comp)
    if(gene != F){
        exp.comp = compress(exp[,gene], window = window, step = step)
    }
    else{
        print(paste0("Compressing lineage ", lineage, " and fitting curves"))
        exp.comp = pbapply(exp, 2, compress, window = window, step = step)
    }
    if(gene != F){
        exp_data.sel = cbind(pt.comp, exp.comp)
        exp_data.sel = as.data.frame(exp_data.sel)
        colnames(exp_data.sel) <- c("pseudotime", "expression")
        exp_data.sel$pseudotime <- as.numeric(as.character(exp_data.sel$pseudotime))
        exp_data.sel$expression <- as.numeric(as.character(exp_data.sel$expression))
        d = as.data.frame(seq(from=0, to=max.pt, by = max.pt/(N-1)))
        colnames(d) <- c("pseudotime")
        tryCatch({fit = speedglm(model, data = exp_data.sel, family = family, acc=1e-3, model=FALSE, y=FALSE)
        fit = predict(fit, newdata = d, type='response')
        }, error=function(cond) {fit = as.data.frame(rep(0, N))})
        exp = as.data.frame(cbind(exp.comp, fit, d))
        colnames(exp) <- c("expression", "expectation", "pseudotime")
        exp$expression[exp$expression < 0] <- 0
        exp$expectation[exp$expectation < 0] <- 0
    }
    else{
        d = as.data.frame(seq(from=0, to=max.pt, by = max.pt/(N-1)))
        fit = pbapply(exp.comp, 2, fit.m3, pt = d, max.pt = max(d), N = N)
        fit = apply(fit, 2, as.numeric)
        return(list("expression" = exp.comp, "expectation" = fit, "pseudotime" = d))
    }
    exp$expression[exp$expression < 0] <- 0
    exp$expectation[exp$expectation < 0] <- 0
    if(cores != F){
        stopCluster(cl)
    }
    return(exp)
}

fit.m3 <- function(exp.sel, pt, max.pt, model = "expression ~ splines::ns(pseudotime, df=3)", N = 500){
    require(speedglm)
    family = stats::quasipoisson()
    exp_data.sel = cbind(pt, exp.sel)
    colnames(exp_data.sel) <- c("pseudotime","expression")
    exp_data.sel = as.data.frame(exp_data.sel)
    exp_data.sel$pseudotime <- as.numeric(as.character(exp_data.sel$pseudotime))
    exp_data.sel$expression <- as.numeric(as.character(exp_data.sel$expression))
    tryCatch({fit = speedglm(model, data = exp_data.sel, family = family, acc=1e-3, model=FALSE, y=FALSE)
    d = as.data.frame(seq(from=0, to=max.pt, by = max.pt/(N-1)))
    colnames(d) <- c("pseudotime")
    fit = stats::predict(fit, newdata=d, type="response")
    return(fit)
    }, error=function(cond) {return(rep("NA", N))})
}

as_matrix <- function(mat){
    tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
    
    row_pos <- mat@i+1
    col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
    val <- mat@x
    
    for (i in seq_along(val)){
        tmp[row_pos[i],col_pos[i]] <- val[i]
    }

    row.names(tmp) <- mat@Dimnames[[1]]
    colnames(tmp) <- mat@Dimnames[[2]]
    return(tmp)
}

compress <- function(df, window, step){
    df.comp = SlidingWindow("mean", df, window, step)
}


graph_test_lm_aging <- function (cds, neighbor_graph = c("knn", "principal_graph"),
    reduction_method = "UMAP", k = 25, method = c("Moran_I"),
    alternative = "greater", expression_family = "quasipoisson",
    cores = 1, verbose = FALSE) 
{
    neighbor_graph <- match.arg(neighbor_graph)
    lw <- monocle3:::calculateLW(cds, k = k, verbose = verbose, neighbor_graph = neighbor_graph, 
        reduction_method = reduction_method)
    if (verbose) {
        message("Performing Moran's I test: ...")
    }
    exprs_mat <- SingleCellExperiment::counts(cds)[, attr(lw, 
        "region.id"), drop = FALSE]
    sz <- size_factors(cds)[attr(lw, "region.id")]
    wc <- spdep::spweights.constants(lw, zero.policy = TRUE, 
        adjust.n = TRUE)
    test_res <- pbmcapply::pbmclapply(row.names(exprs_mat), FUN = function(x, 
        sz, alternative, method, expression_family) {
        exprs_val <- exprs_mat[x, ]
        if (expression_family %in% c("uninormal", "binomialff")) {
            exprs_val <- exprs_val
        }
        else {
            exprs_val <- log10(exprs_val/sz + 0.1)
        }
        df = cbind(as.data.frame(exprs_val), colData(cds)$Batch, colData(cds)$Sex, colData(cds)$PMI, log(colData(cds)$nFeature_RNA))
        colnames(df) <- c("exp", "Batch", "Sex", "PMI", "log_n_genes")
        test_res <- tryCatch({
            if (method == "Moran_I") {
                mt <- suppressWarnings(monocle3:::my.moran.test(df, lw, wc, alternative = alternative))
                data.frame(status = "OK", p_value = mt$p.value, 
                  morans_test_statistic = mt$statistic, morans_I = mt$estimate[["Moran I statistic"]])
            }
            else if (method == "Geary_C") {
                gt <- suppressWarnings(my.geary.test(exprs_val, 
                  lw, wc, alternative = alternative))
                data.frame(status = "OK", p_value = gt$p.value, 
                  geary_test_statistic = gt$statistic, geary_C = gt$estimate[["Geary C statistic"]])
            }
        }, error = function(e) {
            data.frame(status = "FAIL", p_value = NA, morans_test_statistic = NA, 
                morans_I = NA)
        })
    }, sz = sz, alternative = alternative, method = method, expression_family = expression_family, 
        mc.cores = cores, ignore.interactive = TRUE)
    if (verbose) {
        message("returning results: ...")
    }
    test_res <- do.call(rbind.data.frame, test_res)
    row.names(test_res) <- row.names(cds)
    test_res <- merge(test_res, rowData(cds), by = "row.names")
    row.names(test_res) <- test_res[, 1]
    test_res[, 1] <- NULL
    test_res$q_value <- 1
    test_res$q_value[which(test_res$status == "OK")] <- stats::p.adjust(subset(test_res, 
        status == "OK")[, "p_value"], method = "BH")
    test_res$status = as.character(test_res$status)
    test_res[row.names(cds), ]
}

my.moran.test_lm_aging <- function (x, listw, wc, alternative = "greater", randomisation = TRUE)
{
    zero.policy = TRUE
    adjust.n = TRUE
    na.action = stats::na.fail
    drop.EI2 = FALSE
    xname <- deparse(substitute(x))
    wname <- deparse(substitute(listw))
    NAOK <- deparse(substitute(na.action)) == "na.pass"
    x <- na.action(x)
    na.act <- attr(x, "na.action")
    if (!is.null(na.act)) {
        subset <- !(1:length(listw$neighbours) %in% na.act)
        listw <- subset(listw, subset, zero.policy = zero.policy)
    }
    n <- length(listw$neighbours)
    S02 <- wc$S0 * wc$S0
    model <- lm(exp ~ Batch + Sex + as.numeric(PMI) + as.numeric(log_n_genes), data = x)
    res <- spdep::lm.morantest(model, listw, zero.policy = zero.policy, alternative = alternative)
    statistic = as.numeric(res[1])
    names(statistic) <- "Moran I statistic standard deviate"
    PrI = as.numeric(res[2])
    vec <- c(res[3]$estimate[1], res[3]$estimate[2], res[3]$estimate[3])
    names(vec) <- c("Moran I statistic", "Expectation", "Variance")
    method <- paste("Moran I test under", ifelse(randomisation, 
        "randomisation", "normality"))
    res <- list(statistic = statistic, p.value = PrI, estimate = vec)
    if (!is.null(na.act)) 
        attr(res, "na.action") <- na.act
    class(res) <- "htest"
    res
}