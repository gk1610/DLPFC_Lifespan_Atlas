source("trajectory_reconstruction.r")

modified_moran_test_for_each_lineage <- function(file_cds, lineage, pg_start, factor, N, cores, file_label){
	cds <- readRDS(file_cds)

    # modify functions of graph_test, my.moran.test
    assignInNamespace(x = "graph_test", value = graph_test_lm_aging, ns = "monocle3")
    assignInNamespace(x = "my.moran.test", value = my.moran.test_lm_aging, ns = "monocle3")

	# run the test
	print(lineage)
	cds.sub = get_lineage_object(cds, lineage, pg_start, N = N)
	data = counts(cds.sub)
	rows=rownames(data)[rowSums(data> 0) > factor*ncol(data)]
	cds.sel = cds.sub[rows]
	print(paste0("Testing ", ncol(cds.sel)," cells and ", nrow(cds.sel), " genes"))
	cds_pr_test_res <- monocle3:::graph_test(cds.sel, neighbor_graph="principal_graph", verbose = T, cores = cores)
	save(cds_pr_test_res, file = paste0(file_label, "_Moran_", lineage,".R"))
	write.table(subset(cds_pr_test_res, q_value < 0.05), paste0(file_label, "_traDEG_", lineage,".txt"), sep ="\t", quote = F)
}

# Parameters
args <- commandArgs(T)
file_cds <- args[1]
lineage <- args[2]
pg_start <- as.numeric(args[3])
factor <- as.numeric(args[4]) # 0.05
N <- as.numeric(args[5]) # 10000
cores <- as.numeric(args[6]) # 40
file_label <- args[7]

modified_moran_test_for_each_lineage(file_cds, lineage, pg_start, factor, N, cores, file_label)