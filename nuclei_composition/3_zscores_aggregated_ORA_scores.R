# Pathway Activity Analysis and Heatmap Visualization
# ====================================================
# 
# Overview:
#   This script performs comparative pathway analysis across subclasses using aggregated
#   ORA (Over-Representation Analysis) scores. It identifies brain-related pathways with
#   significant differences between developmental trajectory clusters and generates a
#   publication-quality heatmap.
#
# Input requirements:
#   - ORA_aggregated_Scores_wls.RData: Aggregated ORA scores by subclass
#   - Supplementary_Table3.csv: Metadata with subclass-to-cluster mapping
#
# Outputs:
#   - Statistical test results (t-tests with FDR correction)
#   - Heatmap PDF visualization (normalized pathway activity)

# Load required libraries
library(tidyverse)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(viridis)

# ====== DATA LOADING ======

# Load aggregated ORA scores (from single-nucleus analysis)
load("ORA_aggregated_Scores_wls.RData")

# Convert list to data frame and set subclass rownames
agg_df <- bind_cols(agg_list)
rownames(agg_df) <- agg_df[, 1]
agg_df <- agg_df[, -1]

# Load subclass metadata including developmental trajectory clusters to get which subclass is log increasing or decreasing
supp_table <- read.csv("Supplementary_Table3.csv")

# ====== PATHWAY FILTERING ======
# Filter for brain-related pathways using keyword matching

brain_keywords <- c(
  "brain", "neuron", "neuro", "glia", "synapse", "axon",
  "dendrite", "astrocyte", "myelin", "cns", "nerve",
  "cortex", "prefrontal", "dlpfc", "pfc", "frontal",
  "microglia", "oligodendrocyte", "excitatory", "inhibitory"
)

# Create regex pattern for case-insensitive matching
pattern <- paste(brain_keywords, collapse = "|")

# Filter pathways
brain_related_pathways <- tibble(Pathway = colnames(agg_df)) %>%
  filter(str_detect(tolower(Pathway), pattern))

agg_final_df <- agg_df[, colnames(agg_df) %in% unique(brain_related_pathways$Pathway)]

# ====== DATA PREPARATION FOR STATISTICAL TESTING ======

z_df <- as.data.frame(agg_final_df)

# Ensure rownames are in a column for merging
z_df$Subclass <- rownames(z_df)

# Remove duplicate columns if present
z_df <- z_df[, !duplicated(colnames(z_df))]

# Merge with subclass metadata to get cluster assignments
z_df <- z_df %>%
  left_join(supp_table[, c("Subclass", "cluster")], by = "Subclass")

# Validate merge
if (any(is.na(z_df$cluster))) {
  warning("Some subclasses missing cluster assignments. Check Supplementary_Table3.csv")
}

# Reshape to long format for statistical testing
z_df_long <- reshape2::melt(z_df, id = c("Subclass", "cluster"))
colnames(z_df_long)[3:4] <- c("Pathway", "Score")
z_df_long$Score <- as.numeric(as.character(z_df_long$Score))

# ====== STATISTICAL TESTING: T-TESTS BY PATHWAY ======
# Compare pathway activity between developmental trajectory clusters
# (log_increasing vs. log_decreasing)

t_test_results <- z_df_long %>%
  mutate(cluster = factor(cluster, 
                          levels = c("log_decreasing", "log_increasing"))) %>%
  group_by(Pathway) %>%
  # Filter pathways with representation in both clusters and variance
  filter(n_distinct(cluster) == 2) %>%
  summarize(
    n_log_increasing = sum(cluster == "log_increasing"),
    n_log_decreasing = sum(cluster == "log_decreasing"),
    var_log_increasing = var(Score[cluster == "log_increasing"]),
    var_log_decreasing = var(Score[cluster == "log_decreasing"]),
    t_stat = tryCatch(
      t.test(Score ~ cluster)$statistic,
      error = function(e) NA
    ),
    p_value = tryCatch(
      t.test(Score ~ cluster)$p.value,
      error = function(e) NA
    ),
    mean_diff = tryCatch(
      diff(tapply(Score, cluster, mean)),
      error = function(e) NA
    ),
    .groups = "drop"
  ) %>%
  filter(!is.na(p_value)) %>%  # Remove failed tests
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),  # Benjamini-Hochberg FDR
    direction = case_when(
      mean_diff > 0 ~ "log_increasing",
      mean_diff < 0 ~ "log_decreasing",
      TRUE ~ "no_change"
    )
  ) %>%
  arrange(p_adj)

# ====== SELECT SIGNIFICANT PATHWAYS ======
# Top 15 pathways in each direction (FDR < 0.05)

top_pathways_positive <- t_test_results %>%
  filter(p_adj < 0.05, mean_diff > 0) %>%
  arrange(desc(mean_diff)) %>%
  slice_head(n = 15)

top_pathways_negative <- t_test_results %>%
  filter(p_adj < 0.05, mean_diff < 0) %>%
  arrange(mean_diff) %>%  # More negative = stronger effect
  slice_head(n = 15)

# Combine top pathways
top_pathways_combined <- bind_rows(top_pathways_positive, top_pathways_negative)

# ====== PREPARE HEATMAP DATA ======

# Get subclass order from original aggregated data
subclass_order <- rownames(agg_final_df)

# Create annotation dataframe with cluster assignments
annotation_df <- supp_table %>%
  filter(Subclass %in% subclass_order) %>%
  distinct(Subclass, cluster) %>%
  arrange(match(Subclass, subclass_order))

# Define cluster colors for annotation
cluster_colors <- c(
  "log_increasing" = "#000000",  # black
  "log_decreasing" = "#808080"   # gray
)

rownames(annotation_df) <- annotation_df$Subclass
annotation_df <- as.data.frame(annotation_df)

# Extract matrix for top pathways only
mat <- t(agg_final_df[, colnames(agg_final_df) %in% top_pathways_combined$Pathway])

# Min-max normalize each row (pathway) to [0, 1]
minmax_normalize <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}
norm_mat <- t(apply(mat, 1, minmax_normalize))

# Clean pathway names for display
clean_pathway_names <- function(pathway_names) {
  pathway_names <- gsub("^GOBP_", "", pathway_names)  # Remove GO prefix
  pathway_names <- gsub("_", " ", pathway_names)      # Replace underscores
  pathway_names <- tools::toTitleCase(tolower(pathway_names))  # Title case
  return(pathway_names)
}

rownames(norm_mat) <- clean_pathway_names(rownames(norm_mat))

# Ensure annotation order matches matrix columns
annotation_df <- annotation_df[match(colnames(norm_mat), annotation_df$Subclass), ]

# ====== CREATE HEATMAP ANNOTATION ======

# Add cell class information if available from metadata
# (Modify add_class2_metadata() call if using custom metadata)
annotation_df <- annotation_df %>%
  left_join(
    supp_table %>% select(Subclass, class),
    by = "Subclass"
  )

# Define color map for subclasses (modify as needed)
subclass_colors <- setNames(
  scales::hue_pal()(length(unique(annotation_df$Subclass))),
  unique(annotation_df$Subclass)
)

# Create ComplexHeatmap annotation object
top_annotation <- HeatmapAnnotation(
  df = annotation_df[, c("cluster", "Subclass")],
  col = list(
    cluster = cluster_colors,
    Subclass = subclass_colors
  ),
  annotation_name_side = "left",
  gap = unit(2, "mm")
)

# ====== GENERATE HEATMAP ======

norm_mat <- as.matrix(norm_mat)

# Validate data before plotting
cat("Matrix range:", min(norm_mat, na.rm = TRUE), "to", max(norm_mat, na.rm = TRUE), "\n")
cat("Pathways:", nrow(norm_mat), "| Subclasses:", ncol(norm_mat), "\n")

# Create PDF output
pdf("Extended_figure1f_heatmap.pdf", 
    width = 15, 
    height = 20)

heatmap_plot <- Heatmap(
  norm_mat,
  name = "Normalized\nPathway Activity",
  col = circlize::colorRamp2(c(0, 0.5, 1), viridis::plasma(3)),
  
  # Annotations
  top_annotation = top_annotation,
  
  # Clustering
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  
  # Display options
  show_column_names = FALSE,  # Subclass names in annotation
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 7),
  
  # Legend
  heatmap_legend_param = list(
    title = "Normalized\nActivity Score",
    at = c(0, 0.5, 1),
    legend_height = unit(4, "cm")
  ),
  
  # Layout
  use_raster = TRUE,  # Improves rendering speed for large heatmaps
  raster_quality = 2
)

draw(heatmap_plot)
dev.off()

