# Aggregate ORA Pathway Scores by Cell Subclass
# =============================================
#
# Overview:
#   This script aggregates DLPFC lifespan ORA brain-specific pathway activity scores by subclass.
#   It merges subclass metadata and computes mean pathway activity per subclass for downstream analysis.
#
# Input files:
#   - Aging_nuclei_metadata.csv: Single-nucleus metadata (cellID, subclass)
#   - ORA_brain_related_GO_pathways.csv: Brain-related pathway activity scores (nuclei × pathways)
#
# Output:
#   - aggregated_brain_pathway_scores.RData: Mean pathway activity by subclass (26 subclasses × 374 pathways)

library(tidyverse)
library(data.table)

# ====== LOAD LIFESPAN METADATA ======

# Load single-nucleus metadata
# Expected columns: cellID (or V1 if unnamed), Age, Sex, subclass
nuclei_metadata <- fread("Aging_nuclei_metadata.csv")  # syn62147251 h5ad
nuclei_metadata <- as.data.frame(nuclei_metadata)

# Handle unnamed first column (common with fread)
if ("V1" %in% colnames(nuclei_metadata)) {
  colnames(nuclei_metadata)[colnames(nuclei_metadata) == "V1"] <- "cellID"
}

# Validate required columns exist
required_cols <- c("cellID", "subclass")
missing_cols <- setdiff(required_cols, colnames(nuclei_metadata))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

cat("Nuclei metadata loaded:\n")
cat("  Dimensions:", nrow(nuclei_metadata), "nuclei ×", ncol(nuclei_metadata), "columns\n")
cat("  Columns:", paste(colnames(nuclei_metadata), collapse = ", "), "\n")

# ====== FILTER LOW-ABUNDANCE SUBCLASSES ======

# Remove subclasses with insufficient nuclei
nuclei_metadata <- nuclei_metadata %>%
  filter(!subclass %in% c("EN_L5_ET"))  # n < 10, insufficient for statistical testing

cat("\nAfter filtering:\n")
cat("  Total nuclei:", nrow(nuclei_metadata), "\n")
cat("  Unique subclasses:", n_distinct(nuclei_metadata$subclass), "\n")

subclass_counts <- nuclei_metadata %>% count(subclass) %>% arrange(n)
cat("\n  Subclass sizes:\n")
print(subclass_counts)

# ====== LOAD ORA PATHWAY SCORES ======

# Load ORA activity scores for brain-related pathways
# Rows: individual nuclei (cellID)
# Columns: pathway activity scores (should start with GOBP_)
pathway_scores <- fread("ORA_brain_related_GO_pathways.csv")
pathway_scores <- as.data.frame(pathway_scores)

# Handle unnamed first column
if ("V1" %in% colnames(pathway_scores)) {
  rownames(pathway_scores) <- pathway_scores$V1
  pathway_scores <- pathway_scores[, -1]
}

cat("\nPathway scores loaded:\n")
cat("  Dimensions:", nrow(pathway_scores), "nuclei ×", ncol(pathway_scores), "pathways\n")

# Validate pathway naming
pathway_prefixes <- substr(colnames(pathway_scores), 1, 5)
if (!any(pathway_prefixes == "GOBP_")) {
  warning("No pathways with GOBP_ prefix found. Check pathway naming convention.")
}

# ====== ALIGN NUCLEI ACROSS DATASETS ======

# Keep only nuclei present in BOTH pathway scores and metadata
keep_nuclei <- intersect(rownames(pathway_scores), nuclei_metadata$cellID)

cat("\nAlignment check:\n")
cat("  Pathway score nuclei:", nrow(pathway_scores), "\n")
cat("  Metadata nuclei:", nrow(nuclei_metadata), "\n")
cat("  Nuclei in both datasets:", length(keep_nuclei), "\n")

if (length(keep_nuclei) == 0) {
  stop("No matching nuclei found between pathway scores and metadata. Check cellID column names.")
}

# Subset both datasets to matching nuclei
nuclei_metadata_aligned <- nuclei_metadata %>%
  filter(cellID %in% keep_nuclei) %>%
  arrange(cellID)

pathway_scores_aligned <- pathway_scores[keep_nuclei, ]

# Validate alignment
stopifnot(
  identical(rownames(pathway_scores_aligned), nuclei_metadata_aligned$cellID),
  nrow(pathway_scores_aligned) == nrow(nuclei_metadata_aligned)
)

cat("  ✓ Data alignment validated\n")

# ====== MERGE PATHWAY SCORES WITH SUBCLASS LABELS ======

# Create data frame with nuclei IDs, subclass labels, and pathway scores
merged_df <- data.frame(
  cellID = nuclei_metadata_aligned$cellID,
  subclass = nuclei_metadata_aligned$subclass,
  pathway_scores_aligned
)

cat("\nMerged data:\n")
cat("  Dimensions:", nrow(merged_df), "nuclei ×", 
    ncol(merged_df) - 2, "pathways + 2 metadata columns\n")

# ====== AGGREGATE BY SUBCLASS ======

# Compute mean pathway activity for each subclass
# This produces one row per subclass with mean activity across all nuclei in that subclass
agg_df <- merged_df %>%
  group_by(subclass) %>%
  summarize(across(starts_with("GOBP_"), mean, na.rm = TRUE), 
            .groups = "drop") %>%
  as.data.frame()

# Set subclass as rownames for downstream analysis
rownames(agg_df) <- agg_df$subclass
agg_df <- agg_df[, -1]  # Remove subclass column (now in rownames)

cat("\nAggregated pathway scores:\n")
cat("  Dimensions:", nrow(agg_df), "subclasses ×", ncol(agg_df), "pathways\n")
cat("  Subclasses included:\n")
print(rownames(agg_df))


# ====== SAVE RESULTS ======

# Save aggregated scores for downstream analysis
save(
  agg_df,
  file = "aggregated_brain_pathway_scores.RData"
)

# Also save in CSV format for ease of access
write.csv(agg_df, "aggregated_brain_pathway_scores.csv")
