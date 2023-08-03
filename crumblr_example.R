library(dreamlet)
library(crumblr)

## Example code for "no duplicates" matrices

pb=readRDS("/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/pseudobulk/AGING_2023-06-09_01_45_PB_SubID_subtype.RDS")
## this is your normalized cell counts matrices
cobj = crumblr(cellCounts(pb))

colData(pb)$scaled_Age=as.numeric(scale(colData(pb)$Age))
colData(pb)$scaled_PMI=as.numeric(scale(colData(pb)$PMI))
metadata(pb)$aggr_means$log_n_counts=as.numeric(log(metadata(pb)$aggr_means$n_counts))

