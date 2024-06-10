.libPaths(c("/sc/arion/projects/psychAD/aging/kiran/RLib_4_3.2",.libPaths()))
library(variancePartition)
library(dreamlet)

args = commandArgs(trailingOnly=TRUE)
brain_bank="all"
input_groups=args[1]

################################################################################################################################
pb=readRDS("/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/pseudobulk/AGING_2024-02-01_22_23_PB_SubID_subclass.RDS")
colData(pb)$SubID=rownames(colData(pb))

#### keep final samples
samples_to_keep=read.table("/sc/arion/projects/CommonMind/aging/resources/AGING_2024-02-01_22_23_processAssays_SubID_subclass.txt")
pb_subset=pb[,colData(pb)$SubID %in% samples_to_keep$V1]

colData(pb_subset)$groups="NA"
colData(pb_subset)$groups[colData(pb_subset)$Age<1]="Childhood"
colData(pb_subset)$groups[colData(pb_subset)$Age>=1 & colData(pb_subset)$Age<12]="Childhood"
colData(pb_subset)$groups[colData(pb_subset)$Age>=12 & colData(pb_subset)$Age<20]="Childhood"
colData(pb_subset)$groups[colData(pb_subset)$Age>=20 & colData(pb_subset)$Age<40]="Young_Adulthood"
colData(pb_subset)$groups[colData(pb_subset)$Age>=40 & colData(pb_subset)$Age<60]="Middle_Adulthood"
colData(pb_subset)$groups[colData(pb_subset)$Age>=60]="Late_Adulthood"
colData(pb_subset)$scaled_age=scale(colData(pb_subset)$Age)

for (input_groups in c("Childhood","Young_Adulthood","Middle_Adulthood","Late_Adulthood")) {

pb_subset_final=pb_subset[,colData(pb_subset)$groups==input_groups]

if (input_groups == "Neonatal" | input_groups=="Childhood" | input_groups=="Adolescence") {

form ="~ scaled_age + Sex + scale(PMI) + log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes"
sex_specific_form = " ~ as.numeric(Sex)*scaled_age + scale(PMI) + log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes"

} else {

form ="~ scaled_age + Source + Sex + scale(PMI) + log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes"
sex_specific_form = " ~ as.numeric(Sex)*scaled_age + Source + scale(PMI) + log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes"
}

# age specific analysis
form=as.formula(form)
res.proc = processAssays(pb_subset_final,as.formula(form),assays=grep("EN_L5_ET",assayNames(pb_subset_final),value=TRUE,invert=TRUE))
vp.list =  fitVarPart(res.proc,as.formula(form))

# sex*age specific analysis

sex_specific_form=as.formula(sex_specific_form)
res.proc = processAssays(pb_subset_final,as.formula(sex_specific_form),assays=grep("EN_L5_ET",assayNames(pb_subset_final),value=TRUE,invert=TRUE))
sex_specific_vp.list =  fitVarPart(res.proc,as.formula(sex_specific_form))

save(sex_specific_vp.list,vp.list,sex_specific_form,form,file=paste0("/sc/arion/projects/CommonMind/aging/analysis/varpart/four_groups/varpart_age_sex_specific_",input_groups,".RDATA"))


}
