.libPaths(c("/sc/arion/projects/psychAD/aging/kiran/RLib_4_3.2/",.libPaths()))
library(zellkonverter)
library(dreamlet)
library(variancePartition)
library(HDF5Array)
library(dplyr)
library(SingleCellExperiment)


args = commandArgs(trailingOnly=TRUE)
brain_bank=args[1]
input_groups=args[2]
model_name="mito_ribo_model"

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

#### dreamlet calculation
colData(pb_subset)$scaled_age=scale(colData(pb_subset)$Age)
pb_subset_final=pb_subset[,colData(pb_subset)$groups==input_groups]

if (input_groups == "Neonatal" | input_groups=="Childhood" | input_groups=="Adolescence") {

form=" ~ scaled_age + Sex + scale(PMI) + log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes"
form=as.formula(form)

} else {

form="~ scaled_age + Source + Sex + scale(PMI) + log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes"
form=as.formula(form)

}

res.proc = processAssays(pb_subset_final,as.formula(form),assays=grep("EN_L5_ET",assayNames(pb_subset_final),value=TRUE,invert=TRUE))
res.dl = dreamlet(res.proc, form)
metadata=as.data.frame(colData(res.proc))

res_mash=run_mash(res.dl, coef="scaled_age")
res.dvar <- diffVar(res.dl)

save(res_mash,res.dvar,metadata,res.dl,input_groups,form,file=paste0("/sc/arion/projects/psychAD/aging/kiran/analysis/dreamlet/four_groups/subclass/four_groups_specific_",input_groups,".RDATA"))


