.libPaths(c("/sc/arion/projects/CommonMind/kiran/RLib_4_3.2",.libPaths()))
library(zellkonverter)
library(dreamlet)
library(variancePartition)
library(HDF5Array)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
brain_bank=args[1]
model_name=args[2]

################################################################################################################################
if (brain_bank=="all") {

model_list=list("linear"="~ Age + Sex + Source + scale(PMI) + log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes",
"log_Age"="~ log2(Age+1) + Sex + Source + scale(PMI) + log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes",

"poly_with_two_df"= "~ poly(Age,df=2) + Source + Sex + scale(PMI)+ log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes",
"poly_with_two_df_log_Age"= "~ poly(log2(Age+1),df=2) + Source + Sex + scale(PMI) + log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes",

"poly_with_three_df"= "~ poly(Age,df=3) + Source + Sex + scale(PMI)+ log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes",
"poly_with_three_df_log_Age"= "~ poly(log2(Age+1),df=3) + Source + Sex + scale(PMI) + log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes",

"poly_with_four_df"= "~ poly(Age,df=4) + Source + Sex + scale(PMI)+ log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes",
"poly_with_four_df_log_Age"= "~ poly(log2(Age+1),df=4) + Source + Sex + scale(PMI) + log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes",

"poly_with_five_df"= "~ poly(Age,df=5) + Source + Sex + scale(PMI)+ log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes",
"poly_with_five_df_log_Age"= "~ poly(log2(Age+1),df=5) + Source + Sex + scale(PMI) + log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes",

"poly_with_six_df"= "~ poly(Age,df=6) + Source + Sex + scale(PMI)+ log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes",
"poly_with_six_df_log_Age"= "~ poly(log2(Age+1),df=6) + Source + Sex + scale(PMI) + log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes")

} else {

model_list=list("linear"="~ Age + (1 | prep) + (1|Sex) + scale_PMI",
"linear_piecewise"=" ~ Age +  I((Age-2)*(Age >= 2)) +  I((Age-12)*(Age >= 12)) + I((Age-20)*(Age >= 20)) + I((Age-40)*(Age >= 40)) + I((Age-60)*(Age >= 60)) + (1 | prep) + (1|Sex) + scale_PMI",
"linear_piecewise_sex_specific"= " ~ Age +  I((Age-2)*(Age >= 2)) +  I((Age-12)*(Age >= 12)) + I((Age-20)*(Age >= 20)) + I((Age-40)*(Age >= 40)) + I((Age-60)*(Age >= 60)) + (1 | prep) + scale_PMI")

}

### collecting covariates for model
pb=readRDS("/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/pseudobulk/AGING_2024-02-01_22_23_PB_SubID_subclass.RDS")
form = as.formula("~ Source + Sex + scale(PMI) + log(n_genes) + percent_mito + mito_genes + mito_ribo + ribo_genes")
colData(pb)$SubID=rownames(colData(pb))

#### keep final samples
samples_to_keep=read.table("/sc/arion/projects/CommonMind/aging/resources/AGING_2024-02-01_22_23_processAssays_SubID_subclass.txt")
pb_subset=pb[,colData(pb)$SubID %in% samples_to_keep$V1]

#### dreamlet calculation

form=as.formula(model_list[[model_name]])
res.proc = processAssays(pb_subset,form, assays=grep("EN_L5_ET",assayNames(pb_subset),value=TRUE,invert=TRUE))
print(dim(res.proc[[1]]$E))
res.dl = dreamlet(res.proc, form)
metadata=as.data.frame(colData(res.proc))

n_samples=dim(res.proc[[1]]$E)[2]
bic <- -2 * res.dl[[1]]$logLik + res.dl[[1]]$edf * log(n_samples)
aic = -2 * res.dl[[1]]$logLik  + 2 * res.dl[[1]]$edf


coef_list=grep("poly",coefNames(res.dl),value=TRUE)

res_mash=list()
for ( i in (1:length(coef_list))) {
res_mash[[i]]=run_mash(res.dl, coef=coef_list[i])
}


save(res_mash,aic,bic,n_samples,metadata,res.dl,form,file=paste0("/sc/arion/projects/psychAD/aging/kiran/analysis/lifespan/subclass/all_models/",brain_bank,"_",model_name,".RDATA"))

files=list.files("/sc/arion/projects/psychAD/aging/kiran/analysis/lifespan/subclass/all_models/",pattern="*.RDATA",full.names=TRUE)
files=files[grep("BIC_all_models.RDATA",files,invert=TRUE)]

BIC_list=list()
for ( ii in (1:length(files))) {
load(files[ii])
df=lapply(1:length(res.dl),function(i) data.frame("BIC"=res.dl[[i]]$BIC,"celltype"=names(res.dl)[i],"genes"=names(res.dl[[i]]$BIC)))
BIC_list[[ii]]=do.call(rbind,df)
BIC_list[[ii]]$model_type=sub(".*_all_*","",gsub(".RDATA","",basename(files[ii])))
}
save(BIC_list,file="/sc/arion/projects/psychAD/aging/kiran/analysis/lifespan/subclass/all_models/BIC_all_models.RDATA")


  
