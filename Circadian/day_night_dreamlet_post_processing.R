.libPaths(c("/sc/arion/projects/CommonMind/kiran/RLib_4_3",.libPaths()))
suppressPackageStartupMessages({
library(zellkonverter)
library(basilisk)
library(dreamlet)
library(crumblr)
library(foreach)
library(doParallel)
library(data.table)
library(SingleCellExperiment)
library(dreamlet)
library(zenith)
library(DelayedArray)
library(GSEABase)
library(scater)
library(tidyverse)
library(ggplot2)
library(cowplot) 
library(org.Hs.eg.db)
library(R.utils)
})


files=list.files("/sc/arion/projects/psychAD/aging/sleep_patterns/analysis/dreamlet/",pattern="*.RDATA",full.names=TRUE)
models_list=c("model1","model2","model3")
age_group_list=c("adulthood_old_model","old_model","late_adulthood_model","early_adulthood","adulthood_model")

coef_names=c("Diff","Diff","cat_Day_NightNight:scale_age")
diff_genes=list()

ct = 1

for (ii in (1:length(models_list))){

files_subset1=files[grep(models_list[ii],files)]

for (iii in (1:length(age_group_list))){

if (age_group_list[iii]=="adulthood_model") {
files_subset2=files_subset1[grep(age_group_list[iii],files_subset1)]
files_subset2=files_subset2[grep("late|early",files_subset2,invert=TRUE)]

} else if (age_group_list[iii]=="old_model") {
files_subset2=files_subset1[grep(age_group_list[iii],files_subset1)]
files_subset2=files_subset2[grep("adulthood_old",files_subset2,invert=TRUE)]

} else {
files_subset2=files_subset1[grep(age_group_list[iii],files_subset1)]
}


for (i in (1:length(files_subset2))){

load(files_subset2[i])
df <- topTable(res.dl_model, coef = coef_names[ii],number=Inf)

diff_genes[[ct]]=data.frame("celltype"=unique(df$assay),
	"no_diff_genes"= uniqueN(df$ID[df$adj.P.Val<.05]),
	"up_diff_genes"= uniqueN(df$ID[df$adj.P.Val<.05 & df$logFC>0]),
	"down_diff_genes"= uniqueN(df$ID[df$adj.P.Val<.05 & df$logFC<0]),
	"age_group"= age_group_list[iii],
	"model"=models_list[ii],
	"ngenes"=nrow(df))

ct = ct  + 1

}

}

}


diff_genes1=do.call(rbind,diff_genes)


### example of ggplot to see model2 results

g1=ggplot(subset(diff_genes1,age_group=="adulthood_old_model" & model=="model2"),aes(factor(celltype),no_diff_genes,fill=factor(model)))+geom_bar(stat="identity")
g1+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


### example of pathway analysis

go.gs <- get_GeneOntology("CC", to = "SYMBOL")
load("/sc/arion/projects/psychAD/aging/sleep_patterns/analysis/dreamlet//subclass_day_night_EN_L3_5_IT_1_adulthood_old_model2.RDATA")
res_zenith <- zenith_gsa(res.dl_model, coef = "Diff", go.gs)
head(res_zenith,20)
plotZenithResults(res_zenith, 15, 1)
