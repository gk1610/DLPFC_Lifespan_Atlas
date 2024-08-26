
## this code created the data frame which has nl coef2 and coef1 from dreamlet usign the best model (poly(log2(Age + 1), df = 2))

data_dir="lifespan_trends/best_model/" 
dir.create(data_dir)

setwd("lifespan_trends/all_models/")
model_name="poly_with_two_df_log_Age" # optimal model using BIC 

## get dreamlet files from all tested models 

files=list.files(pattern=glob2rx("*nl.RDATA"),full.names=TRUE)
files=files[grep(model_name,files)]

lifespan_dream_results_list=list()
lifespan_dream_coefs_list=list()

for (i in (1:length(files))){

load(files[i])
lifespan_dream=as.data.frame(topTable(res.dl,coef=grep("poly",coefNames(res.dl),value=TRUE),num=Inf))
colnames(lifespan_dream)[grep("poly",colnames(lifespan_dream))]=paste0("coef_",1:2)
lifespan_dream$celltype=lifespan_dream$assay
lifespan_dream$celltype_ID=paste0(lifespan_dream$assay,"_",lifespan_dream$ID)
lifespan_dream_all=lifespan_dream

### separate linear and non-linear coefs and their p-values

lifespan_dream1=as.data.frame(topTable(res.dl,coef="poly(log2(Age + 1), df = 2)1",num=Inf))
colnames(lifespan_dream1)[grep("P.Value",colnames(lifespan_dream1))]=paste0("P.Value_1")
colnames(lifespan_dream1)[grep("adj.P.Val",colnames(lifespan_dream1))]=paste0("adj.P.Val_1")
lifespan_dream1$celltype_ID=paste0(lifespan_dream1$assay,"_",lifespan_dream1$ID)


lifespan_dream2=as.data.frame(topTable(res.dl,coef="poly(log2(Age + 1), df = 2)2",num=Inf))
colnames(lifespan_dream2)[grep("P.Value",colnames(lifespan_dream2))]=paste0("P.Value_2")
colnames(lifespan_dream2)[grep("adj.P.Val",colnames(lifespan_dream2))]=paste0("adj.P.Val_2")
lifespan_dream2$celltype_ID=paste0(lifespan_dream2$assay,"_",lifespan_dream2$ID)

lifespan_dream=merge(merge(lifespan_dream_all,lifespan_dream1[,c("celltype_ID","adj.P.Val_1","P.Value_1")],by="celltype_ID"),lifespan_dream2[,c("celltype_ID","adj.P.Val_2","P.Value_2")],by="celltype_ID")

for (ii in (1:length(res.dl))){
lifespan_dream_coefs[[ii]]=as.data.frame(res.dl[[ii]]$coefficients[,grep("poly|Intercept",colnames(res.dl[[ii]]$coefficients))])
lifespan_dream_coefs[[ii]]$ID=rownames(res.dl[[ii]]$coefficients)
lifespan_dream_coefs[[ii]]$assay=names(res.dl)[ii]

}

lifespan_dream_coefs1=do.call(rbind,lifespan_dream_coefs)
lifespan_dream_coefs=lifespan_dream_coefs1
colnames(lifespan_dream_coefs)[1:3]=c("Intercept",paste0("coef_",1:2))
lifespan_dream_coefs$celltype=lifespan_dream_coefs$assay
lifespan_dream_coefs$celltype_ID=paste0(lifespan_dream_coefs$celltype,"_",lifespan_dream_coefs$ID)
lifespan_dream_coefs_list[[i]]=lifespan_dream_coefs
lifespan_dream_results_list[[i]]=lifespan_dream

}

lifespan_dream_coefs_df=do.call(rbind,lifespan_dream_coefs_list)
lifespan_dream_coefs_df=as.data.frame(lifespan_dream_coefs_df)

lifespan_dream_results_df=do.call(rbind,lifespan_dream_results_list)
lifespan_dream_results_df=as.data.frame(lifespan_dream_results_df)
rownames(lifespan_dream_coefs_df)=lifespan_dream_coefs_df$celltype_ID

## we save this file as .feather for clustering in python as it is faster

library(feather)
X = model.matrix(~ poly(coef_2, df=1) + poly(coef_1, df=1),data=lifespan_dream_coefs_df)
X = as.data.frame(X)
X$celltype_ID=rownames(X)
feather::write_feather(X,paste0(data_dir,"/best_model_coefs.feather"))

save(lifespan_dream_coefs_df,lifespan_dream_results_df,file=paste0(data_dir,"/lifespan_coefs_dream_results_best_model.RDATA"))

