source("/sc/arion/projects/psychAD/aging/Circadian/scripts/cosinor.R")
args = commandArgs(trailingOnly=TRUE)
celltype="OPC"

# here res_mat is residualized data and metadata saved as RDATA object
load(paste0("/sc/arion/projects/psychAD/aging/Circadian/analysis/residuals/channel_",celltype,"_adulthood_old.RDATA"))

data=res_mat
pheno=metadata_celltype

## you can use TOD from here
metadata_daynight=read.csv("/sc/arion/projects/psychAD/aging/Circadian/metadata_aging_TOD_DayNight.csv")
metadata_daynight=metadata_daynight[metadata_daynight$Day_Night!=0,]

### subset metadata and res_mat for samples with TOD

pheno_TOD=merge(metadata_daynight,pheno,by="SubID")
res_mat_with_TOD=res_mat[,match(pheno_TOD$Channel,colnames(res_mat))]

## check if all samples and colnames of res_mat_with_TOD matches
identical(as.character(pheno_TOD$Channel),colnames(res_mat_with_TOD))

#ACTUAL RUNNING OF THE COSINOR CODE
out.list = mclapply(1:nrow(res_mat_with_TOD), function(i){ 
  if(i%%1000==0) print(i)
  res = one_cosinor_OLS(tod = as.numeric(pheno_TOD$TOD), y = as.numeric(res_mat_with_TOD[i,]))
  res.onerow = as.data.frame(list(offset = res$M$est, offset.ll = res$M$CI[1], offset.ul = res$M$CI[2], 
                                  A = res$A$est, A.sd = res$A$sd, 
                                  phase = res$phase$est, phase.sd = res$phase$sd, 
                                  pvalue = res$test$pval, R2 = res$test$R2,gene=rownames(res_mat_with_TOD)[i]))
  return(res.onerow)
}, mc.cores = 1) 

out.tab = do.call(rbind.data.frame, out.list)
out.tab$qvalue = p.adjust(out.tab$pvalue, "BH")
out.tab$peak = 24 - out.tab$phase*24/(2*pi)
row.names(out.tab) = row.names(cell.data)

save(out.tab,file=paste0("/sc/arion/projects/psychAD/aging/Circadian/analysis/residuals/channel_",celltype,"_adulthood_old.RDATA"))

