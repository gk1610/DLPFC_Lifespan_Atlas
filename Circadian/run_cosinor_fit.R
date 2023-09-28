source("/sc/arion/projects/psychAD/aging/Circadian/scripts/cosinor.R")
args = commandArgs(trailingOnly=TRUE)
celltype=args[1]

# here res_mat is residualized data and metadata saved as RDATA object
load("/sc/arion/projects/psychAD/aging/Circadian/analysis/residuals/channel_OPC_adulthood_old.RDATA")

data=res_mat
pheno=metadata_celltype

## you can use TOD from here
metadata_daynight=read.csv("/sc/arion/projects/psychAD/aging/sleep_patterns/metadata_aging_TOD_DayNight.csv")
metadata_daynight=metadata_daynight[metadata_daynight$Day_Night!=0,]

### subsetting data for only samples with TOD values

pheno_subset=pheno[match(metadata_$SubID,colnames())]
pb_subset=subset(pb, ,SubID %in% keep_subids)

#END DATA SET UP 

#ACTUAL RUNNING OF THE COSINOR CODE
out.list = mclapply(1:nrow(cell_Com.data), function(i){ 
  if(i%%1000==0) print(i)
  res = one_cosinor_OLS(tod = as.numeric(pheno_cell_Com$TOD), y = as.numeric(cell_Com.data[i,]))
  res.onerow = as.data.frame(list(offset = res$M$est, offset.ll = res$M$CI[1], offset.ul = res$M$CI[2], 
                                  A = res$A$est, A.sd = res$A$sd, 
                                  phase = res$phase$est, phase.sd = res$phase$sd, 
                                  pvalue = res$test$pval, R2 = res$test$R2))
  return(res.onerow)
}, mc.cores = 1) 

out.tab = do.call(rbind.data.frame, out.list)
out.tab$qvalue = p.adjust(out.tab$pvalue, "BH")
out.tab$peak = 24 - out.tab$phase*24/(2*pi)
row.names(out.tab) = row.names(cell.data)
setwd("~/Desktop/Panos Collab/v0_rhythmicity/Astro//")
write.csv(out.tab,paste0("observed_Astro_Com.csv"))
